#include <osmium/io/reader.hpp>
#include <osmium/io/pbf_input.hpp>
#include <osmium/handler.hpp>
#include <osmium/osm.hpp>
#include <osmium/visitor.hpp>
#include <osmium/util/progress_bar.hpp>
#include <osmium/memory/buffer.hpp>

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <iostream>
#include <string>
#include <filesystem>
#include <limits>
#include <algorithm>

#include "gdal.h"
#include "ogrsf_frmts.h"

// =========================
// Constantes
// =========================
static constexpr uint64_t DIST_INF = std::numeric_limits<uint64_t>::max();
static constexpr uint32_t NO_NODE  = std::numeric_limits<uint32_t>::max();

// =========================
// Coordonnées E7
// =========================
struct CoordE7 {
    int32_t lat_e7 = 0;
    int32_t lon_e7 = 0;
};

static double e7_to_deg(int32_t x) {
    return static_cast<double>(x) / 1e7;
}

// =========================
// Arc temporaire (ordre de génération)
// =========================
struct ArcTmp {
    uint32_t u = 0;        // graph node id
    uint32_t v = 0;        // graph node id
    uint32_t len_cm = 0;
    uint8_t  hw_class = 0;
    int64_t  way_id = 0;   // OSM way id
};

// =========================
// Restriction relation OSM (via=node uniquement)
// =========================
struct RestrRel {
    int64_t from_way = 0;
    int64_t to_way = 0;
    int64_t via_node_osm = 0; // OSM node id
    std::string kind;         // restriction / restriction:motorcar
};

// =========================
// Haversine
// =========================
static double haversine(double lat1, double lon1, double lat2, double lon2) {
    constexpr double R = 6371000.0;
    constexpr double DEG = M_PI / 180.0;

    double dlat = (lat2 - lat1) * DEG;
    double dlon = (lon2 - lon1) * DEG;
    lat1 *= DEG;
    lat2 *= DEG;

    double a = std::sin(dlat/2)*std::sin(dlat/2)
             + std::sin(dlon/2)*std::sin(dlon/2)*std::cos(lat1)*std::cos(lat2);
    return 2 * R * std::atan2(std::sqrt(a), std::sqrt(1-a));
}

// =========================
// Apply avec progress
// =========================
template <typename Handler>
void apply_with_progress(const std::string& filename, const std::string& label, Handler& handler) {
    osmium::io::Reader reader(filename);
    const uint64_t total = reader.file_size();

    std::cerr << label << "\n";
    osmium::ProgressBar bar(total, true);

    while (osmium::memory::Buffer buffer = reader.read()) {
        osmium::apply(buffer, handler);
        bar.update(reader.offset());
    }

    bar.done();
    reader.close();
}

// =========================
// PASS 1 : compter usage nodes highway + endpoints
// =========================
class NodeUseCounter : public osmium::handler::Handler {
public:
    std::unordered_map<osmium::object_id_type, uint8_t> use; // saturé à 2
    std::unordered_set<osmium::object_id_type> endpoints;

    void way(const osmium::Way& w) {
        if (!w.tags().has_key("highway")) return;
        const auto& nodes = w.nodes();
        if (nodes.size() < 2) return;

        endpoints.insert(nodes.front().ref());
        endpoints.insert(nodes.back().ref());

        for (const auto& n : nodes) {
            uint8_t& c = use[n.ref()];
            if (c < 2) ++c;
        }
    }
};

// =========================
// PASS 2 : lire relations restriction (via=node)
// =========================
class RestrictionCollector : public osmium::handler::Handler {
public:
    std::vector<RestrRel>& restrictions;
    std::unordered_set<osmium::object_id_type>& via_nodes; // OSM node ids to force as graph nodes

    uint64_t ignored_no_tag = 0;
    uint64_t ignored_bad_members = 0;
    uint64_t ignored_via_way = 0;

    RestrictionCollector(std::vector<RestrRel>& r, std::unordered_set<osmium::object_id_type>& via)
        : restrictions(r), via_nodes(via) {}

    void relation(const osmium::Relation& rel) {
        const char* type = rel.tags()["type"];
        if (!type || std::string(type) != "restriction") return;

        const char* r_motorcar = rel.tags()["restriction:motorcar"];
        const char* r_value = r_motorcar ? r_motorcar : rel.tags()["restriction"];
        if (!r_value) { ++ignored_no_tag; return; }

        int64_t from_way = 0;
        int64_t to_way = 0;
        int64_t via_node = 0;
        bool via_is_way = false;

        for (const auto& m : rel.members()) {
            const std::string role = m.role();
            if (role == "from" && m.type() == osmium::item_type::way) {
                if (from_way == 0) from_way = m.ref();
            } else if (role == "to" && m.type() == osmium::item_type::way) {
                if (to_way == 0) to_way = m.ref();
            } else if (role == "via") {
                if (m.type() == osmium::item_type::node) {
                    if (via_node == 0) via_node = m.ref();
                } else if (m.type() == osmium::item_type::way) {
                    via_is_way = true;
                }
            }
        }

        if (via_is_way) { ++ignored_via_way; return; }
        if (from_way == 0 || to_way == 0 || via_node == 0) { ++ignored_bad_members; return; }

        restrictions.push_back({from_way, to_way, via_node, std::string(r_value)});
        via_nodes.insert(static_cast<osmium::object_id_type>(via_node));
    }
};

// =========================
// PASS 3 : collecter coords + node_index des nodes de graphe
//  - graph node si (use>=2) ou endpoint ou via restriction
// =========================
class NodeCollector : public osmium::handler::Handler {
public:
    const std::unordered_map<osmium::object_id_type, uint8_t>& use;
    const std::unordered_set<osmium::object_id_type>& endpoints;
    const std::unordered_set<osmium::object_id_type>& via_nodes;

    std::unordered_map<osmium::object_id_type, CoordE7>& coords_all;         // OSM node id -> coord
    std::unordered_map<osmium::object_id_type, uint32_t>& node_index;        // OSM node id -> graph node id

    NodeCollector(
        const std::unordered_map<osmium::object_id_type, uint8_t>& u,
        const std::unordered_set<osmium::object_id_type>& e,
        const std::unordered_set<osmium::object_id_type>& v,
        std::unordered_map<osmium::object_id_type, CoordE7>& ca,
        std::unordered_map<osmium::object_id_type, uint32_t>& ni
    ) : use(u), endpoints(e), via_nodes(v), coords_all(ca), node_index(ni) {}

    void node(const osmium::Node& n) {
        if (!n.location().valid()) return;
        const auto id = n.id();

        // on ne stocke les coords que si le node est utile (apparait dans un way highway OU est un via)
        if (use.find(id) == use.end() && via_nodes.find(id) == via_nodes.end()) return;

        coords_all[id] = CoordE7{
            static_cast<int32_t>(std::llround(n.location().lat() * 1e7)),
            static_cast<int32_t>(std::llround(n.location().lon() * 1e7))
        };

        bool is_graph_node = false;

        auto it = use.find(id);
        if (it != use.end() && it->second >= 2) is_graph_node = true;
        if (endpoints.find(id) != endpoints.end()) is_graph_node = true;
        if (via_nodes.find(id) != via_nodes.end()) is_graph_node = true;

        if (is_graph_node && node_index.find(id) == node_index.end()) {
            node_index[id] = static_cast<uint32_t>(node_index.size());
        }
    }
};

// =========================
// Utilitaires oneway
// =========================
static bool tag_is_oneway_yes(const char* v) {
    if (!v) return false;
    const std::string s(v);
    return s == "yes" || s == "true" || s == "1";
}
static bool tag_is_oneway_minus1(const char* v) {
    if (!v) return false;
    return std::string(v) == "-1";
}

// =========================
// PASS 4 : découper ways, sens, arcs, et stocker géométrie par arc (ordre génération)
// =========================
struct EdgeGeom {
    std::vector<std::pair<double,double>> lonlat; // [lon,lat]
    std::string highway;
    std::string name;
    std::string ref;
    std::string maxspeed;
    std::string oneway;
    std::string junction;
};

class WayToArcs : public osmium::handler::Handler {
public:
    const std::unordered_map<osmium::object_id_type, uint32_t>& node_index;
    const std::unordered_map<osmium::object_id_type, CoordE7>& coords_all;

    std::vector<ArcTmp>& arcs;
    std::vector<EdgeGeom>& arc_geom; // même ordre que arcs

    WayToArcs(
        const std::unordered_map<osmium::object_id_type, uint32_t>& ni,
        const std::unordered_map<osmium::object_id_type, CoordE7>& ca,
        std::vector<ArcTmp>& a,
        std::vector<EdgeGeom>& g
    ) : node_index(ni), coords_all(ca), arcs(a), arc_geom(g) {}

    void way(const osmium::Way& w) {
        if (!w.tags().has_key("highway")) return;

        const auto& ns = w.nodes();
        if (ns.size() < 2) return;

        bool forward = true;
        bool backward = true;

        const char* oneway = w.tags()["oneway"];
        const char* junction = w.tags()["junction"];

        if (junction && std::string(junction) == "roundabout") {
            forward = true;
            backward = false;
        } else if (tag_is_oneway_yes(oneway)) {
            forward = true;
            backward = false;
        } else if (tag_is_oneway_minus1(oneway)) {
            forward = false;
            backward = true;
        } else {
            forward = true;
            backward = true;
        }

        auto get_lonlat = [&](osmium::object_id_type id, std::pair<double,double>& out) -> bool {
            auto it = coords_all.find(id);
            if (it == coords_all.end()) return false;
            const double lat = e7_to_deg(it->second.lat_e7);
            const double lon = e7_to_deg(it->second.lon_e7);
            out = {lon, lat};
            return true;
        };

        osmium::object_id_type start_osm = 0;
        bool has_start = false;
        double acc_m = 0.0;
        std::vector<std::pair<double,double>> segment;

        const auto id0 = ns.front().ref();
        if (node_index.find(id0) != node_index.end()) {
            start_osm = id0;
            has_start = true;
            segment.clear();
            std::pair<double,double> p0;
            if (!get_lonlat(id0, p0)) has_start = false;
            else segment.push_back(p0);
        }

        for (size_t i = 1; i < ns.size(); ++i) {
            const auto a_id = ns[i-1].ref();
            const auto b_id = ns[i].ref();

            auto ca_it = coords_all.find(a_id);
            auto cb_it = coords_all.find(b_id);
            if (ca_it == coords_all.end() || cb_it == coords_all.end()) {
                acc_m = 0.0;
                has_start = false;
                segment.clear();
                if (node_index.find(b_id) != node_index.end()) {
                    start_osm = b_id;
                    has_start = true;
                    std::pair<double,double> pb;
                    if (get_lonlat(b_id, pb)) segment.push_back(pb);
                    else { has_start = false; segment.clear(); }
                }
                continue;
            }

            const auto& A = ca_it->second;
            const auto& B = cb_it->second;
            acc_m += haversine(
                e7_to_deg(A.lat_e7), e7_to_deg(A.lon_e7),
                e7_to_deg(B.lat_e7), e7_to_deg(B.lon_e7)
            );

            if (has_start) {
                std::pair<double,double> pb;
                if (!get_lonlat(b_id, pb)) {
                    has_start = false;
                    acc_m = 0.0;
                    segment.clear();
                } else {
                    segment.push_back(pb);
                }
            } else {
                if (node_index.find(b_id) != node_index.end()) {
                    start_osm = b_id;
                    has_start = true;
                    acc_m = 0.0;
                    segment.clear();
                    std::pair<double,double> pb;
                    if (get_lonlat(b_id, pb)) segment.push_back(pb);
                    else { has_start = false; segment.clear(); }
                }
                continue;
            }

            if (node_index.find(b_id) != node_index.end()) {
                if (has_start && start_osm != b_id && segment.size() >= 2) {
                    auto su_it = node_index.find(start_osm);
                    auto sv_it = node_index.find(b_id);
                    if (su_it != node_index.end() && sv_it != node_index.end()) {
                        const uint32_t u = su_it->second;
                        const uint32_t v = sv_it->second;
                        const uint32_t len_cm = static_cast<uint32_t>(std::llround(acc_m * 100.0));

                        auto fill_geom = [&](EdgeGeom& eg, const char* dir) {
                            eg.lonlat = segment;
                            if (dir && std::string(dir) == "backward") {
                                std::reverse(eg.lonlat.begin(), eg.lonlat.end());
                            }
                            const char* highway  = w.tags()["highway"];
                            const char* name     = w.tags()["name"];
                            const char* ref      = w.tags()["ref"];
                            const char* maxspeed = w.tags()["maxspeed"];
                            const char* oneway_t = w.tags()["oneway"];
                            const char* junction_t = w.tags()["junction"];
                            eg.highway  = highway ? highway : "";
                            eg.name     = name ? name : "";
                            eg.ref      = ref ? ref : "";
                            eg.maxspeed = maxspeed ? maxspeed : "";
                            eg.oneway   = oneway_t ? oneway_t : "";
                            eg.junction = junction_t ? junction_t : "";
                        };

                        if (forward) {
                            arcs.push_back({u, v, len_cm, 0, w.id()});
                            arc_geom.emplace_back();
                            fill_geom(arc_geom.back(), "forward");
                        }
                        if (backward) {
                            arcs.push_back({v, u, len_cm, 0, w.id()});
                            arc_geom.emplace_back();
                            fill_geom(arc_geom.back(), "backward");
                        }
                    }
                }

                start_osm = b_id;
                has_start = true;
                acc_m = 0.0;
                segment.clear();
                std::pair<double,double> pb;
                if (get_lonlat(b_id, pb)) segment.push_back(pb);
                else { has_start = false; segment.clear(); }
            }
        }
    }
};

// =========================
// CSR builder + mapping (ID d'arc = index CSR)
//   csr_to_tmp[p] = index dans arcs/arc_geom (ordre génération)
// =========================
static void build_csr(
    uint32_t N,
    const std::vector<ArcTmp>& arcs,
    std::vector<uint32_t>& offsets,
    std::vector<uint32_t>& to,
    std::vector<uint32_t>& len_cm,
    std::vector<uint8_t>& cls,
    std::vector<int64_t>& wid,
    std::vector<uint32_t>& csr_to_tmp
) {
    const uint32_t M = static_cast<uint32_t>(arcs.size());
    std::vector<uint32_t> deg(N, 0);
    for (const auto& a : arcs) ++deg[a.u];

    offsets.assign(N + 1, 0);
    for (uint32_t i = 0; i < N; ++i) offsets[i + 1] = offsets[i] + deg[i];

    to.assign(M, 0);
    len_cm.assign(M, 0);
    cls.assign(M, 0);
    wid.assign(M, 0);
    csr_to_tmp.assign(M, 0);

    auto cur = offsets;
    for (uint32_t tmp_i = 0; tmp_i < M; ++tmp_i) {
        const auto& a = arcs[tmp_i];
        const uint32_t p = cur[a.u]++; // arc id = p
        to[p] = a.v;
        len_cm[p] = a.len_cm;
        cls[p] = a.hw_class;
        wid[p] = a.way_id;
        csr_to_tmp[p] = tmp_i;
    }
}

// incoming arcs index (rev) : pour chaque node v, liste des arcs entrants (indices CSR)
static void build_reverse_arc_index(
    uint32_t N,
    const std::vector<uint32_t>& offsets,
    const std::vector<uint32_t>& to,
    std::vector<uint32_t>& rev_off,
    std::vector<uint32_t>& rev_arc
) {
    const uint32_t M = static_cast<uint32_t>(to.size());
    std::vector<uint32_t> indeg(N, 0);

    for (uint32_t u = 0; u < N; ++u) {
        for (uint32_t p = offsets[u]; p < offsets[u + 1]; ++p) {
            ++indeg[to[p]];
        }
    }

    rev_off.assign(N + 1, 0);
    for (uint32_t i = 0; i < N; ++i) rev_off[i + 1] = rev_off[i] + indeg[i];

    rev_arc.assign(M, 0);
    auto cur = rev_off;
    for (uint32_t u = 0; u < N; ++u) {
        for (uint32_t p = offsets[u]; p < offsets[u + 1]; ++p) {
            const uint32_t v = to[p];
            const uint32_t pos = cur[v]++;
            rev_arc[pos] = p;
        }
    }
}

static bool is_only_restriction(const std::string& k) { return k.rfind("only_", 0) == 0; }
static bool is_no_restriction(const std::string& k) { return k.rfind("no_", 0) == 0; }

// transitions interdites (arc_in, arc_out) + index vers restriction source
struct ForbiddenTurn {
    uint32_t in_arc = 0;     // CSR arc id
    uint32_t out_arc = 0;    // CSR arc id
    uint32_t restr_idx = 0;  // index dans rels
    uint32_t via = 0;        // graph node id
};

// Construction transitions interdites : IDs = indices CSR
static void build_forbidden_turns(
    const std::vector<RestrRel>& rels,
    const std::unordered_map<osmium::object_id_type, uint32_t>& node_index,
    const std::vector<uint32_t>& offsets,
    const std::vector<int64_t>& wid,
    const std::vector<uint32_t>& rev_off,
    const std::vector<uint32_t>& rev_arc,
    std::vector<ForbiddenTurn>& forbidden
) {
    forbidden.clear();
    forbidden.reserve(rels.size() * 8);

    std::unordered_set<uint64_t> seen;
    seen.reserve(rels.size() * 16);

    for (uint32_t ri = 0; ri < static_cast<uint32_t>(rels.size()); ++ri) {
        const auto& rr = rels[ri];

        auto it_v = node_index.find(static_cast<osmium::object_id_type>(rr.via_node_osm));
        if (it_v == node_index.end()) continue;
        const uint32_t via = it_v->second;

        std::vector<uint32_t> in_arcs;   // CSR arc ids
        std::vector<uint32_t> out_to;    // CSR arc ids
        std::vector<uint32_t> out_all;   // CSR arc ids

        for (uint32_t p = offsets[via]; p < offsets[via + 1]; ++p) {
            out_all.push_back(p);
            if (wid[p] == rr.to_way) out_to.push_back(p);
        }

        for (uint32_t q = rev_off[via]; q < rev_off[via + 1]; ++q) {
            const uint32_t p = rev_arc[q];
            if (wid[p] == rr.from_way) in_arcs.push_back(p);
        }

        if (in_arcs.empty() || out_all.empty()) continue;

        if (is_no_restriction(rr.kind)) {
            if (out_to.empty()) continue;
            for (uint32_t a : in_arcs) {
                for (uint32_t b : out_to) {
                    const uint64_t key = (static_cast<uint64_t>(a) << 32) | b;
                    if (seen.insert(key).second) forbidden.push_back({a, b, ri, via});
                }
            }
        } else if (is_only_restriction(rr.kind)) {
            if (out_to.empty()) continue;
            std::unordered_set<uint32_t> allowed(out_to.begin(), out_to.end());
            for (uint32_t a : in_arcs) {
                for (uint32_t b : out_all) {
                    if (allowed.find(b) != allowed.end()) continue;
                    const uint64_t key = (static_cast<uint64_t>(a) << 32) | b;
                    if (seen.insert(key).second) forbidden.push_back({a, b, ri, via});
                }
            }
        }
    }
}

// =========================
// Graph file v2 (avec restrictions) - IDs cohérents
// Format little-endian:
//  char[4] "GRF2"
//  uint32 N
//  uint32 M
//  offsets uint32*(N+1)
//  to      uint32*M
//  len_cm  uint32*M
//  cls     uint8*M
//  way_id  int64*M
//  dist    uint64*N
//  label   uint32*N
//  pred    uint32*N
//  uint32 R
//  pairs   (uint32 in_arc, uint32 out_arc) * R
// =========================
static void write_graph_v2_with_restrictions(
    const std::string& path,
    uint32_t N,
    const std::vector<uint32_t>& offsets,
    const std::vector<uint32_t>& to,
    const std::vector<uint32_t>& len_cm,
    const std::vector<uint8_t>& cls,
    const std::vector<int64_t>& wid,
    const std::vector<ForbiddenTurn>& forbidden
) {
    const uint32_t M = static_cast<uint32_t>(to.size());
    const uint32_t R = static_cast<uint32_t>(forbidden.size());

    std::vector<uint64_t> dist(N, DIST_INF);
    std::vector<uint32_t> label(N, NO_NODE);
    std::vector<uint32_t> pred(N, NO_NODE);

    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("cannot open output graph file");

    const char magic[4] = {'G','R','F','2'};
    out.write(magic, 4);

    out.write((char*)&N, 4);
    out.write((char*)&M, 4);

    out.write((char*)offsets.data(), offsets.size() * 4);
    out.write((char*)to.data(),      to.size() * 4);
    out.write((char*)len_cm.data(),  len_cm.size() * 4);
    out.write((char*)cls.data(),     cls.size());
    out.write((char*)wid.data(),     wid.size() * 8);

    out.write((char*)dist.data(),    dist.size() * 8);
    out.write((char*)label.data(),   label.size() * 4);
    out.write((char*)pred.data(),    pred.size() * 4);

    out.write((char*)&R, 4);
    for (const auto& ft : forbidden) {
        out.write((char*)&ft.in_arc, 4);
        out.write((char*)&ft.out_arc, 4);
    }
}

// =========================
// FlatGeobuf exports (nodes, edges, restrictions) - IDs cohérents
// =========================
static std::string base_name_no_ext(const std::string& graph_path) {
    std::filesystem::path p(graph_path);
    auto stem = p.stem().string();
    auto parent = p.parent_path();
    if (parent.empty()) return stem;
    return (parent / stem).string();
}

static void build_graph_node_arrays(
    uint32_t N,
    const std::unordered_map<osmium::object_id_type, uint32_t>& node_index,
    const std::unordered_map<osmium::object_id_type, CoordE7>& coords_all,
    std::vector<int64_t>& osm_node_id_by_idx,
    std::vector<double>& lon_by_idx,
    std::vector<double>& lat_by_idx
) {
    osm_node_id_by_idx.assign(N, -1);
    lon_by_idx.assign(N, 0.0);
    lat_by_idx.assign(N, 0.0);

    for (const auto& kv : node_index) {
        const osmium::object_id_type osm_id = kv.first;
        const uint32_t idx = kv.second;

        auto itc = coords_all.find(osm_id);
        if (itc == coords_all.end()) continue;

        osm_node_id_by_idx[idx] = static_cast<int64_t>(osm_id);
        lat_by_idx[idx] = e7_to_deg(itc->second.lat_e7);
        lon_by_idx[idx] = e7_to_deg(itc->second.lon_e7);
    }
}

static void export_nodes_fgb(
    const std::string& path,
    const std::vector<int64_t>& osm_node_id_by_idx,
    const std::vector<double>& lon_by_idx,
    const std::vector<double>& lat_by_idx
) {
    GDALAllRegister();
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("FlatGeobuf");
    if (!drv) throw std::runtime_error("GDAL FlatGeobuf driver not available");

    GDALDataset* ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (!ds) throw std::runtime_error("Cannot create FlatGeobuf file for nodes");

    OGRSpatialReference srs;
    srs.SetWellKnownGeogCS("WGS84");

    OGRLayer* layer = ds->CreateLayer("node", &srs, wkbPoint, nullptr);
    if (!layer) { GDALClose(ds); throw std::runtime_error("Cannot create nodes layer"); }

    layer->CreateField(new OGRFieldDefn("id", OFTInteger));            // graph node id
    layer->CreateField(new OGRFieldDefn("osm_node_id", OFTInteger64)); // OSM id

    OGRFeatureDefn* defn = layer->GetLayerDefn();

    for (int id = 0; id < static_cast<int>(osm_node_id_by_idx.size()); ++id) {
        if (osm_node_id_by_idx[id] < 0) continue;

        OGRFeature* feat = OGRFeature::CreateFeature(defn);
        feat->SetField("id", id);
        feat->SetField("osm_node_id", static_cast<GIntBig>(osm_node_id_by_idx[id]));

        OGRPoint pt(lon_by_idx[id], lat_by_idx[id]);
        feat->SetGeometry(&pt);

        if (layer->CreateFeature(feat) != OGRERR_NONE) {
            OGRFeature::DestroyFeature(feat);
            GDALClose(ds);
            throw std::runtime_error("Failed writing node feature");
        }
        OGRFeature::DestroyFeature(feat);
    }

    GDALClose(ds);
}

// Edges export : ID d'arc = index CSR (p)
static void export_edges_fgb(
    const std::string& path,
    uint32_t N,
    const std::vector<uint32_t>& offsets,
    const std::vector<uint32_t>& to,
    const std::vector<uint32_t>& len_cm,
    const std::vector<int64_t>& wid,
    const std::vector<uint32_t>& csr_to_tmp,
    const std::vector<EdgeGeom>& geom_by_tmp
) {
    GDALAllRegister();
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("FlatGeobuf");
    if (!drv) throw std::runtime_error("GDAL FlatGeobuf driver not available");

    GDALDataset* ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (!ds) throw std::runtime_error("Cannot create FlatGeobuf file for edges");

    OGRSpatialReference srs;
    srs.SetWellKnownGeogCS("WGS84");

    OGRLayer* layer = ds->CreateLayer("edge", &srs, wkbLineString, nullptr);
    if (!layer) { GDALClose(ds); throw std::runtime_error("Cannot create edges layer"); }

    layer->CreateField(new OGRFieldDefn("id", OFTInteger));      // arc id (CSR index)
    layer->CreateField(new OGRFieldDefn("u", OFTInteger));       // graph node id
    layer->CreateField(new OGRFieldDefn("v", OFTInteger));       // graph node id
    layer->CreateField(new OGRFieldDefn("len_cm", OFTInteger));
    layer->CreateField(new OGRFieldDefn("way_id", OFTInteger64));

    layer->CreateField(new OGRFieldDefn("highway", OFTString));
    layer->CreateField(new OGRFieldDefn("name", OFTString));
    layer->CreateField(new OGRFieldDefn("ref", OFTString));
    layer->CreateField(new OGRFieldDefn("maxspeed", OFTString));
    layer->CreateField(new OGRFieldDefn("oneway", OFTString));
    layer->CreateField(new OGRFieldDefn("junction", OFTString));

    OGRFeatureDefn* defn = layer->GetLayerDefn();
    const uint32_t M = static_cast<uint32_t>(to.size());

    for (uint32_t u = 0; u < N; ++u) {
        for (uint32_t p = offsets[u]; p < offsets[u + 1]; ++p) {
            const uint32_t v = to[p];
            const uint32_t tmp_i = csr_to_tmp[p];
            const auto& g = geom_by_tmp[tmp_i];
            if (g.lonlat.size() < 2) continue;

            OGRFeature* feat = OGRFeature::CreateFeature(defn);

            feat->SetField("id", static_cast<int>(p));
            feat->SetField("u", static_cast<int>(u));
            feat->SetField("v", static_cast<int>(v));
            feat->SetField("len_cm", static_cast<int>(len_cm[p]));
            feat->SetField("way_id", static_cast<GIntBig>(wid[p]));

            feat->SetField("highway", g.highway.c_str());
            feat->SetField("name", g.name.c_str());
            feat->SetField("ref", g.ref.c_str());
            feat->SetField("maxspeed", g.maxspeed.c_str());
            feat->SetField("oneway", g.oneway.c_str());
            feat->SetField("junction", g.junction.c_str());

            OGRLineString ls;
            ls.setNumPoints(static_cast<int>(g.lonlat.size()));
            for (int k = 0; k < static_cast<int>(g.lonlat.size()); ++k) {
                ls.setPoint(k, g.lonlat[k].first, g.lonlat[k].second);
            }
            feat->SetGeometry(&ls);

            if (layer->CreateFeature(feat) != OGRERR_NONE) {
                OGRFeature::DestroyFeature(feat);
                GDALClose(ds);
                throw std::runtime_error("Failed writing edge feature");
            }

            OGRFeature::DestroyFeature(feat);
        }
    }

    (void)M;
    GDALClose(ds);
}

// Restrictions export : via_node = graph node id ; in_arc/out_arc = CSR arc ids
static void export_restrictions_fgb(
    const std::string& path,
    const std::vector<ForbiddenTurn>& forbidden,
    const std::vector<RestrRel>& rels,
    const std::vector<int64_t>& osm_node_id_by_idx,
    const std::vector<double>& lon_by_idx,
    const std::vector<double>& lat_by_idx
) {
    GDALAllRegister();
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("FlatGeobuf");
    if (!drv) throw std::runtime_error("GDAL FlatGeobuf driver not available");

    GDALDataset* ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (!ds) throw std::runtime_error("Cannot create FlatGeobuf file for restrictions");

    OGRSpatialReference srs;
    srs.SetWellKnownGeogCS("WGS84");

    OGRLayer* layer = ds->CreateLayer("restriction", &srs, wkbPoint, nullptr);
    if (!layer) { GDALClose(ds); throw std::runtime_error("Cannot create restrictions layer"); }

    layer->CreateField(new OGRFieldDefn("id", OFTInteger));
    layer->CreateField(new OGRFieldDefn("kind", OFTString));
    layer->CreateField(new OGRFieldDefn("from_way", OFTInteger64));
    layer->CreateField(new OGRFieldDefn("to_way", OFTInteger64));
    layer->CreateField(new OGRFieldDefn("via_node", OFTInteger));        // graph node id
    layer->CreateField(new OGRFieldDefn("via_osm", OFTInteger64));        // OSM node id (debug)
    layer->CreateField(new OGRFieldDefn("in_arc", OFTInteger));          // CSR arc id
    layer->CreateField(new OGRFieldDefn("out_arc", OFTInteger));         // CSR arc id

    OGRFeatureDefn* defn = layer->GetLayerDefn();

    for (int i = 0; i < static_cast<int>(forbidden.size()); ++i) {
        const auto& ft = forbidden[i];
        if (ft.restr_idx >= rels.size()) continue;
        const auto& rr = rels[ft.restr_idx];

        const uint32_t via = ft.via;
        if (via >= osm_node_id_by_idx.size()) continue;
        if (osm_node_id_by_idx[via] < 0) continue;

        const double lon = lon_by_idx[via];
        const double lat = lat_by_idx[via];

        OGRFeature* feat = OGRFeature::CreateFeature(defn);

        feat->SetField("id", i);
        feat->SetField("kind", rr.kind.c_str());
        feat->SetField("from_way", static_cast<GIntBig>(rr.from_way));
        feat->SetField("to_way", static_cast<GIntBig>(rr.to_way));
        feat->SetField("via_node", static_cast<int>(via));
        feat->SetField("via_osm", static_cast<GIntBig>(rr.via_node_osm));
        feat->SetField("in_arc", static_cast<int>(ft.in_arc));
        feat->SetField("out_arc", static_cast<int>(ft.out_arc));

        OGRPoint pt(lon, lat);
        feat->SetGeometry(&pt);

        if (layer->CreateFeature(feat) != OGRERR_NONE) {
            OGRFeature::DestroyFeature(feat);
            GDALClose(ds);
            throw std::runtime_error("Failed writing restriction feature");
        }

        OGRFeature::DestroyFeature(feat);
    }

    GDALClose(ds);
}

// =========================
// MAIN
// usage:
//   build_graph input.osm.pbf output.graph [--fgb]
// sorties :
//   output.graph (GRF2)
//   <base>_node.fgb, <base>_edge.fgb, <base>_restriction.fgb (si --fgb)
// =========================
int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: build_graph input.osm.pbf output.graph [--fgb]\n";
        return 1;
    }

    const std::string in = argv[1];
    const std::string out_graph = argv[2];

    bool emit_fgb = false;
    for (int i = 3; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--fgb" || a == "-fgb" || a == "--gfb" || a == "-gfb") emit_fgb = true;
    }

    // PASS 1: node usage
    NodeUseCounter counter;
    apply_with_progress(in, "PASS 1: counting node usage", counter);

    // PASS 2: restrictions (collect via nodes)
    std::vector<RestrRel> rels;
    rels.reserve(200000);
    std::unordered_set<osmium::object_id_type> via_nodes;
    via_nodes.reserve(200000);

    RestrictionCollector rc(rels, via_nodes);
    apply_with_progress(in, "PASS 2: reading restrictions (via=node only)", rc);

    // PASS 3: coords + node_index (force via nodes as graph nodes)
    std::unordered_map<osmium::object_id_type, CoordE7> coords_all;
    coords_all.reserve(counter.use.size() + via_nodes.size());

    std::unordered_map<osmium::object_id_type, uint32_t> node_index;
    node_index.reserve((counter.use.size() / 4) + via_nodes.size());

    NodeCollector nc(counter.use, counter.endpoints, via_nodes, coords_all, node_index);
    apply_with_progress(in, "PASS 3: collecting coords + graph nodes", nc);

    const uint32_t N = static_cast<uint32_t>(node_index.size());

    // PASS 4: arcs
    std::vector<ArcTmp> arcs;
    std::vector<EdgeGeom> arc_geom;
    arcs.reserve(1000000);
    arc_geom.reserve(1000000);

    WayToArcs wta(node_index, coords_all, arcs, arc_geom);
    apply_with_progress(in, "PASS 4: building arcs (oneway + roundabout applied)", wta);

    // CSR + mapping arc-id
    std::vector<uint32_t> offsets;
    std::vector<uint32_t> to;
    std::vector<uint32_t> len_cm;
    std::vector<uint8_t> cls;
    std::vector<int64_t> wid;
    std::vector<uint32_t> csr_to_tmp;

    build_csr(N, arcs, offsets, to, len_cm, cls, wid, csr_to_tmp);
    const uint32_t M = static_cast<uint32_t>(to.size());

    // reverse arcs
    std::vector<uint32_t> rev_off;
    std::vector<uint32_t> rev_arc;
    build_reverse_arc_index(N, offsets, to, rev_off, rev_arc);

    // forbidden transitions
    std::vector<ForbiddenTurn> forbidden;
    build_forbidden_turns(rels, node_index, offsets, wid, rev_off, rev_arc, forbidden);

    std::cerr << "N=" << N
              << " M=" << M
              << " rels=" << rels.size()
              << " forbidden=" << forbidden.size()
              << " ignored(via_way)=" << rc.ignored_via_way
              << " ignored(no_tag)=" << rc.ignored_no_tag
              << " ignored(bad_members)=" << rc.ignored_bad_members
              << "\n";

    // graph file
    write_graph_v2_with_restrictions(out_graph, N, offsets, to, len_cm, cls, wid, forbidden);

    // exports fgb (3 fichiers)
    if (emit_fgb) {
        const std::string base = base_name_no_ext(out_graph);
        const std::string node_path  = base + "_node.fgb";
        const std::string edge_path  = base + "_edge.fgb";
        const std::string restr_path = base + "_restriction.fgb";

        std::vector<int64_t> osm_node_id_by_idx;
        std::vector<double> lon_by_idx;
        std::vector<double> lat_by_idx;

        build_graph_node_arrays(N, node_index, coords_all, osm_node_id_by_idx, lon_by_idx, lat_by_idx);

        export_nodes_fgb(node_path, osm_node_id_by_idx, lon_by_idx, lat_by_idx);
        export_edges_fgb(edge_path, N, offsets, to, len_cm, wid, csr_to_tmp, arc_geom);
        export_restrictions_fgb(restr_path, forbidden, rels, osm_node_id_by_idx, lon_by_idx, lat_by_idx);
    }

    return 0;
}
