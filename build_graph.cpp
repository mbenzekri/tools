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

static double e7_to_deg(int32_t x) { return static_cast<double>(x) / 1e7; }

// =========================
// Arc temporaire
// =========================
struct ArcTmp {
    uint32_t u = 0;
    uint32_t v = 0;
    uint32_t len_cm = 0;
    uint8_t  hw_class = 0;
    int64_t  way_id = 0;
};

// =========================
// Restriction relation OSM (via=node uniquement)
// =========================
struct RestrRel {
    int64_t from_way = 0;
    int64_t to_way = 0;
    int64_t via_node = 0;
    std::string kind;
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
// PASS 2 : collecter coords + node_index des nodes de graphe
// =========================
class NodeCollector : public osmium::handler::Handler {
public:
    const std::unordered_map<osmium::object_id_type, uint8_t>& use;
    const std::unordered_set<osmium::object_id_type>& endpoints;

    std::unordered_map<osmium::object_id_type, CoordE7>& coords_all;
    std::unordered_map<osmium::object_id_type, uint32_t>& node_index;

    NodeCollector(
        const std::unordered_map<osmium::object_id_type, uint8_t>& u,
        const std::unordered_set<osmium::object_id_type>& e,
        std::unordered_map<osmium::object_id_type, CoordE7>& ca,
        std::unordered_map<osmium::object_id_type, uint32_t>& ni
    ) : use(u), endpoints(e), coords_all(ca), node_index(ni) {}

    void node(const osmium::Node& n) {
        if (!n.location().valid()) return;
        const auto id = n.id();
        if (use.find(id) == use.end()) return;

        coords_all[id] = CoordE7{
            static_cast<int32_t>(std::llround(n.location().lat() * 1e7)),
            static_cast<int32_t>(std::llround(n.location().lon() * 1e7))
        };

        bool is_graph_node = false;
        auto it = use.find(id);
        if (it != use.end() && it->second >= 2) is_graph_node = true;
        if (endpoints.find(id) != endpoints.end()) is_graph_node = true;

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
// PASS 3 : découper ways, sens, arcs, et stocker géométrie par arc (pour export FGB)
// =========================
struct EdgeGeom {
    // géométrie complète de l’arc : liste [lon,lat]
    std::vector<std::pair<double,double>> lonlat;
    // attributs
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

    explicit WayToArcs(
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

        osmium::object_id_type start_node_id = 0;
        bool has_start = false;
        double acc_m = 0.0;

        std::vector<std::pair<double,double>> segment;

        const auto id0 = ns.front().ref();
        if (node_index.find(id0) != node_index.end()) {
            start_node_id = id0;
            has_start = true;
            segment.clear();
            std::pair<double,double> p0;
            if (!get_lonlat(id0, p0)) {
                has_start = false;
            } else {
                segment.push_back(p0);
            }
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
                    start_node_id = b_id;
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
                    start_node_id = b_id;
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
                if (has_start && start_node_id != b_id && segment.size() >= 2) {
                    auto su_it = node_index.find(start_node_id);
                    auto sv_it = node_index.find(b_id);
                    if (su_it != node_index.end() && sv_it != node_index.end()) {
                        const uint32_t u = su_it->second;
                        const uint32_t v = sv_it->second;
                        const uint32_t len_cm = static_cast<uint32_t>(std::llround(acc_m * 100.0));

                        auto fill_geom = [&](EdgeGeom& eg, const osmium::Way& wref, const char* dir) {
                            eg.lonlat = segment;
                            if (dir && std::string(dir) == "backward") {
                                std::reverse(eg.lonlat.begin(), eg.lonlat.end());
                            }
                            const char* highway  = wref.tags()["highway"];
                            const char* name     = wref.tags()["name"];
                            const char* ref      = wref.tags()["ref"];
                            const char* maxspeed = wref.tags()["maxspeed"];
                            const char* oneway   = wref.tags()["oneway"];
                            const char* junction = wref.tags()["junction"];
                            eg.highway  = highway ? highway : "";
                            eg.name     = name ? name : "";
                            eg.ref      = ref ? ref : "";
                            eg.maxspeed = maxspeed ? maxspeed : "";
                            eg.oneway   = oneway ? oneway : "";
                            eg.junction = junction ? junction : "";
                        };

                        if (forward) {
                            arcs.push_back({u, v, len_cm, 0, w.id()});
                            arc_geom.emplace_back();
                            fill_geom(arc_geom.back(), w, "forward");
                        }
                        if (backward) {
                            arcs.push_back({v, u, len_cm, 0, w.id()});
                            arc_geom.emplace_back();
                            fill_geom(arc_geom.back(), w, "backward");
                        }
                    }
                }

                start_node_id = b_id;
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
// PASS 4 : lire relations restriction (via=node)
// =========================
class RestrictionCollector : public osmium::handler::Handler {
public:
    std::vector<RestrRel>& restrictions;

    uint64_t ignored_no_tag = 0;
    uint64_t ignored_bad_members = 0;
    uint64_t ignored_via_way = 0;

    explicit RestrictionCollector(std::vector<RestrRel>& r) : restrictions(r) {}

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
    }
};

// =========================
// CSR builder
// =========================
static void build_csr(
    uint32_t N,
    const std::vector<ArcTmp>& arcs,
    std::vector<uint32_t>& offsets,
    std::vector<uint32_t>& to,
    std::vector<uint32_t>& len_cm,
    std::vector<uint8_t>& cls,
    std::vector<int64_t>& wid
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

    auto cur = offsets;
    for (uint32_t i = 0; i < M; ++i) {
        const auto& a = arcs[i];
        const uint32_t p = cur[a.u]++;
        to[p] = a.v;
        len_cm[p] = a.len_cm;
        cls[p] = a.hw_class;
        wid[p] = a.way_id;
    }
}

// incoming arcs index (rev)
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
            rev_arc[pos] = p; // arc index
        }
    }
}

static bool is_only_restriction(const std::string& k) { return k.rfind("only_", 0) == 0; }
static bool is_no_restriction(const std::string& k) { return k.rfind("no_", 0) == 0; }

// transitions interdites (arc_in, arc_out) + index vers restriction source
struct ForbiddenTurn {
    uint32_t in_arc = 0;
    uint32_t out_arc = 0;
    uint32_t restr_idx = 0;
};

static void build_forbidden_turns(
    const std::vector<RestrRel>& rels,
    const std::unordered_map<osmium::object_id_type, uint32_t>& node_index,
    const std::vector<uint32_t>& offsets,
    const std::vector<uint32_t>& to,
    const std::vector<int64_t>& wid,
    const std::vector<uint32_t>& rev_off,
    const std::vector<uint32_t>& rev_arc,
    std::vector<ForbiddenTurn>& forbidden
) {
    forbidden.clear();
    forbidden.reserve(rels.size() * 4);

    std::unordered_set<uint64_t> seen;
    seen.reserve(rels.size() * 8);

    for (uint32_t ri = 0; ri < static_cast<uint32_t>(rels.size()); ++ri) {
        const auto& rr = rels[ri];

        auto it_v = node_index.find(rr.via_node);
        if (it_v == node_index.end()) continue;
        const uint32_t via = it_v->second;

        std::vector<uint32_t> in_arcs;
        std::vector<uint32_t> out_to;
        std::vector<uint32_t> out_all;

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
                    if (seen.insert(key).second) forbidden.push_back({a, b, ri});
                }
            }
        } else if (is_only_restriction(rr.kind)) {
            if (out_to.empty()) continue;
            std::unordered_set<uint32_t> allowed(out_to.begin(), out_to.end());
            for (uint32_t a : in_arcs) {
                for (uint32_t b : out_all) {
                    if (allowed.find(b) != allowed.end()) continue;
                    const uint64_t key = (static_cast<uint64_t>(a) << 32) | b;
                    if (seen.insert(key).second) forbidden.push_back({a, b, ri});
                }
            }
        }
    }
}

// =========================
// Graph file v2 (avec restrictions)
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
// FlatGeobuf export (nodes, edges, restrictions)
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
    if (!drv) throw std::runtime_error("GDAL FlatGeobuf driver not available"); // :contentReference[oaicite:1]{index=1}

    GDALDataset* ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (!ds) throw std::runtime_error("Cannot create FlatGeobuf file for nodes");

    OGRSpatialReference srs;
    srs.SetWellKnownGeogCS("WGS84");

    OGRLayer* layer = ds->CreateLayer("node", &srs, wkbPoint, nullptr);
    if (!layer) { GDALClose(ds); throw std::runtime_error("Cannot create nodes layer"); }

    OGRFieldDefn f_id("id", OFTInteger);
    OGRFieldDefn f_osm("osm_node_id", OFTInteger64);
    layer->CreateField(&f_id);
    layer->CreateField(&f_osm);

    OGRFeatureDefn* defn = layer->GetLayerDefn();

    for (int i = 0; i < static_cast<int>(osm_node_id_by_idx.size()); ++i) {
        if (osm_node_id_by_idx[i] < 0) continue;

        OGRFeature* feat = OGRFeature::CreateFeature(defn);
        feat->SetField("id", i);
        feat->SetField("osm_node_id", static_cast<GIntBig>(osm_node_id_by_idx[i]));

        OGRPoint pt(lon_by_idx[i], lat_by_idx[i]);
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

static void export_edges_fgb(
    const std::string& path,
    const std::vector<ArcTmp>& arcs,
    const std::vector<EdgeGeom>& geom
) {
    GDALAllRegister();
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("FlatGeobuf");
    if (!drv) throw std::runtime_error("GDAL FlatGeobuf driver not available"); // :contentReference[oaicite:2]{index=2}

    GDALDataset* ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (!ds) throw std::runtime_error("Cannot create FlatGeobuf file for edges");

    OGRSpatialReference srs;
    srs.SetWellKnownGeogCS("WGS84");

    OGRLayer* layer = ds->CreateLayer("edge", &srs, wkbLineString, nullptr);
    if (!layer) { GDALClose(ds); throw std::runtime_error("Cannot create edges layer"); }

    layer->CreateField(new OGRFieldDefn("id", OFTInteger));
    layer->CreateField(new OGRFieldDefn("u", OFTInteger));
    layer->CreateField(new OGRFieldDefn("v", OFTInteger));
    layer->CreateField(new OGRFieldDefn("len_cm", OFTInteger));
    layer->CreateField(new OGRFieldDefn("way_id", OFTInteger64));

    layer->CreateField(new OGRFieldDefn("highway", OFTString));
    layer->CreateField(new OGRFieldDefn("name", OFTString));
    layer->CreateField(new OGRFieldDefn("ref", OFTString));
    layer->CreateField(new OGRFieldDefn("maxspeed", OFTString));
    layer->CreateField(new OGRFieldDefn("oneway", OFTString));
    layer->CreateField(new OGRFieldDefn("junction", OFTString));

    OGRFeatureDefn* defn = layer->GetLayerDefn();

    for (int i = 0; i < static_cast<int>(arcs.size()); ++i) {
        const auto& a = arcs[i];
        const auto& g = geom[i];

        if (g.lonlat.size() < 2) continue;

        OGRFeature* feat = OGRFeature::CreateFeature(defn);
        feat->SetField("id", i);
        feat->SetField("u", static_cast<int>(a.u));
        feat->SetField("v", static_cast<int>(a.v));
        feat->SetField("len_cm", static_cast<int>(a.len_cm));
        feat->SetField("way_id", static_cast<GIntBig>(a.way_id));

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

    GDALClose(ds);
}

static void export_restrictions_fgb(
    const std::string& path,
    const std::vector<ForbiddenTurn>& forbidden,
    const std::vector<RestrRel>& rels,
    const std::unordered_map<osmium::object_id_type, uint32_t>& node_index,
    const std::unordered_map<osmium::object_id_type, CoordE7>& coords_all
) {
    GDALAllRegister();
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("FlatGeobuf");
    if (!drv) throw std::runtime_error("GDAL FlatGeobuf driver not available"); // :contentReference[oaicite:3]{index=3}

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
    layer->CreateField(new OGRFieldDefn("via_node", OFTInteger64));
    layer->CreateField(new OGRFieldDefn("in_arc", OFTInteger));
    layer->CreateField(new OGRFieldDefn("out_arc", OFTInteger));

    OGRFeatureDefn* defn = layer->GetLayerDefn();

    for (int i = 0; i < static_cast<int>(forbidden.size()); ++i) {
        const auto& ft = forbidden[i];
        if (ft.restr_idx >= rels.size()) continue;

        const auto& rr = rels[ft.restr_idx];

        auto itc = coords_all.find(rr.via_node);
        if (itc == coords_all.end()) continue;

        const double lat = e7_to_deg(itc->second.lat_e7);
        const double lon = e7_to_deg(itc->second.lon_e7);

        OGRFeature* feat = OGRFeature::CreateFeature(defn);
        feat->SetField("id", i);
        feat->SetField("kind", rr.kind.c_str());
        feat->SetField("from_way", static_cast<GIntBig>(rr.from_way));
        feat->SetField("to_way", static_cast<GIntBig>(rr.to_way));
        feat->SetField("via_node", static_cast<GIntBig>(rr.via_node));
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
//   build_graph input.osm.pbf output.graph [--gfb] [--restriction-gfb=FILE]
// sorties:
//   output.graph (format GRF2 avec restrictions)
//   <base>_node.gfb, <base>_edge.gfb (si --gfb)
//   alsace_restriction.gfb (si --restriction-gfb=...)
// =========================
int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: build_graph input.osm.pbf output.graph [--gfb] [--restriction-gfb=FILE]\n";
        return 1;
    }

    const std::string in = argv[1];
    const std::string out_graph = argv[2];

    // Objectif : --gfb génère automatiquement 3 fichiers :
    //   <base>_node.gfb
    //   <base>_edge.gfb
    //   <base>_restriction.gfb

    bool emit_gfb = false;
    for (int i = 3; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--fgb" || a == "-fgb") emit_gfb = true;
    }

    NodeUseCounter counter;
    apply_with_progress(in, "PASS 1: counting node usage", counter);

    std::unordered_map<osmium::object_id_type, CoordE7> coords_all;
    coords_all.reserve(counter.use.size());

    std::unordered_map<osmium::object_id_type, uint32_t> node_index;
    node_index.reserve(counter.use.size() / 4);

    NodeCollector nc(counter.use, counter.endpoints, coords_all, node_index);
    apply_with_progress(in, "PASS 2: collecting coords + graph nodes", nc);

    std::vector<ArcTmp> arcs;
    std::vector<EdgeGeom> arc_geom;
    arcs.reserve(1'000'000);
    arc_geom.reserve(1'000'000);

    WayToArcs wta(node_index, coords_all, arcs, arc_geom);
    apply_with_progress(in, "PASS 3: building arcs (oneway + roundabout applied)", wta);

    std::vector<RestrRel> rels;
    rels.reserve(200000);
    RestrictionCollector rc(rels);
    apply_with_progress(in, "PASS 4: reading restrictions (via=node only)", rc);

    const uint32_t N = static_cast<uint32_t>(node_index.size());

    std::vector<uint32_t> offsets;
    std::vector<uint32_t> to;
    std::vector<uint32_t> len_cm;
    std::vector<uint8_t> cls;
    std::vector<int64_t> wid;
    build_csr(N, arcs, offsets, to, len_cm, cls, wid);

    std::vector<uint32_t> rev_off;
    std::vector<uint32_t> rev_arc;
    build_reverse_arc_index(N, offsets, to, rev_off, rev_arc);

    std::vector<ForbiddenTurn> forbidden;
    build_forbidden_turns(rels, node_index, offsets, to, wid, rev_off, rev_arc, forbidden);

    std::cerr << "N=" << N
              << " M=" << to.size()
              << " restrictions=" << rels.size()
              << " forbidden_pairs=" << forbidden.size()
              << " ignored(via_way)=" << rc.ignored_via_way
              << "\n";

    write_graph_v2_with_restrictions(out_graph, N, offsets, to, len_cm, cls, wid, forbidden);

    if (emit_gfb) {
        const std::string base = base_name_no_ext(out_graph);

        const std::string node_path = base + "_node.gfb";
        const std::string edge_path = base + "_edge.gfb";
        const std::string restr_path = base + "_restriction.gfb";

        std::vector<int64_t> osm_node_id_by_idx;
        std::vector<double> lon_by_idx;
        std::vector<double> lat_by_idx;
        build_graph_node_arrays(N, node_index, coords_all, osm_node_id_by_idx, lon_by_idx, lat_by_idx);

        export_nodes_fgb(node_path, osm_node_id_by_idx, lon_by_idx, lat_by_idx);
        export_edges_fgb(edge_path, arcs, arc_geom);
        export_restrictions_fgb(restr_path, forbidden, rels, node_index, coords_all);
    }

    return 0;
}
