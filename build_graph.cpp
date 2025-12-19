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

// =======================================================
// Constantes
// =======================================================
static constexpr uint64_t DIST_INF = UINT64_MAX;
static constexpr uint32_t NO_NODE  = UINT32_MAX;

// =======================================================
// Coordonnées stockées en entiers (lat/lon * 1e7)
// =======================================================
struct CoordE7 {
    int32_t lat_e7;
    int32_t lon_e7;
};

struct ArcTmp {
    uint32_t u;
    uint32_t v;
    uint32_t len_cm;
    uint8_t  hw_class;
    int64_t  way_id;
};

// =======================================================
// Haversine
// =======================================================
static double haversine(double lat1, double lon1,
                        double lat2, double lon2) {
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

static double e7_to_deg(int32_t x) {
    return static_cast<double>(x) / 1e7;
}

// =======================================================
// Appliquer un handler avec barre de progression
// =======================================================
template <typename Handler>
void apply_with_progress(const std::string& filename,
                         const std::string& label,
                         Handler& handler) {
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

// =======================================================
// PASS 1 : compter usage des nodes + extrémités
// =======================================================
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

// =======================================================
// PASS 2 : coords de tous nodes utiles + index des nodes de graphe
// =======================================================
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

// =======================================================
// GeoJSONSeq writers (1 Feature / ligne)
// =======================================================
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

static void write_nodes_geojsonseq(
    const std::string& path,
    const std::vector<int64_t>& osm_node_id_by_idx,
    const std::vector<double>& lon_by_idx,
    const std::vector<double>& lat_by_idx
) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open node geojsonseq output");

    out.setf(std::ios::fixed);
    out.precision(7);

    for (size_t i = 0; i < osm_node_id_by_idx.size(); ++i) {
        if (osm_node_id_by_idx[i] < 0) continue;
        out << "{\"type\":\"Feature\",\"properties\":{\"id\":" << i
            << ",\"osm_node_id\":" << osm_node_id_by_idx[i]
            << "},\"geometry\":{\"type\":\"Point\",\"coordinates\":["
            << lon_by_idx[i] << "," << lat_by_idx[i] << "]}}\n";
    }
}

// échappement JSON minimal pour les valeurs string (guillemets + backslash)
static std::string json_escape(const char* s) {
    if (!s) return "";
    std::string out;
    for (const char* p = s; *p; ++p) {
        if (*p == '\"' || *p == '\\') out.push_back('\\');
        out.push_back(*p);
    }
    return out;
}

// =======================================================
// PASS 3 : découper ways et générer arcs
// + si GeoJSONSeq activé : écrire edges avec géométrie complète
// =======================================================
class WayToArcs : public osmium::handler::Handler {
public:
    const std::unordered_map<osmium::object_id_type, uint32_t>& node_index;
    const std::unordered_map<osmium::object_id_type, CoordE7>& coords_all;
    std::vector<ArcTmp>& arcs;

    // GeoJSONSeq edge stream (optionnel)
    std::ofstream* edge_seq = nullptr;
    const std::vector<double>* lon_by_idx = nullptr;
    const std::vector<double>* lat_by_idx = nullptr;

    WayToArcs(
        const std::unordered_map<osmium::object_id_type, uint32_t>& ni,
        const std::unordered_map<osmium::object_id_type, CoordE7>& ca,
        std::vector<ArcTmp>& a
    ) : node_index(ni), coords_all(ca), arcs(a) {}

    void enable_geojsonseq_edges(std::ofstream& out,
                                const std::vector<double>& lon,
                                const std::vector<double>& lat) {
        edge_seq = &out;
        lon_by_idx = &lon;
        lat_by_idx = &lat;
        edge_seq->setf(std::ios::fixed);
        edge_seq->precision(7);
    }

    void emit_edge_feature(uint64_t fid,
                           const ArcTmp& arc,
                           const std::vector<std::pair<double,double>>& coords,
                           const osmium::Way& w) {
        if (!edge_seq) return;

        const char* highway = w.tags()["highway"];
        const char* name    = w.tags()["name"];
        const char* ref     = w.tags()["ref"];
        const char* maxsp   = w.tags()["maxspeed"];
        const char* oneway  = w.tags()["oneway"];
        const char* junction= w.tags()["junction"];

        *edge_seq << "{\"type\":\"Feature\",\"properties\":{"
                  << "\"id\":" << fid
                  << ",\"u\":" << arc.u
                  << ",\"v\":" << arc.v
                  << ",\"len_cm\":" << arc.len_cm
                  << ",\"way_id\":" << arc.way_id;

        if (highway)  *edge_seq << ",\"highway\":\""  << json_escape(highway)  << "\"";
        if (name)     *edge_seq << ",\"name\":\""     << json_escape(name)     << "\"";
        if (ref)      *edge_seq << ",\"ref\":\""      << json_escape(ref)      << "\"";
        if (maxsp)    *edge_seq << ",\"maxspeed\":\"" << json_escape(maxsp)    << "\"";
        if (oneway)   *edge_seq << ",\"oneway\":\""   << json_escape(oneway)   << "\"";
        if (junction) *edge_seq << ",\"junction\":\"" << json_escape(junction) << "\"";

        *edge_seq << "},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[";

        bool first = true;
        for (const auto& c : coords) {
            if (!first) *edge_seq << ",";
            first = false;
            *edge_seq << "[" << c.first << "," << c.second << "]";
        }

        *edge_seq << "]}}\n";
    }

    void way(const osmium::Way& w) {
        if (!w.tags().has_key("highway")) return;

        const auto& ns = w.nodes();
        if (ns.size() < 2) return;

        osmium::object_id_type start_node_id = 0;
        bool has_start = false;
        double acc_m = 0.0;

        // géométrie complète du tronçon courant (lon,lat)
        std::vector<std::pair<double,double>> segment_coords;

        auto push_coord = [&](osmium::object_id_type osm_node_id) -> bool {
            auto it = coords_all.find(osm_node_id);
            if (it == coords_all.end()) return false;
            double lat = e7_to_deg(it->second.lat_e7);
            double lon = e7_to_deg(it->second.lon_e7);
            segment_coords.emplace_back(lon, lat); // GeoJSON: [lon,lat]
            return true;
        };

        const auto id0 = ns.front().ref();
        if (node_index.find(id0) != node_index.end()) {
            start_node_id = id0;
            has_start = true;
            segment_coords.clear();
            if (!push_coord(id0)) {
                has_start = false;
                segment_coords.clear();
            }
        }

        for (size_t i = 1; i < ns.size(); ++i) {
            const auto a_id = ns[i-1].ref();
            const auto b_id = ns[i].ref();

            auto ca_it = coords_all.find(a_id);
            auto cb_it = coords_all.find(b_id);
            if (ca_it == coords_all.end() || cb_it == coords_all.end()) {
                // cassure: reset
                acc_m = 0.0;
                has_start = false;
                segment_coords.clear();
                if (node_index.find(b_id) != node_index.end()) {
                    start_node_id = b_id;
                    has_start = push_coord(b_id);
                }
                continue;
            }

            const auto& A = ca_it->second;
            const auto& B = cb_it->second;
            acc_m += haversine(
                e7_to_deg(A.lat_e7), e7_to_deg(A.lon_e7),
                e7_to_deg(B.lat_e7), e7_to_deg(B.lon_e7)
            );

            // accumuler géométrie (tous points)
            if (has_start) {
                // on pousse b_id (y compris points intermédiaires)
                // éviter duplicat: si segment_coords déjà a_id, pousser b_id suffit
                if (!push_coord(b_id)) {
                    acc_m = 0.0;
                    has_start = false;
                    segment_coords.clear();
                }
            } else {
                // si on n'a pas de start, mais b_id est un node de graphe, on repart
                if (node_index.find(b_id) != node_index.end()) {
                    start_node_id = b_id;
                    has_start = true;
                    segment_coords.clear();
                    if (!push_coord(b_id)) {
                        has_start = false;
                        segment_coords.clear();
                    }
                    acc_m = 0.0;
                }
                continue;
            }

            // fermeture au prochain node de graphe
            if (node_index.find(b_id) != node_index.end()) {
                if (start_node_id != b_id) {
                    auto su_it = node_index.find(start_node_id);
                    auto sv_it = node_index.find(b_id);
                    if (su_it != node_index.end() && sv_it != node_index.end()) {
                        uint32_t u = su_it->second;
                        uint32_t v = sv_it->second;
                        uint32_t len_cm = static_cast<uint32_t>(std::llround(acc_m * 100.0));

                        // arcs (double sens)
                        ArcTmp ab{u, v, len_cm, 0, w.id()};
                        ArcTmp ba{v, u, len_cm, 0, w.id()};
                        arcs.push_back(ab);
                        arcs.push_back(ba);

                        // GeoJSONSeq: écrire la géométrie complète pour chaque arc
                        // NB: pour l'arc inverse, on écrit la géométrie inversée
                        const uint64_t fid_ab = static_cast<uint64_t>(arcs.size() - 2);
                        emit_edge_feature(fid_ab, ab, segment_coords, w);

                        if (edge_seq) {
                            std::vector<std::pair<double,double>> rev(segment_coords.rbegin(), segment_coords.rend());
                            const uint64_t fid_ba = static_cast<uint64_t>(arcs.size() - 1);
                            emit_edge_feature(fid_ba, ba, rev, w);
                        }
                    }
                }

                // reset segment: nouveau départ à b_id
                start_node_id = b_id;
                has_start = true;
                acc_m = 0.0;
                segment_coords.clear();
                // b_id devient premier point du prochain tronçon
                if (!push_coord(b_id)) {
                    has_start = false;
                    segment_coords.clear();
                }
            }
        }
    }
};

// =======================================================
// Écriture du graphe binaire complet
// =======================================================
void write_graph(const std::string& path,
                 uint32_t N,
                 const std::vector<ArcTmp>& arcs) {
    uint32_t M = static_cast<uint32_t>(arcs.size());

    std::vector<uint32_t> deg(N, 0);
    for (const auto& a : arcs) ++deg[a.u];

    std::vector<uint32_t> offsets(N+1, 0);
    for (uint32_t i = 0; i < N; ++i)
        offsets[i+1] = offsets[i] + deg[i];

    std::vector<uint32_t> to(M);
    std::vector<uint32_t> len(M);
    std::vector<uint8_t>  cls(M);
    std::vector<int64_t>  wid(M);

    auto cur = offsets;
    for (const auto& a : arcs) {
        uint32_t p = cur[a.u]++;
        to[p]  = a.v;
        len[p] = a.len_cm;
        cls[p] = a.hw_class;
        wid[p] = a.way_id;
    }

    std::vector<uint64_t> dist(N, DIST_INF);
    std::vector<uint32_t> label(N, NO_NODE);
    std::vector<uint32_t> pred(N, NO_NODE);

    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("cannot open output graph file");

    out.write((char*)&N,4);
    out.write((char*)&M,4);
    out.write((char*)offsets.data(), offsets.size()*4);
    out.write((char*)to.data(),      to.size()*4);
    out.write((char*)len.data(),     len.size()*4);
    out.write((char*)cls.data(),     cls.size());
    out.write((char*)wid.data(),     wid.size()*8);
    out.write((char*)dist.data(),  dist.size()*8);
    out.write((char*)label.data(), label.size()*4);
    out.write((char*)pred.data(),  pred.size()*4);
}

// =======================================================
// MAIN
// =======================================================
int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: build_graph input.osm.pbf output.graph [--geojsonseq]\n";
        return 1;
    }

    const std::string in  = argv[1];
    const std::string out_graph = argv[2];

    bool emit_geojsonseq = false;
    for (int i = 3; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--geojsonseq") emit_geojsonseq = true;
    }

    std::cerr << "Input  : " << in << "\n";
    std::cerr << "Output : " << out_graph << "\n";
    if (emit_geojsonseq) std::cerr << "GeoJSONSeq: enabled (full geometry + attributes)\n";

    // PASS 1
    NodeUseCounter counter;
    apply_with_progress(in, "PASS 1: counting node usage", counter);

    // PASS 2
    std::unordered_map<osmium::object_id_type, CoordE7> coords_all;
    coords_all.reserve(counter.use.size());

    std::unordered_map<osmium::object_id_type, uint32_t> node_index;
    node_index.reserve(counter.use.size() / 4);

    NodeCollector nc(counter.use, counter.endpoints, coords_all, node_index);
    apply_with_progress(in, "PASS 2: collecting coords + graph nodes", nc);

    // Préparer GeoJSONSeq nodes (nécessite arrays index->coord)
    std::vector<int64_t> osm_node_id_by_idx;
    std::vector<double> lon_by_idx;
    std::vector<double> lat_by_idx;

    std::ofstream edge_out;
    if (emit_geojsonseq) {
        const std::string base = base_name_no_ext(out_graph);
        const std::string node_path = base + "_node.geojsonseq";
        const std::string edge_path = base + "_edge.geojsonseq";

        build_graph_node_arrays(
            static_cast<uint32_t>(node_index.size()),
            node_index, coords_all,
            osm_node_id_by_idx, lon_by_idx, lat_by_idx
        );

        std::cerr << "Writing nodes GeoJSONSeq: " << node_path << "\n";
        write_nodes_geojsonseq(node_path, osm_node_id_by_idx, lon_by_idx, lat_by_idx);

        std::cerr << "Writing edges GeoJSONSeq: " << edge_path << "\n";
        edge_out.open(edge_path);
        if (!edge_out) throw std::runtime_error("cannot open edge geojsonseq output");
    }

    // PASS 3
    std::vector<ArcTmp> arcs;
    arcs.reserve(1'000'000);

    WayToArcs wta(node_index, coords_all, arcs);
    if (emit_geojsonseq) {
        wta.enable_geojsonseq_edges(edge_out, lon_by_idx, lat_by_idx);
    }
    apply_with_progress(in, "PASS 3: building arcs", wta);

    if (emit_geojsonseq) {
        edge_out.close();
        std::cerr << "GeoJSONSeq export done.\n";
    }

    const uint32_t N = static_cast<uint32_t>(node_index.size());
    const uint32_t M = static_cast<uint32_t>(arcs.size());

    std::cerr << "\n==== GRAPH STATS ====\n";
    std::cerr << "N (graph nodes) : " << N << "\n";
    std::cerr << "M (directed arcs): " << M << "\n";
    std::cerr << "=====================\n";

    write_graph(out_graph, N, arcs);
    std::cerr << "Graph written to: " << out_graph << "\n";

    return 0;
}
