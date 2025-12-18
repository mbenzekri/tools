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

// =======================================================
// Structure temporaire d'arc
// =======================================================
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
// Appliquer un handler avec barre de progression (correct libosmium)
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
// PASS 1 : compter usage des nodes + mémoriser extrémités de ways
// =======================================================
class NodeUseCounter : public osmium::handler::Handler {
public:
    // compteur saturé à 2 : 0/1/2 (2 = partagé)
    std::unordered_map<osmium::object_id_type, uint8_t> use;
    std::unordered_set<osmium::object_id_type> endpoints;

    void way(const osmium::Way& w) {
        if (!w.tags().has_key("highway")) return;

        const auto& nodes = w.nodes();
        if (nodes.size() < 2) return;

        // extrémités
        endpoints.insert(nodes.front().ref());
        endpoints.insert(nodes.back().ref());

        // usage nodes
        for (const auto& n : nodes) {
            uint8_t& c = use[n.ref()];
            if (c < 2) ++c;
        }
    }
};

// =======================================================
// PASS 2 : collecter coordonnées de tous nodes utiles
//          + créer node_index uniquement pour "nodes de graphe"
// =======================================================
class NodeCollector : public osmium::handler::Handler {
public:
    const std::unordered_map<osmium::object_id_type, uint8_t>& use;
    const std::unordered_set<osmium::object_id_type>& endpoints;

    // coordonnées pour tous nodes utilisés par highways
    std::unordered_map<osmium::object_id_type, CoordE7>& coords_all;

    // index compact 0..N-1 uniquement pour nodes de graphe
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

        // On ne garde que les nodes présents dans les ways highway (use map)
        if (use.find(id) == use.end()) return;

        // stocker coord pour calcul des longueurs
        // (OSM stocke coord en double ; on stocke en e7 pour compacité)
        const double lat = n.location().lat();
        const double lon = n.location().lon();
        CoordE7 ce7{
            static_cast<int32_t>(std::llround(lat * 1e7)),
            static_cast<int32_t>(std::llround(lon * 1e7))
        };
        coords_all[id] = ce7;

        // node de graphe si partagé (use>=2) OU extrémité
        bool is_graph_node = false;
        auto it = use.find(id);
        if (it != use.end() && it->second >= 2) is_graph_node = true;
        if (endpoints.find(id) != endpoints.end()) is_graph_node = true;

        if (is_graph_node) {
            if (node_index.find(id) == node_index.end()) {
                uint32_t idx = static_cast<uint32_t>(node_index.size());
                node_index[id] = idx;
            }
        }
    }
};

// =======================================================
// PASS 3 : découper ways aux nodes de graphe, générer arcs
// =======================================================
class WayToArcs : public osmium::handler::Handler {
public:
    const std::unordered_map<osmium::object_id_type, uint32_t>& node_index;
    const std::unordered_map<osmium::object_id_type, CoordE7>& coords_all;
    std::vector<ArcTmp>& arcs;

    WayToArcs(
        const std::unordered_map<osmium::object_id_type, uint32_t>& ni,
        const std::unordered_map<osmium::object_id_type, CoordE7>& ca,
        std::vector<ArcTmp>& a
    ) : node_index(ni), coords_all(ca), arcs(a) {}

    void way(const osmium::Way& w) {
        if (!w.tags().has_key("highway")) return;

        const auto& ns = w.nodes();
        if (ns.size() < 2) return;

        // On veut couper uniquement quand on rencontre un node_index (node de graphe)
        // start_node_id = dernier node de graphe rencontré
        osmium::object_id_type start_node_id = 0;
        bool has_start = false;
        double acc_m = 0.0;

        // initialiser au premier node s'il est un node de graphe
        {
            auto id0 = ns.front().ref();
            if (node_index.find(id0) != node_index.end()) {
                start_node_id = id0;
                has_start = true;
            }
        }

        // parcourir segments consécutifs
        for (size_t i = 1; i < ns.size(); ++i) {
            const auto a_id = ns[i-1].ref();
            const auto b_id = ns[i].ref();

            auto ca_it = coords_all.find(a_id);
            auto cb_it = coords_all.find(b_id);
            if (ca_it == coords_all.end() || cb_it == coords_all.end()) {
                // coord manquante => on ne peut pas calculer la longueur
                // on coupe la continuité
                acc_m = 0.0;
                has_start = false;
                if (node_index.find(b_id) != node_index.end()) {
                    start_node_id = b_id;
                    has_start = true;
                }
                continue;
            }

            const auto& A = ca_it->second;
            const auto& B = cb_it->second;
            acc_m += haversine(
                e7_to_deg(A.lat_e7), e7_to_deg(A.lon_e7),
                e7_to_deg(B.lat_e7), e7_to_deg(B.lon_e7)
            );

            // si b_id est un node de graphe, on ferme un tronçon start -> b
            if (node_index.find(b_id) != node_index.end()) {
                if (has_start && start_node_id != b_id) {
                    auto su_it = node_index.find(start_node_id);
                    auto sv_it = node_index.find(b_id);
                    if (su_it != node_index.end() && sv_it != node_index.end()) {
                        uint32_t u = su_it->second;
                        uint32_t v = sv_it->second;

                        uint32_t len_cm = static_cast<uint32_t>(std::llround(acc_m * 100.0));

                        // pour l'instant : double sens (pas de oneway)
                        arcs.push_back({u, v, len_cm, 0, w.id()});
                        arcs.push_back({v, u, len_cm, 0, w.id()});
                    }
                }

                // reset
                start_node_id = b_id;
                has_start = true;
                acc_m = 0.0;
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
    if (!out) throw std::runtime_error("cannot open output file");

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
        std::cerr << "usage: build_graph input.osm.pbf output.graph\n";
        return 1;
    }

    const std::string in  = argv[1];
    const std::string out = argv[2];

    std::cerr << "Input file : " << in << "\n";

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

    // PASS 3
    std::vector<ArcTmp> arcs;
    arcs.reserve(1'000'000); // réserve indicative, ça grossira si besoin

    WayToArcs wta(node_index, coords_all, arcs);
    apply_with_progress(in, "PASS 3: building arcs", wta);

    uint32_t N = static_cast<uint32_t>(node_index.size());
    uint32_t M = static_cast<uint32_t>(arcs.size());

    std::cerr << "\n==== GRAPH STATS ====\n";
    std::cerr << "N (graph nodes) : " << N << "\n";
    std::cerr << "M (directed arcs): " << M << "\n";
    std::cerr << "=====================\n";

    write_graph(out, N, arcs);

    std::cerr << "Graph written to: " << out << "\n";
    return 0;
}

