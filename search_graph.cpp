#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <sstream>
#include <filesystem>

#include "gdal.h"
#include "ogrsf_frmts.h"

static constexpr uint64_t DIST_INF = std::numeric_limits<uint64_t>::max();
static constexpr uint32_t NO_NODE  = std::numeric_limits<uint32_t>::max();

static void read_exact(std::ifstream& in, void* dst, std::size_t n) {
    in.read(reinterpret_cast<char*>(dst), static_cast<std::streamsize>(n));
    if (!in || static_cast<std::size_t>(in.gcount()) != n) {
        throw std::runtime_error("read_exact failed (file truncated or IO error)");
    }
}

struct Graph {
    uint32_t N = 0;
    uint32_t M = 0;

    std::vector<uint32_t> offsets; // N+1
    std::vector<uint32_t> to;      // M
    std::vector<uint32_t> len_cm;  // M
    std::vector<uint8_t>  cls;     // M
    std::vector<int64_t>  way_id;  // M

    std::vector<uint64_t> dist;    // N
    std::vector<uint32_t> label;   // N
    std::vector<uint32_t> pred;    // N

    // GRF2 only
    uint32_t R = 0;
    std::vector<std::pair<uint32_t,uint32_t>> forbidden; // R
};

static bool starts_with_magic_grf2(std::ifstream& in) {
    char magic[4];
    read_exact(in, magic, 4);
    in.seekg(0, std::ios::beg);
    return magic[0]=='G' && magic[1]=='R' && magic[2]=='F' && magic[3]=='2';
}

static Graph read_graph(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open graph file");

    Graph g;

    if (starts_with_magic_grf2(in)) {
        char magic[4];
        read_exact(in, magic, 4); // "GRF2"

        read_exact(in, &g.N, 4);
        read_exact(in, &g.M, 4);

        g.offsets.resize(static_cast<std::size_t>(g.N) + 1);
        g.to.resize(g.M);
        g.len_cm.resize(g.M);
        g.cls.resize(g.M);
        g.way_id.resize(g.M);
        g.dist.resize(g.N);
        g.label.resize(g.N);
        g.pred.resize(g.N);

        read_exact(in, g.offsets.data(), g.offsets.size() * 4);
        read_exact(in, g.to.data(),      g.to.size() * 4);
        read_exact(in, g.len_cm.data(),  g.len_cm.size() * 4);
        read_exact(in, g.cls.data(),     g.cls.size());
        read_exact(in, g.way_id.data(),  g.way_id.size() * 8);

        read_exact(in, g.dist.data(),    g.dist.size() * 8);
        read_exact(in, g.label.data(),   g.label.size() * 4);
        read_exact(in, g.pred.data(),    g.pred.size() * 4);

        read_exact(in, &g.R, 4);
        g.forbidden.resize(g.R);
        for (uint32_t i = 0; i < g.R; ++i) {
            uint32_t a = 0, b = 0;
            read_exact(in, &a, 4);
            read_exact(in, &b, 4);
            g.forbidden[i] = {a, b};
        }
    } else {
        read_exact(in, &g.N, 4);
        read_exact(in, &g.M, 4);

        g.offsets.resize(static_cast<std::size_t>(g.N) + 1);
        g.to.resize(g.M);
        g.len_cm.resize(g.M);
        g.cls.resize(g.M);
        g.way_id.resize(g.M);
        g.dist.resize(g.N);
        g.label.resize(g.N);
        g.pred.resize(g.N);

        read_exact(in, g.offsets.data(), g.offsets.size() * 4);
        read_exact(in, g.to.data(),      g.to.size() * 4);
        read_exact(in, g.len_cm.data(),  g.len_cm.size() * 4);
        read_exact(in, g.cls.data(),     g.cls.size());
        read_exact(in, g.way_id.data(),  g.way_id.size() * 8);

        read_exact(in, g.dist.data(),    g.dist.size() * 8);
        read_exact(in, g.label.data(),   g.label.size() * 4);
        read_exact(in, g.pred.data(),    g.pred.size() * 4);

        g.R = 0;
    }

    return g;
}

static std::string base_name_no_ext(const std::string& graph_path) {
    std::filesystem::path p(graph_path);
    auto stem = p.stem().string();
    auto parent = p.parent_path();
    if (parent.empty()) return stem;
    return (parent / stem).string();
}

static std::vector<int64_t> read_node_map_text(const std::string& path, uint32_t N) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open node map file");

    std::vector<int64_t> osm_by_graph(N, -1);
    std::string line;
    std::size_t line_no = 0;

    while (std::getline(in, line)) {
        ++line_no;
        if (line.empty() || line[0] == '#') continue;
        for (char& c : line) {
            if (c == ',' || c == ';') c = ' ';
        }
        std::stringstream ss(line);
        long long a = 0, b = 0;
        if (!(ss >> a >> b)) continue;

        uint32_t graph_id = NO_NODE;
        int64_t osm_id = -1;

        if (a >= 0 && static_cast<uint64_t>(a) < N) {
            graph_id = static_cast<uint32_t>(a);
            osm_id = static_cast<int64_t>(b);
        } else if (b >= 0 && static_cast<uint64_t>(b) < N) {
            graph_id = static_cast<uint32_t>(b);
            osm_id = static_cast<int64_t>(a);
        } else {
            std::cerr << "Warning: bad node-map line " << line_no << " (skipped)\n";
            continue;
        }

        osm_by_graph[graph_id] = osm_id;
    }

    return osm_by_graph;
}

static std::vector<int64_t> read_node_map_fgb(const std::string& path, uint32_t N) {
    GDALAllRegister();
    GDALDataset* ds = static_cast<GDALDataset*>(GDALOpenEx(path.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr));
    if (!ds) throw std::runtime_error("cannot open node fgb file");

    OGRLayer* layer = ds->GetLayerByName("node");
    if (!layer) layer = ds->GetLayer(0);
    if (!layer) {
        GDALClose(ds);
        throw std::runtime_error("no layer in node fgb file");
    }

    std::vector<int64_t> osm_by_graph(N, -1);
    layer->ResetReading();

    OGRFeature* feat = nullptr;
    while ((feat = layer->GetNextFeature()) != nullptr) {
        const int graph_id = feat->GetFieldAsInteger("id");
        const GIntBig osm_id = feat->GetFieldAsInteger64("osm_node_id");
        if (graph_id >= 0 && static_cast<uint32_t>(graph_id) < N) {
            osm_by_graph[static_cast<uint32_t>(graph_id)] = static_cast<int64_t>(osm_id);
        }
        OGRFeature::DestroyFeature(feat);
    }

    GDALClose(ds);
    return osm_by_graph;
}

static std::vector<int64_t> read_node_map_auto(const std::string& graph_path, uint32_t N, const std::string& fgb_path) {
    if (!fgb_path.empty()) {
        return read_node_map_fgb(fgb_path, N);
    }

    const std::string base = base_name_no_ext(graph_path);
    const std::string default_fgb = base + "_node.fgb";
    if (std::filesystem::exists(default_fgb)) {
        return read_node_map_fgb(default_fgb, N);
    }

    throw std::runtime_error("node map not found (use --node-map or --node-fgb)");
}

static std::vector<int64_t> read_osm_node_list(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open input nodes file");

    std::vector<int64_t> ids;
    std::string line;
    std::size_t line_no = 0;

    while (std::getline(in, line)) {
        ++line_no;
        if (line.empty() || line[0] == '#') continue;
        for (char& c : line) {
            if (c == ',' || c == ';') c = ' ';
        }

        std::stringstream ss(line);
        std::string id_str;
        if (!(ss >> id_str)) continue;
        if (!id_str.empty() && id_str[0] == '#') continue;

        try {
            ids.push_back(std::stoll(id_str));
        } catch (...) {
            std::cerr << "Warning: bad node id on line " << line_no << " (skipped)\n";
        }
    }

    return ids;
}

struct RadixHeap {
    struct Item {
        uint64_t key = 0;
        uint32_t value = 0;
    };

    static constexpr int B = 65; // 0..64 for uint64 keys

    uint64_t last = 0;
    std::size_t size_ = 0;
    std::vector<std::vector<Item>> buckets;

    RadixHeap() : buckets(B) {}

    void clear() {
        for (auto& b : buckets) b.clear();
        size_ = 0;
        last = 0;
    }

    bool empty() const { return size_ == 0; }

    int bucket_index(uint64_t key) const {
        if (key == last) return 0;
        const uint64_t x = key ^ last;
        return 64 - __builtin_clzll(x);
    }

    void push(uint64_t key, uint32_t value) {
        if (key < last) throw std::runtime_error("radix heap key < last");
        const int b = bucket_index(key);
        buckets[b].push_back({key, value});
        ++size_;
    }

    Item pop() {
        if (buckets[0].empty()) {
            int i = 1;
            while (i < B && buckets[i].empty()) ++i;
            if (i == B) throw std::runtime_error("pop from empty radix heap");

            uint64_t new_last = buckets[i][0].key;
            for (const auto& it : buckets[i]) {
                if (it.key < new_last) new_last = it.key;
            }
            last = new_last;

            for (const auto& it : buckets[i]) {
                const int b = bucket_index(it.key);
                buckets[b].push_back(it);
            }
            buckets[i].clear();
        }

        Item it = buckets[0].back();
        buckets[0].pop_back();
        --size_;
        return it;
    }
};

static bool build_path(uint32_t s, uint32_t t, const std::vector<uint32_t>& pred, std::vector<uint32_t>& out) {
    out.clear();
    uint32_t v = t;
    while (v != NO_NODE) {
        out.push_back(v);
        if (v == s) break;
        v = pred[v];
    }
    if (out.empty() || out.back() != s) {
        out.clear();
        return false;
    }
    std::reverse(out.begin(), out.end());
    return true;
}

struct Args {
    std::string graph_path;
    std::string nodes_path;
    std::string output_path;
    std::string node_map_path;
    std::string node_fgb_path;
    bool out_binary = false;
    bool out_binary_set = false;
};

static Args parse_args(int argc, char** argv) {
    if (argc < 4) {
        throw std::runtime_error(
            "usage: search_graph graph.grf2 nodes_osm.txt output.[txt|bin] "
            "[--text|--bin] [--node-map map.txt] [--node-fgb nodes.fgb]"
        );
    }

    Args a;
    a.graph_path = argv[1];
    a.nodes_path = argv[2];
    a.output_path = argv[3];

    for (int i = 4; i < argc; ++i) {
        const std::string opt = argv[i];
        if (opt == "--text") {
            a.out_binary = false;
            a.out_binary_set = true;
        } else if (opt == "--bin" || opt == "--binary") {
            a.out_binary = true;
            a.out_binary_set = true;
        } else if (opt == "--node-map") {
            if (i + 1 >= argc) throw std::runtime_error("--node-map needs a file path");
            a.node_map_path = argv[++i];
        } else if (opt == "--node-fgb") {
            if (i + 1 >= argc) throw std::runtime_error("--node-fgb needs a file path");
            a.node_fgb_path = argv[++i];
        } else {
            throw std::runtime_error("unknown argument: " + opt);
        }
    }

    if (!a.out_binary_set) {
        const std::string ext = std::filesystem::path(a.output_path).extension().string();
        if (ext == ".bin") a.out_binary = true;
    }

    return a;
}

int main(int argc, char** argv) {
    try {
        const Args args = parse_args(argc, argv);
        Graph g = read_graph(args.graph_path);

        std::vector<int64_t> osm_by_graph;
        if (!args.node_map_path.empty()) {
            osm_by_graph = read_node_map_text(args.node_map_path, g.N);
        } else {
            osm_by_graph = read_node_map_auto(args.graph_path, g.N, args.node_fgb_path);
        }

        std::unordered_map<int64_t, uint32_t> graph_by_osm;
        graph_by_osm.reserve(g.N);
        for (uint32_t i = 0; i < g.N; ++i) {
            if (osm_by_graph[i] >= 0) graph_by_osm[osm_by_graph[i]] = i;
        }

        const std::vector<int64_t> osm_input = read_osm_node_list(args.nodes_path);
        if (osm_input.empty()) throw std::runtime_error("input nodes file is empty");

        std::vector<uint8_t> is_target(g.N, 0);
        std::vector<uint32_t> targets;
        targets.reserve(osm_input.size());

        std::size_t missing = 0;
        for (int64_t osm_id : osm_input) {
            auto it = graph_by_osm.find(osm_id);
            if (it == graph_by_osm.end()) {
                ++missing;
                continue;
            }
            const uint32_t gid = it->second;
            if (!is_target[gid]) {
                is_target[gid] = 1;
                targets.push_back(gid);
            }
        }

        if (missing > 0) {
            std::cerr << "Warning: " << missing << " node ids not found in graph mapping\n";
        }
        if (targets.size() < 2) throw std::runtime_error("need at least two valid nodes to route");

        std::ofstream out;
        if (args.out_binary) out.open(args.output_path, std::ios::binary);
        else out.open(args.output_path);
        if (!out) throw std::runtime_error("cannot open output file");

        std::vector<uint64_t>& dist = g.dist;
        std::vector<uint32_t>& pred = g.pred;
        std::vector<uint32_t> touched;
        touched.reserve(100000);

        std::vector<uint32_t> settled_epoch(g.N, 0);
        uint32_t epoch = 1;

        RadixHeap heap;
        std::vector<uint32_t> path;
        path.reserve(1024);

        bool warned_missing_osm = false;
        auto to_osm_id = [&](uint32_t gid) -> int64_t {
            if (gid < osm_by_graph.size() && osm_by_graph[gid] >= 0) return osm_by_graph[gid];
            if (!warned_missing_osm) {
                warned_missing_osm = true;
                std::cerr << "Warning: missing osm ids in mapping (fallback to graph ids)\n";
            }
            return static_cast<int64_t>(gid);
        };

        for (uint32_t s : targets) {
            heap.clear();
            touched.clear();
            const uint32_t remaining_init = static_cast<uint32_t>(targets.size() - 1);
            uint32_t remaining = remaining_init;

            dist[s] = 0;
            pred[s] = NO_NODE;
            touched.push_back(s);
            heap.push(0, s);

            while (!heap.empty() && remaining > 0) {
                const auto it = heap.pop();
                const uint64_t du = it.key;
                const uint32_t u = it.value;

                if (du != dist[u]) continue;
                if (settled_epoch[u] == epoch) continue;
                settled_epoch[u] = epoch;

                if (u != s && is_target[u]) {
                    if (remaining > 0) --remaining;
                    if (remaining == 0) break;
                }

                for (uint32_t p = g.offsets[u]; p < g.offsets[u + 1]; ++p) {
                    const uint32_t v = g.to[p];
                    const uint64_t nd = du + static_cast<uint64_t>(g.len_cm[p]);
                    if (nd < dist[v]) {
                        if (dist[v] == DIST_INF) touched.push_back(v);
                        dist[v] = nd;
                        pred[v] = u;
                        heap.push(nd, v);
                    }
                }
            }

            const int64_t source_osm = to_osm_id(s);
            for (uint32_t t : targets) {
                if (t == s) continue;
                const int64_t target_osm = to_osm_id(t);
                if (dist[t] == DIST_INF) {
                    if (args.out_binary) {
                        const uint64_t d = DIST_INF;
                        const uint32_t len = 0;
                        out.write(reinterpret_cast<const char*>(&source_osm), sizeof(int64_t));
                        out.write(reinterpret_cast<const char*>(&target_osm), sizeof(int64_t));
                        out.write(reinterpret_cast<const char*>(&d), sizeof(uint64_t));
                        out.write(reinterpret_cast<const char*>(&len), sizeof(uint32_t));
                    } else {
                        out << source_osm << ";" << target_osm << ";INF;\n";
                    }
                    continue;
                }

                if (!build_path(s, t, pred, path)) {
                    if (args.out_binary) {
                        const uint64_t d = DIST_INF;
                        const uint32_t len = 0;
                        out.write(reinterpret_cast<const char*>(&source_osm), sizeof(int64_t));
                        out.write(reinterpret_cast<const char*>(&target_osm), sizeof(int64_t));
                        out.write(reinterpret_cast<const char*>(&d), sizeof(uint64_t));
                        out.write(reinterpret_cast<const char*>(&len), sizeof(uint32_t));
                    } else {
                        out << source_osm << ";" << target_osm << ";INF;\n";
                    }
                    continue;
                }

                if (args.out_binary) {
                    const uint64_t d = dist[t];
                    const uint32_t len = static_cast<uint32_t>(path.size());
                    out.write(reinterpret_cast<const char*>(&source_osm), sizeof(int64_t));
                    out.write(reinterpret_cast<const char*>(&target_osm), sizeof(int64_t));
                    out.write(reinterpret_cast<const char*>(&d), sizeof(uint64_t));
                    out.write(reinterpret_cast<const char*>(&len), sizeof(uint32_t));
                    for (uint32_t gid : path) {
                        const int64_t node_osm = to_osm_id(gid);
                        out.write(reinterpret_cast<const char*>(&node_osm), sizeof(int64_t));
                    }
                } else {
                    out << source_osm << ";" << target_osm << ";" << dist[t] << ";";
                    for (std::size_t i = 0; i < path.size(); ++i) {
                        if (i > 0) out << ",";
                        out << to_osm_id(path[i]);
                    }
                    out << "\n";
                }
            }

            for (uint32_t v : touched) {
                dist[v] = DIST_INF;
                pred[v] = NO_NODE;
            }

            ++epoch;
            if (epoch == 0) {
                std::fill(settled_epoch.begin(), settled_epoch.end(), 0);
                epoch = 1;
            }
        }

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
}

