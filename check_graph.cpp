#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>

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
    if (!in) throw std::runtime_error("cannot open file");

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
        // Ancien format (sans magic) :
        // N, M, offsets, to, len_cm, cls, way_id, dist, label, pred
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

static void stats(const Graph& g) {
    std::cout << "Read header: N=" << g.N << " M=" << g.M << "\n";
    std::cout << "CSR OK\n";

    if (g.offsets.size() != static_cast<std::size_t>(g.N) + 1) {
        throw std::runtime_error("offsets size mismatch");
    }
    if (g.offsets.back() != g.M) {
        std::cerr << "Warning: offsets[N] != M (offsets[N]=" << g.offsets.back() << ")\n";
    }

    uint32_t min_deg = UINT32_MAX, max_deg = 0;
    uint64_t sum_deg = 0;

    for (uint32_t u = 0; u < g.N; ++u) {
        uint32_t deg = g.offsets[u+1] - g.offsets[u];
        min_deg = std::min(min_deg, deg);
        max_deg = std::max(max_deg, deg);
        sum_deg += deg;
    }

    double avg_deg = (g.N ? static_cast<double>(sum_deg) / static_cast<double>(g.N) : 0.0);
    std::cout << "Out-degree: min=" << min_deg << " max=" << max_deg << " avg=" << avg_deg << "\n";

    uint32_t min_len = UINT32_MAX, max_len = 0;
    uint64_t sum_len = 0;
    uint64_t zero_cnt = 0;

    for (uint32_t i = 0; i < g.M; ++i) {
        uint32_t L = g.len_cm[i];
        if (L == 0) ++zero_cnt;
        min_len = std::min(min_len, L);
        max_len = std::max(max_len, L);
        sum_len += L;
    }

    double avg_len_cm = (g.M ? static_cast<double>(sum_len) / static_cast<double>(g.M) : 0.0);
    std::cout << "Lengths (cm): min=" << min_len << " max=" << max_len
              << " avg=" << avg_len_cm << " zero_count=" << zero_cnt << "\n";
    std::cout << "Lengths (m):  min=" << (min_len/100.0) << " max=" << (max_len/100.0)
              << " avg=" << (avg_len_cm/100.0) << "\n";

    std::cout << "Restrictions (forbidden transitions): " << g.R << "\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: check_graph file.graph\n";
        return 1;
    }

    try {
        Graph g = read_graph(argv[1]);
        stats(g);
        std::cout << "OK\n";
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }

    return 0;
}
