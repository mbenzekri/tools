#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <stdexcept>
#include <string>

static void read_exact(std::ifstream& in, void* dst, std::size_t n) {
    in.read(reinterpret_cast<char*>(dst), static_cast<std::streamsize>(n));
    if (!in || static_cast<std::size_t>(in.gcount()) != n) {
        throw std::runtime_error("read_exact failed (file truncated or IO error)");
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: check_graph graph_file\n";
        return 1;
    }
    const std::string path = argv[1];

    std::ifstream in(path, std::ios::binary);
    if (!in) {
        std::cerr << "cannot open: " << path << "\n";
        return 1;
    }

    uint32_t N = 0, M = 0;
    read_exact(in, &N, 4);
    read_exact(in, &M, 4);

    std::cerr << "Read header: N=" << N << " M=" << M << "\n";

    std::vector<uint32_t> offsets(static_cast<size_t>(N) + 1);
    read_exact(in, offsets.data(), offsets.size() * sizeof(uint32_t));

    std::vector<uint32_t> to(M);
    std::vector<uint32_t> len_cm(M);
    std::vector<uint8_t>  cls(M);
    std::vector<int64_t>  wid(M);

    read_exact(in, to.data(),      to.size()      * sizeof(uint32_t));
    read_exact(in, len_cm.data(),  len_cm.size()  * sizeof(uint32_t));
    read_exact(in, cls.data(),     cls.size()     * sizeof(uint8_t));
    read_exact(in, wid.data(),     wid.size()     * sizeof(int64_t));

    // dist/label/pred existent dans le fichier : on les saute (mais on vérifie que ça existe bien)
    const uint64_t skip_dist_bytes  = static_cast<uint64_t>(N) * sizeof(uint64_t);
    const uint64_t skip_label_bytes = static_cast<uint64_t>(N) * sizeof(uint32_t);
    const uint64_t skip_pred_bytes  = static_cast<uint64_t>(N) * sizeof(uint32_t);

    in.seekg(static_cast<std::streamoff>(skip_dist_bytes + skip_label_bytes + skip_pred_bytes), std::ios::cur);
    if (!in) throw std::runtime_error("failed to skip dist/label/pred (file too small?)");

    // Optional: ensure EOF or additional bytes
    // in.peek();

    // =========================
    // 1) Vérifs CSR
    // =========================
    if (offsets[0] != 0) {
        std::cerr << "CSR ERROR: offsets[0] != 0 (" << offsets[0] << ")\n";
        return 2;
    }
    if (offsets[N] != M) {
        std::cerr << "CSR ERROR: offsets[N] != M (offsets[N]=" << offsets[N] << ")\n";
        return 2;
    }

    for (uint32_t i = 0; i < N; ++i) {
        if (offsets[i] > offsets[i+1]) {
            std::cerr << "CSR ERROR: offsets not monotonic at i=" << i
                      << " (" << offsets[i] << " > " << offsets[i+1] << ")\n";
            return 2;
        }
    }

    uint64_t bad_to = 0;
    for (uint32_t e = 0; e < M; ++e) {
        if (to[e] >= N) ++bad_to;
    }
    if (bad_to) {
        std::cerr << "CSR ERROR: " << bad_to << " arcs have to[e] >= N\n";
        return 2;
    }
    std::cerr << "CSR OK\n";

    // =========================
    // 2) Stats degré sortant
    // =========================
    uint32_t min_deg = std::numeric_limits<uint32_t>::max();
    uint32_t max_deg = 0;
    uint64_t sum_deg = 0;

    for (uint32_t u = 0; u < N; ++u) {
        uint32_t deg = offsets[u+1] - offsets[u];
        if (deg < min_deg) min_deg = deg;
        if (deg > max_deg) max_deg = deg;
        sum_deg += deg;
    }
    double avg_deg = (N > 0) ? static_cast<double>(sum_deg) / static_cast<double>(N) : 0.0;

    std::cerr << "Out-degree: min=" << min_deg << " max=" << max_deg << " avg=" << avg_deg << "\n";

    // =========================
    // 3) Stats longueurs
    // =========================
    uint32_t min_len = std::numeric_limits<uint32_t>::max();
    uint32_t max_len = 0;
    long double sum_len = 0.0L;
    uint64_t zero_len = 0;

    for (uint32_t e = 0; e < M; ++e) {
        uint32_t L = len_cm[e];
        if (L == 0) ++zero_len;
        if (L < min_len) min_len = L;
        if (L > max_len) max_len = L;
        sum_len += static_cast<long double>(L);
    }
    long double avg_len_cm = (M > 0) ? (sum_len / static_cast<long double>(M)) : 0.0L;

    std::cerr << "Lengths (cm): min=" << min_len << " max=" << max_len
              << " avg=" << static_cast<double>(avg_len_cm)
              << " zero_count=" << zero_len << "\n";

    std::cerr << "Lengths (m):  min=" << (min_len / 100.0)
              << " max=" << (max_len / 100.0)
              << " avg=" << (static_cast<double>(avg_len_cm) / 100.0) << "\n";

    // =========================
    // 4) Connectivité (non-orienté)
    //    On construit l'adjacence inverse (rev) pour parcourir en non-orienté.
    // =========================
    std::vector<uint32_t> indeg(N, 0);
    for (uint32_t u = 0; u < N; ++u) {
        for (uint32_t p = offsets[u]; p < offsets[u+1]; ++p) {
            uint32_t v = to[p];
            ++indeg[v];
        }
    }

    std::vector<uint32_t> rev_off(N + 1, 0);
    for (uint32_t i = 0; i < N; ++i) rev_off[i+1] = rev_off[i] + indeg[i];

    std::vector<uint32_t> rev_to(M);
    auto cur = rev_off;
    for (uint32_t u = 0; u < N; ++u) {
        for (uint32_t p = offsets[u]; p < offsets[u+1]; ++p) {
            uint32_t v = to[p];
            uint32_t pos = cur[v]++;
            rev_to[pos] = u;
        }
    }

    std::vector<uint8_t> visited(N, 0);
    uint32_t components = 0;
    uint32_t largest = 0;

    std::queue<uint32_t> q;
    for (uint32_t s = 0; s < N; ++s) {
        if (visited[s]) continue;
        ++components;
        uint32_t size = 0;
        visited[s] = 1;
        q.push(s);

        while (!q.empty()) {
            uint32_t u = q.front();
            q.pop();
            ++size;

            // voisins sortants
            for (uint32_t p = offsets[u]; p < offsets[u+1]; ++p) {
                uint32_t v = to[p];
                if (!visited[v]) {
                    visited[v] = 1;
                    q.push(v);
                }
            }
            // voisins entrants (pour non-orienté)
            for (uint32_t p = rev_off[u]; p < rev_off[u+1]; ++p) {
                uint32_t v = rev_to[p];
                if (!visited[v]) {
                    visited[v] = 1;
                    q.push(v);
                }
            }
        }

        if (size > largest) largest = size;
    }

    double largest_pct = (N > 0) ? (100.0 * static_cast<double>(largest) / static_cast<double>(N)) : 0.0;
    std::cerr << "Connectivity (undirected view): components=" << components
              << " largest_component=" << largest << " (" << largest_pct << "%)\n";

    std::cerr << "OK\n";
    return 0;
}

