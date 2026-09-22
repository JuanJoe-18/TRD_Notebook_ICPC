/**
 * @file GraphTemplateICPC.cpp
 * @brief Plantilla de algoritmos y estructuras para grafos.
 * @details Incluye DSU, BFS 0-1, Dijkstra (con camino), SCC, ciclos, Euler, LCA, 2-SAT, Flujos y utilidades de Grilla.
 * @note Ajusta la indexacion y elimina las estructuras que no uses.
 */
//   ____ ___  ____  _____   ____  _   _
//  / ___/ _ \|  _ \| ____| / ___|| | | |
// | |  | | | | | | |  _|   \___ \| | | |
// | |__| |_| | |_| | |___   ___) | |_| |
//  \____\___/|____/|_____| |____/ \___/
//
//                    GRAPH TEMPLATE

#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
const ll LINF = 1e18;

// ==========================================
// UTILIDADES PARA GRILLAS (Grid / Matrices)
// ==========================================
const int dx4[4] = {1, -1, 0, 0};
const int dy4[4] = {0, 0, 1, -1};
const int dx8[8] = {1, 1, 0, -1, -1, -1, 0, 1};
const int dy8[8] = {0, 1, 1, 1, 0, -1, -1, -1};

// Helper para validar límites de grilla en O(1)
#define isValid(x, y, R, C) ((x) >= 0 && (x) < (R) && (y) >= 0 && (y) < (C))


// ==========================================
// ESTRUCTURA MAESTRA DE GRAFOS (ICPC)
// ==========================================
template <typename T = ll>
struct Edge {
    int from, to;
    T weight;
    int id; // Util para caminos Eulerianos o puentes
};

// ==========================================
// DISJOINT SET UNION (DSU)
// ==========================================
struct DSU {
    vector<int> p, sz;
    vector<vector<int>> members;
    int num_components;

    DSU(int n) : p(n), sz(n, 1), members(n), num_components(n) {
        iota(p.begin(), p.end(), 0);
        for (int i = 0; i < n; i++) members[i].push_back(i);
    }

    int find(int x) {
        return p[x] == x ? x : p[x] = find(p[x]);
    }

    bool unite(int a, int b) {
        a = find(a);
        b = find(b);
        if (a == b) return false;
        if (sz[a] < sz[b]) swap(a, b);
        p[b] = a;
        sz[a] += sz[b];
        members[a].reserve(members[a].size() + members[b].size());
        members[a].insert(members[a].end(), members[b].begin(), members[b].end());
        vector<int>().swap(members[b]);
        num_components--;
        return true;
    }

    const vector<int>& get_component_nodes(int x) { return members[find(x)]; }
    int component_size(int x) { return sz[find(x)]; }
    
    vector<vector<int>> get_components(int first_node = 0) const {
        vector<vector<int>> result;
        for (int root = first_node; root < (int)p.size(); root++) {
            if (p[root] == root && !members[root].empty()) {
                result.push_back(members[root]);
            }
        }
        return result;
    }

    int component_count(int first_node = 0) const {
        if (first_node == 0) return num_components;
        int count = 0;
        for (int root = first_node; root < (int)p.size(); root++) {
            if (p[root] == root && !members[root].empty()) count++;
        }
        return count;
    }
};

template <typename T = ll>
struct Graph {
    int n;
    vector<vector<int>> adj; 
    vector<Edge<T>> edges;  

    Graph(int _n) : n(_n) {
        adj.resize(n + 1);
    }

    void add_directed_edge(int from, int to, T weight = 1, int id = -1) {
        adj[from].push_back(edges.size());
        edges.push_back({from, to, weight, id});
    }

    void add_undirected_edge(int u, int v, T weight = 1, int id = -1) {
        add_directed_edge(u, v, weight, id);
        add_directed_edge(v, u, weight, id);
    }

    // ---------------------------------------------------------
    // ALGORITMOS BÁSICOS INTEGRADOS
    // ---------------------------------------------------------

    // 1. Dijkstra O(E log V) Estándar
    vector<T> dijkstra(int src) {
        vector<T> dist(n + 1, LINF);
        priority_queue<pair<T, int>, vector<pair<T, int>>, greater<>> pq;
        dist[src] = 0;
        pq.push({0, src});

        while (!pq.empty()) {
            auto [d, u] = pq.top(); pq.pop();
            if (d > dist[u]) continue;
            for (int id : adj[u]) {
                auto& edge = edges[id];
                if (dist[u] + edge.weight < dist[edge.to]) {
                    dist[edge.to] = dist[u] + edge.weight;
                    pq.push({dist[edge.to], edge.to});
                }
            }
        }
        return dist;
    }

    // 2. Dijkstra O(E log V) con Reconstrucción de Camino
    // Retorna {arreglo_distancias, arreglo_padres}
    pair<vector<T>, vector<int>> dijkstra_path(int src) {
        vector<T> dist(n + 1, LINF);
        vector<int> p(n + 1, -1);
        priority_queue<pair<T, int>, vector<pair<T, int>>, greater<>> pq;
        dist[src] = 0;
        pq.push({0, src});

        while (!pq.empty()) {
            auto [d, u] = pq.top(); pq.pop();
            if (d > dist[u]) continue;
            for (int id : adj[u]) {
                auto& edge = edges[id];
                if (dist[u] + edge.weight < dist[edge.to]) {
                    dist[edge.to] = dist[u] + edge.weight;
                    p[edge.to] = u; // Guardar el nodo desde donde llegamos al óptimo
                    pq.push({dist[edge.to], edge.to});
                }
            }
        }
        return {dist, p};
    }

    // Método de utilidad para reconstruir el camino tras llamar a dijkstra_path
    vector<int> restore_path(int target, const vector<int>& p) {
        vector<int> path;
        for (int v = target; v != -1; v = p[v]) {
            path.push_back(v);
        }
        reverse(path.begin(), path.end());
        return path;
    }

    // 3. BFS 0-1 en O(V + E)
    // Para grafos donde los pesos son estrictamente 0 o 1. No requiere Priority Queue.
    vector<T> zero_one_bfs(int src) {
        vector<T> dist(n + 1, LINF); 
        deque<int> dq;
        dist[src] = 0;
        dq.push_back(src);

        while (!dq.empty()) {
            int u = dq.front(); 
            dq.pop_front();

            for (int id : adj[u]) {
                auto& edge = edges[id];
                if (dist[u] + edge.weight < dist[edge.to]) {
                    dist[edge.to] = dist[u] + edge.weight;
                    if (edge.weight == 1) {
                        dq.push_back(edge.to);  // Costo 1: Atrás de la fila
                    } else {
                        dq.push_front(edge.to); // Costo 0: Adelante de la fila
                    }
                }
            }
        }
        return dist;
    }

    // 4. Verificación de Grafo Bipartito (2-Coloring)
    // Retorna {es_bipartito, arreglo_de_colores}. Los colores son 0 o 1.
    pair<bool, vector<int>> is_bipartite() {
        vector<int> color(n + 1, -1);
        bool bipartite = true;

        for (int i = 1; i <= n; i++) {
            if (color[i] != -1) continue;
            queue<int> q;
            q.push(i);
            color[i] = 0;

            while (!q.empty()) {
                int u = q.front(); q.pop();
                for (int id : adj[u]) {
                    int v = edges[id].to;
                    if (color[v] == -1) {
                        color[v] = color[u] ^ 1;
                        q.push(v);
                    } else if (color[v] == color[u]) {
                        bipartite = false;
                    }
                }
            }
        }
        return {bipartite, color};
    }

    // 5. Floyd-Warshall O(V^3) - All-Pairs Shortest Path
    vector<vector<T>> floyd_warshall() {
        vector<vector<T>> dist(n + 1, vector<T>(n + 1, LINF));
        for (int i = 1; i <= n; i++) dist[i][i] = 0;
        for (auto& e : edges) dist[e.from][e.to] = min(dist[e.from][e.to], e.weight);
        for (int k = 1; k <= n; k++) {
            for (int i = 1; i <= n; i++) {
                for (int j = 1; j <= n; j++) {
                    if (dist[i][k] < LINF && dist[k][j] < LINF) {
                        dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
                    }
                }
            }
        }
        return dist;
    }

    // 6. Ordenamiento Topológico (Kahn's Algorithm)
    vector<int> topo_sort() {
        vector<int> in_degree(n + 1, 0), order;
        for (auto& e : edges) in_degree[e.to]++;
        queue<int> q;
        for (int i = 1; i <= n; i++) if (in_degree[i] == 0) q.push(i);
        while (!q.empty()) {
            int u = q.front(); q.pop();
            order.push_back(u);
            for (int id : adj[u]) {
                if (--in_degree[edges[id].to] == 0) q.push(edges[id].to);
            }
        }
        return (order.size() == n) ? order : vector<int>(); // Vacio si hay ciclo
    }

    // 7. Tarjan para Componentes Fuertemente Conexas (SCC)
    vector<int> get_scc() {
        vector<int> val(n + 1, 0), comp(n + 1, -1);
        vector<int> z; int timer = 0, ncomps = 0;

        auto dfs = [&](auto& self, int u) -> int {
            int low = val[u] = ++timer, x;
            z.push_back(u);
            for (int id : adj[u]) {
                int v = edges[id].to;
                if (comp[v] < 0) low = min(low, val[v] ?: self(self, v));
            }
            if (low == val[u]) {
                do {
                    x = z.back(); z.pop_back();
                    comp[x] = ncomps;
                } while (x != u);
                ncomps++;
            }
            return val[u] = low;
        };
        for (int i = 1; i <= n; i++) if (!val[i]) dfs(dfs, i);
        return comp;
    }

    // 8. Bellman-Ford O(V E)
    tuple<bool, vector<T>, vector<int>, vector<int>> bellman_ford(int src) {
        vector<T> dist(n + 1, LINF);
        vector<int> p(n + 1, -1);
        dist[src] = 0;
        int x = -1;
        for (int i = 0; i < n; i++) {
            x = -1;
            for (auto& edge : edges) {
                if (dist[edge.from] < LINF) {
                    if (dist[edge.to] > dist[edge.from] + edge.weight) {
                        dist[edge.to] = max(-LINF, dist[edge.from] + edge.weight);
                        p[edge.to] = edge.from;
                        x = edge.to;
                    }
                }
            }
        }
        vector<int> cycle;
        bool ok = true;
        if (x != -1) {
            ok = false;
            int y = x;
            for (int i = 0; i < n; i++) y = p[y];
            for (int cur = y;; cur = p[cur]) {
                cycle.push_back(cur);
                if (cur == y && cycle.size() > 1) break;
            }
            reverse(cycle.begin(), cycle.end());
        }
        return {ok, dist, p, cycle};
    }

    // 9. Encontrar Ciclo Dirigido (Round Trip II)
    vector<int> find_directed_cycle() {
        vector<int> state(n + 1, 0); 
        vector<int> parent(n + 1, -1);
        vector<int> cycle;
        auto dfs = [&](auto& self, int u) -> bool {
            state[u] = 1; 
            for (int id : adj[u]) {
                int v = edges[id].to;
                if (state[v] == 0) {
                    parent[v] = u;
                    if (self(self, v)) return true;
                } else if (state[v] == 1) { 
                    cycle.push_back(v);
                    int curr = u;
                    while (curr != -1 && curr != v) { 
                        cycle.push_back(curr);
                        curr = parent[curr];
                    }
                    cycle.push_back(v);
                    reverse(cycle.begin(), cycle.end());
                    return true;
                }
            }
            state[u] = 2; 
            return false;
        };
        for (int i = 1; i <= n; i++) if (state[i] == 0) if (dfs(dfs, i)) return cycle;
        return cycle;
    }

    // 10. Encontrar Ciclo No Dirigido (Round Trip I)
    vector<int> find_undirected_cycle() {
        vector<int> state(n + 1, 0); 
        vector<int> parent_node(n + 1, -1);
        vector<int> cycle;
        auto dfs = [&](auto& self, int u, int p_edge_idx) -> bool {
            state[u] = 1;
            for (int id : adj[u]) {
                if ((id ^ 1) == p_edge_idx) continue;
                int v = edges[id].to;
                if (state[v] == 0) {
                    parent_node[v] = u;
                    if (self(self, v, id)) return true;
                } else if (state[v] == 1) { 
                    cycle.push_back(v);
                    int curr = u;
                    while (curr != -1 && curr != v) { 
                        cycle.push_back(curr);
                        curr = parent_node[curr];
                    }
                    cycle.push_back(v);
                    reverse(cycle.begin(), cycle.end());
                    return true;
                }
            }
            state[u] = 2;
            return false;
        };
        for (int i = 1; i <= n; i++) if (state[i] == 0) if (dfs(dfs, i, -1)) return cycle;
        return cycle;
    }

    // 11. Kruskal MST
    pair<T, vector<int>> kruskal() {
        vector<int> edge_indices(edges.size());
        iota(edge_indices.begin(), edge_indices.end(), 0);
        sort(edge_indices.begin(), edge_indices.end(), [&](int a, int b) {
            return edges[a].weight < edges[b].weight;
        });
        DSU dsu(n + 1);
        T mst_weight = 0;
        vector<int> mst_edges;
        for (int idx : edge_indices) {
            auto& e = edges[idx];
            if (dsu.unite(e.from, e.to)) {
                mst_weight += e.weight;
                mst_edges.push_back(idx);
                if (mst_edges.size() == (size_t)(n - 1)) break;
            }
        }
        return {mst_weight, mst_edges};
    }

    // 12. Circuito Euleriano
    vector<int> eulerian_circuit(int start_node = 1, bool undirected = true) {
        vector<int> in_deg(n + 1, 0), out_deg(n + 1, 0);
        int max_id = -1;
        for (auto& e : edges) {
            out_deg[e.from]++;
            in_deg[e.to]++;
            max_id = max(max_id, e.id);
        }
        for (int i = 1; i <= n; i++) {
            if (undirected) {
                if (out_deg[i] % 2 != 0) return {};
            } else {
                if (in_deg[i] != out_deg[i]) return {};
            }
        }
        vector<bool> used_edge(max_id + 1, false);
        vector<int> circuit;
        vector<int> head(n + 1, 0); 
        auto dfs = [&](auto& self, int u) -> void {
            while (head[u] < adj[u].size()) {
                int edge_idx = adj[u][head[u]++];
                auto& e = edges[edge_idx];
                if (e.id != -1) {
                    if (!used_edge[e.id]) {
                        used_edge[e.id] = true;
                        self(self, e.to);
                    }
                } else {
                    self(self, e.to);
                }
            }
            circuit.push_back(u);
        };
        dfs(dfs, start_node);
        reverse(circuit.begin(), circuit.end());
        int required_size = undirected ? (edges.size() / 2 + 1) : (edges.size() + 1);
        if (circuit.size() != required_size && edges.size() > 0) return {};
        return circuit;
    }
};

// ==========================================
// GRAFOS FUNCIONALES (Successor Graphs)
// ==========================================
struct FunctionalGraph {
    int n, log_k;
    vector<vector<int>> up;
    vector<int> succ, in_cycle, cycle_id, cycle_pos, cycle_size, dist_to_cycle;
    vector<vector<int>> cycles;

    FunctionalGraph(int _n, int max_k_bits = 60) : n(_n), log_k(max_k_bits) {
        up.assign(n + 1, vector<int>(log_k, 0));
        succ.assign(n + 1, 0);
        in_cycle.assign(n + 1, 0);
        cycle_id.assign(n + 1, -1);
        cycle_pos.assign(n + 1, -1);
        cycle_size.assign(n + 1, 0);
        dist_to_cycle.assign(n + 1, 0);
    }

    void decompose_cycles(const vector<int>& _succ) {
        succ = _succ;
        vector<int> state(n + 1, 0); 
        cycles.clear();
        for (int i = 1; i <= n; i++) {
            if (state[i] != 0) continue;
            int curr = i;
            vector<int> path;
            while (curr >= 1 && curr <= n && state[curr] == 0) {
                state[curr] = 1;
                path.push_back(curr);
                curr = succ[curr];
            }
            if (curr >= 1 && curr <= n && state[curr] == 1) {
                vector<int> cyc;
                bool found_start = false;
                for (int u : path) {
                    if (u == curr) found_start = true;
                    if (found_start) cyc.push_back(u);
                }
                int cid = cycles.size(), c_sz = cyc.size();
                for (int pos = 0; pos < c_sz; pos++) {
                    int u = cyc[pos];
                    in_cycle[u] = 1; cycle_id[u] = cid;
                    cycle_pos[u] = pos; cycle_size[u] = c_sz;
                    dist_to_cycle[u] = 0;
                }
                cycles.push_back(cyc);
            }
            for (int u : path) state[u] = 2;
        }
        for (int i = 1; i <= n; i++) {
            if (in_cycle[i]) continue;
            int curr = i;
            vector<int> path;
            while (curr >= 1 && curr <= n && !in_cycle[curr] && cycle_id[curr] == -1) {
                path.push_back(curr);
                curr = succ[curr];
            }
            int target_cid = (curr >= 1 && curr <= n) ? cycle_id[curr] : -1;
            int target_csz = (curr >= 1 && curr <= n) ? cycle_size[curr] : 0;
            int d = (curr >= 1 && curr <= n) ? dist_to_cycle[curr] : 0;
            for (int j = (int)path.size() - 1; j >= 0; j--) {
                d++; int u = path[j];
                cycle_id[u] = target_cid; cycle_size[u] = target_csz; dist_to_cycle[u] = d;
            }
        }
    }

    void build_binary_lifting(const vector<int>& _succ) {
        succ = _succ;
        for (int i = 1; i <= n; i++) up[i][0] = succ[i];
        for (int j = 1; j < log_k; j++) {
            for (int i = 1; i <= n; i++) {
                int p = up[i][j - 1];
                up[i][j] = (p >= 1 && p <= n) ? up[p][j - 1] : 0;
            }
        }
    }

    int get_kth_successor(int u, ll k) {
        for (int j = 0; j < log_k; j++) {
            if (k & (1LL << j)) {
                u = up[u][j];
                if (u < 1 || u > n) return 0;
            }
        }
        return u;
    }
};

// LOWEST COMMON ANCESTOR (Árboles)
// ==========================================
struct LCA {
    int n, log_n;
    vector<vector<int>> up;
    vector<int> depth;

    LCA(int _n, int _log_n = 20) : n(_n), log_n(_log_n) {
        up.assign(n + 1, vector<int>(log_n, 0));
        depth.assign(n + 1, 0);
    }

    void dfs(int u, int p, const vector<vector<int>>& adj) {
        up[u][0] = p;
        for (int j = 1; j < log_n; j++) up[u][j] = up[up[u][j - 1]][j - 1];
        for (int v : adj[u]) {
            if (v != p) {
                depth[v] = depth[u] + 1;
                dfs(v, u, adj);
            }
        }
    }

    void build(int root, const vector<vector<int>>& adj) { dfs(root, root, adj); }

    int get_lca(int u, int v) {
        if (depth[u] < depth[v]) swap(u, v);
        int k = depth[u] - depth[v];
        for (int j = 0; j < log_n; j++) if (k & (1 << j)) u = up[u][j];
        if (u == v) return u;
        for (int j = log_n - 1; j >= 0; j--) {
            if (up[u][j] != up[v][j]) {
                u = up[u][j]; v = up[v][j];
            }
        }
        return up[u][0];
    }

    int get_dist(int u, int v) {
        return depth[u] + depth[v] - 2 * depth[get_lca(u, v)];
    }
};

// 2-SAT (Satisfactibilidad Booleana)
// ==========================================
struct TwoSat {
    int n;
    Graph<ll> G;
    vector<bool> assignment;

    TwoSat(int _n) : n(_n), G(2 * _n) { assignment.assign(n + 1, false); }
    int get_node(int u, bool is_true) { return is_true ? u : u + n; }

    void add_clause(int u, bool is_u_true, int v, bool is_v_true) {
        G.add_directed_edge(get_node(u, !is_u_true), get_node(v, is_v_true)); 
        G.add_directed_edge(get_node(v, !is_v_true), get_node(u, is_u_true)); 
    }
    void force_value(int u, bool is_true) { add_clause(u, is_true, u, is_true); }
    void add_implication(int u, bool is_u_true, int v, bool is_v_true) { add_clause(u, !is_u_true, v, is_v_true); }
    void add_xor(int u, bool is_u_true, int v, bool is_v_true) {
        add_clause(u, is_u_true, v, is_v_true);
        add_clause(u, !is_u_true, v, !is_v_true);
    }
    void add_equivalence(int u, bool is_u_true, int v, bool is_v_true) {
        add_implication(u, is_u_true, v, is_v_true);
        add_implication(v, is_v_true, u, is_u_true);
    }

    bool solve() {
        vector<int> comp = G.get_scc();
        for (int i = 1; i <= n; i++) {
            if (comp[i] == comp[i + n]) return false;
            assignment[i] = comp[i] < comp[i + n];
        }
        return true;
    }
};

// ==========================================
// EMPAREJAMIENTO BIPARTITO (Kuhn's Algorithm) O(V * E)
// ==========================================
struct BipartiteMatcher {
    int n, m;
    vector<vector<int>> adj;
    vector<int> match; 
    vector<bool> vis;

    BipartiteMatcher(int _n, int _m) : n(_n), m(_m) {
        adj.resize(n + 1);
        match.assign(m + 1, -1);
    }
    void add_edge(int u, int v) { adj[u].push_back(v); }
    bool dfs(int u) {
        for (int v : adj[u]) {
            if (vis[v]) continue;
            vis[v] = true;
            if (match[v] < 0 || dfs(match[v])) {
                match[v] = u; return true;
            }
        }
        return false;
    }
    int solve() {
        int ans = 0;
        for (int i = 1; i <= n; i++) {
            vis.assign(m + 1, false);
            if (dfs(i)) ans++;
        }
        return ans;
    }
};

// ==========================================
// EMPAREJAMIENTO BIPARTITO (Kuhn's Algorithm) O()
// ==========================================
struct HopcroftKarp {
    int n, m;
    vector<vector<int>> adj;
    vector<int> pairU, pairV, dist;
    const int INF = 1e9;

    HopcroftKarp(int n, int m) : n(n), m(m), adj(n + 1), pairU(n + 1, 0), pairV(m + 1, 0), dist(n + 1) {}

    void add_edge(int u, int v) {
        adj[u].push_back(v); // u en [1..n], v en [1..m]
    }

    bool bfs() {
        queue<int> q;
        for (int u = 1; u <= n; u++) {
            if (pairU[u] == 0) { dist[u] = 0; q.push(u); } 
            else { dist[u] = INF; }
        }
        dist[0] = INF;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            if (dist[u] < dist[0]) {
                for (int v : adj[u]) {
                    if (dist[pairV[v]] == INF) {
                        dist[pairV[v]] = dist[u] + 1;
                        q.push(pairV[v]);
                    }
                }
            }
        }
        return dist[0] != INF;
    }

    bool dfs(int u) {
        if (u != 0) {
            for (int v : adj[u]) {
                if (dist[pairV[v]] == dist[u] + 1 && dfs(pairV[v])) {
                    pairV[v] = u;
                    pairU[u] = v;
                    return true;
                }
            }
            dist[u] = INF;
            return false;
        }
        return true;
    }

    int max_matching() {
        int res = 0;
        while (bfs()) {
            for (int u = 1; u <= n; u++) {
                if (pairU[u] == 0 && dfs(u)) res++;
            }
        }
        return res;
    }
};

// ==========================================
// PUENTES, PUNTOS DE ARTICULACIÓN Y BRIDGE TREE
// ==========================================
struct BridgeTree {
    int n, timer, num_comps;
    vector<vector<pair<int, int>>> adj; 
    vector<int> tin, low, comp_id;
    vector<bool> is_bridge, is_articulation;
    vector<vector<int>> tree_adj;

    BridgeTree(int _n, int num_edges) : n(_n) {
        adj.resize(n + 1);
        tin.assign(n + 1, -1); low.assign(n + 1, -1); comp_id.assign(n + 1, 0);
        is_bridge.assign(num_edges, false); is_articulation.assign(n + 1, false);
        timer = 0; num_comps = 0;
    }
    void add_edge(int u, int v, int id) {
        adj[u].push_back({v, id}); adj[v].push_back({u, id});
    }
    void dfs_tarjan(int u, int p = -1) {
        tin[u] = low[u] = ++timer;
        int children = 0;
        for (auto& edge : adj[u]) {
            int v = edge.first, id = edge.second;
            if (v == p) continue;
            if (tin[v] != -1) {
                low[u] = min(low[u], tin[v]);
            } else {
                children++;
                dfs_tarjan(v, u);
                low[u] = min(low[u], low[v]);
                if (low[v] > tin[u]) is_bridge[id] = true;
                if (low[v] >= tin[u] && p != -1) is_articulation[u] = true;
            }
        }
        if (p == -1 && children > 1) is_articulation[u] = true;
    }
    void dfs_comp(int u, int current_comp) {
        comp_id[u] = current_comp;
        for (auto& edge : adj[u]) {
            int v = edge.first, id = edge.second;
            if (comp_id[v] == 0) {
                if (is_bridge[id]) {
                    num_comps++;
                    tree_adj.push_back(vector<int>());
                    tree_adj[current_comp].push_back(num_comps);
                    tree_adj[num_comps].push_back(current_comp);
                    dfs_comp(v, num_comps);
                } else {
                    dfs_comp(v, current_comp);
                }
            }
        }
    }
    void build() {
        for (int i = 1; i <= n; i++) if (tin[i] == -1) dfs_tarjan(i);
        tree_adj.push_back(vector<int>());
        for (int i = 1; i <= n; i++) {
            if (comp_id[i] == 0) {
                num_comps++;
                tree_adj.push_back(vector<int>());
                dfs_comp(i, num_comps);
            }
        }
    }
};

// ==========================================
// MAX FLOW (Dinic's Algorithm) O(V^2 * E)
// ==========================================
struct Dinic {
    struct FlowEdge {
        int v, u;
        ll cap, flow = 0;
        FlowEdge(int _u, int _v, ll _cap) : u(_u), v(_v), cap(_cap) {}
    };
    int n, s, t;
    vector<FlowEdge> edges;
    vector<vector<int>> adj;
    vector<int> level, ptr;

    Dinic(int _n, int _s, int _t) : n(_n), s(_s), t(_t) {
        adj.resize(n + 1); level.resize(n + 1); ptr.resize(n + 1);
    }
    void add_edge(int u, int v, ll cap, bool directed = true) {
        adj[u].push_back(edges.size()); edges.push_back(FlowEdge(u, v, cap));
        adj[v].push_back(edges.size()); edges.push_back(FlowEdge(v, u, directed ? 0 : cap));
    }
    bool bfs() {
        fill(level.begin(), level.end(), -1);
        level[s] = 0; queue<int> q; q.push(s);
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int id : adj[u]) {
                if (edges[id].cap - edges[id].flow < 1) continue;
                if (level[edges[id].v] != -1) continue;
                level[edges[id].v] = level[u] + 1;
                q.push(edges[id].v);
            }
        }
        return level[t] != -1;
    }
    ll dfs(int u, ll pushed) {
        if (pushed == 0) return 0;
        if (u == t) return pushed;
        for (int& cid = ptr[u]; cid < adj[u].size(); ++cid) {
            int id = adj[u][cid];
            int v = edges[id].v;
            ll tr = edges[id].cap - edges[id].flow;
            if (level[u] + 1 != level[v] || tr == 0) continue;
            ll push = dfs(v, min(pushed, tr));
            if (push == 0) continue;
            edges[id].flow += push;
            edges[id ^ 1].flow -= push;
            return push;
        }
        return 0;
    }
    ll max_flow() {
        ll flow = 0;
        while (bfs()) {
            fill(ptr.begin(), ptr.end(), 0);
            while (ll pushed = dfs(s, LINF)) flow += pushed;
        }
        return flow;
    }
};

// ==========================================
// MIN COST MAX FLOW (Dijkstra with Potentials)
// ==========================================
struct MCMF {
    struct Edge {
        int to;
        ll cap, flow, cost;
        int rev;
    };
    int n;
    vector<vector<Edge>> adj;
    vector<ll> dist, pot;
    vector<int> p_node, p_edge;
    const ll LINF = 1e18; // Asegúrate de tener esta constante

    MCMF(int _n) : n(_n) {
        adj.resize(n + 1); dist.resize(n + 1); pot.resize(n + 1, 0);
        p_node.resize(n + 1); p_edge.resize(n + 1);
    }

    void add_edge(int u, int v, ll cap, ll cost) {
        adj[u].push_back({v, cap, 0, cost, (int)adj[v].size()});
        adj[v].push_back({u, 0, 0, -cost, (int)adj[u].size() - 1});
    }

    // Inicializa potenciales. Solo necesario si el grafo inicial tiene aristas negativas.
    // Si sabes que todos los costos iniciales son >= 0, puedes comentar la llamada en solve().
    void init_potentials(int s) {
        fill(pot.begin(), pot.end(), LINF);
        pot[s] = 0;
        bool changed = true;
        for (int i = 0; i < n && changed; i++) {
            changed = false;
            for (int u = 0; u <= n; u++) {
                if (pot[u] == LINF) continue;
                for (auto& e : adj[u]) {
                    if (e.cap > 0 && pot[e.to] > pot[u] + e.cost) {
                        pot[e.to] = pot[u] + e.cost;
                        changed = true;
                    }
                }
            }
        }
    }

    bool dijkstra(int s, int t) {
        fill(dist.begin(), dist.end(), LINF);
        priority_queue<pair<ll, int>, vector<pair<ll, int>>, greater<pair<ll, int>>> pq;
        
        dist[s] = 0;
        pq.push({0, s});
        
        while (!pq.empty()) {
            auto [d, u] = pq.top(); pq.pop();
            
            if (d > dist[u]) continue;
            
            for (int i = 0; i < adj[u].size(); i++) {
                auto& e = adj[u][i];
                // Costo reducido con potenciales de Johnson
                ll reduced_cost = e.cost + pot[u] - pot[e.to];
                
                if (e.cap - e.flow > 0 && dist[e.to] > dist[u] + reduced_cost) {
                    dist[e.to] = dist[u] + reduced_cost;
                    p_node[e.to] = u; p_edge[e.to] = i;
                    pq.push({dist[e.to], e.to});
                }
            }
        }
        return dist[t] != LINF;
    }

    pair<ll, ll> solve(int s, int t) {
        init_potentials(s); 
        ll flow = 0, cost = 0;
        
        while (dijkstra(s, t)) {
            // Actualización del potencial de los nodos
            for (int i = 0; i <= n; i++) {
                if (dist[i] != LINF) pot[i] += dist[i];
            }
            
            ll push = LINF; 
            int curr = t;
            while (curr != s) {
                int p = p_node[curr], idx = p_edge[curr];
                push = min(push, adj[p][idx].cap - adj[p][idx].flow);
                curr = p;
            }
            
            flow += push; 
            curr = t;
            while (curr != s) {
                int p = p_node[curr], idx = p_edge[curr], rev_idx = adj[p][idx].rev;
                adj[p][idx].flow += push; 
                adj[curr][rev_idx].flow -= push;
                cost += push * adj[p][idx].cost; 
                curr = p;
            }
        }
        return {flow, cost};
    }
};

struct DominatorTree {
    int n, t;
    vector<vector<int>> adj, rev, dom;
    vector<int> dfn, id, p, sdom, idom, dsu, best;

    DominatorTree(int n) : n(n), t(0), adj(n + 1), rev(n + 1), dom(n + 1), 
                           dfn(n + 1, 0), id(n + 1), p(n + 1), sdom(n + 1), 
                           idom(n + 1), dsu(n + 1), best(n + 1) {}

    void add_edge(int u, int v) { adj[u].push_back(v); }

    void dfs(int u) {
        dfn[u] = ++t; id[t] = u;
        sdom[t] = best[t] = dsu[t] = t;
        for (int v : adj[u]) {
            if (!dfn[v]) { dfs(v); p[dfn[v]] = dfn[u]; }
            rev[dfn[v]].push_back(dfn[u]);
        }
    }

    int find(int x) {
        if (x == dsu[x]) return x;
        int y = find(dsu[x]);
        if (sdom[best[x]] > sdom[best[dsu[x]]]) best[x] = best[dsu[x]];
        return dsu[x] = y;
    }

    void build(int root) {
        dfs(root);
        vector<vector<int>> bkt(n + 1);
        for (int i = t; i >= 2; i--) {
            for (int u : rev[i]) {
                find(u);
                if (sdom[best[u]] < sdom[i]) sdom[i] = sdom[best[u]];
            }
            bkt[sdom[i]].push_back(i);
            dsu[i] = p[i];
            for (int u : bkt[p[i]]) {
                find(u);
                idom[u] = (sdom[best[u]] == sdom[u]) ? sdom[u] : best[u];
            }
            bkt[p[i]].clear();
        }
        for (int i = 2; i <= t; i++) {
            if (idom[i] != sdom[i]) idom[i] = idom[idom[i]];
            dom[id[idom[i]]].push_back(id[i]);
        }
    }
};


int main() {
    ios_base::sync_with_stdio(0); cin.tie(0);

    // ==========================================
    // EJEMPLOS DE USO NUEVOS
    // ==========================================

    /* 
    Graph<ll> G(5);
    G.add_directed_edge(1, 2, 10);
    G.add_directed_edge(2, 3, 5);
    G.add_directed_edge(1, 4, 2);
    G.add_directed_edge(4, 3, 20);

    // 1. Dijkstra con Reconstrucción
    auto [distancias, padres] = G.dijkstra_path(1);
    cout << "Distancia a 3: " << distancias[3] << "\n"; // Ruta óptima: 1 -> 2 -> 3 (Costo 15)
    
    vector<int> camino = G.restore_path(3, padres);
    cout << "Camino a 3: ";
    for (int nodo : camino) cout << nodo << " "; 
    cout << "\n\n";

    // 2. BFS 0-1 (Solo util para pesos de 0 o 1)
    Graph<ll> G01(3);
    G01.add_directed_edge(1, 2, 0); // Costo 0
    G01.add_directed_edge(2, 3, 1); // Costo 1
    vector<ll> dist01 = G01.zero_one_bfs(1);
    cout << "Distancia BFS 0-1 al nodo 3: " << dist01[3] << "\n\n";

    // 3. Verificación de Grafo Bipartito
    Graph<ll> G_Bip(4);
    G_Bip.add_undirected_edge(1, 2);
    G_Bip.add_undirected_edge(2, 3);
    G_Bip.add_undirected_edge(3, 4);
    G_Bip.add_undirected_edge(4, 1);
    auto [es_bipartito, colores] = G_Bip.is_bipartite();
    cout << "¿Es Bipartito? " << (es_bipartito ? "Si" : "No") << "\n";
    if (es_bipartito) {
        cout << "Color del nodo 1: " << colores[1] << " | Color del nodo 2: " << colores[2] << "\n";
    }
    */

    return 0;
}