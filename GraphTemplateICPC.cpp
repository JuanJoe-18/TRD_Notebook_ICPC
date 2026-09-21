/**
 * @file GraphTemplateICPC.cpp
 * @brief Plantilla de algoritmos y estructuras para grafos.
 * @details Incluye DSU, caminos minimos, SCC, ciclos, Euler, LCA y 2-SAT.
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

    /**
     * @brief Inicializa el DSU con n elementos.
     * @param n Cantidad de elementos (0-indexed o 1-indexed, usar tamaño n+1 si es 1-indexed).
     */
    DSU(int n) : p(n), sz(n, 1), members(n), num_components(n) {
        iota(p.begin(), p.end(), 0);
        for (int i = 0; i < n; i++) {
            members[i].push_back(i);
        }
    }

    /**
     * @brief Encuentra el representante del conjunto al que pertenece x.
     * @param x Elemento a consultar.
     * @return Representante del conjunto.
     */
    int find(int x) {
        return p[x] == x ? x : p[x] = find(p[x]);
    }

    /**
     * @brief Une los conjuntos a los que pertenecen a y b.
     * @param a Primer elemento.
     * @param b Segundo elemento.
     * @return true si estaban en distintos conjuntos y se unieron, false si ya estaban unidos.
     */
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

    /**
     * @brief Devuelve los nodos que pertenecen a la componente de x.
     * @param x Nodo cuya componente se quiere consultar.
     * @return Referencia constante al vector de nodos de la componente.
     */
    const vector<int>& get_component_nodes(int x) {
        return members[find(x)];
    }

    /**
     * @brief Devuelve todas las componentes conexas actuales.
     * @param first_node Primer indice que se debe considerar. Usa 1 con DSU(n+1)
     *        cuando el nodo 0 se deja como dummy.
     * @return Vector con los nodos agrupados por componente.
     */
    vector<vector<int>> get_components(int first_node = 0) const {
        vector<vector<int>> result;
        for (int root = first_node; root < (int)p.size(); root++) {
            if (p[root] == root && !members[root].empty()) {
                result.push_back(members[root]);
            }
        }
        return result;
    }

    /**
     * @brief Devuelve el tamaño de la componente de x.
     */
    int component_size(int x) {
        return sz[find(x)];
    }

    /**
     * @brief Devuelve la cantidad de componentes.
     * @param first_node Usa 1 para ignorar el nodo 0 dummy en DSU(n+1).
     */
    int component_count(int first_node = 0) const {
        if (first_node == 0) return num_components;
        int count = 0;
        for (int root = first_node; root < (int)p.size(); root++) {
            if (p[root] == root && !members[root].empty()) count++;
        }
        return count;
    }
};

// Ejemplo:
// DSU dsu(n + 1); // vertices 1..n; el indice 0 queda sin usar
// dsu.unite(1, 4);
// for (int node : dsu.get_component_nodes(4)) cout << node << ' ';
// for (const auto& component : dsu.get_components(1)) { /* usar component */ }

template <typename T = ll>
struct Graph {
    int n;
    vector<vector<int>> adj; // Guarda indices de las aristas
    vector<Edge<T>> edges;   // Guarda la informacion real de las aristas

    Graph(int _n) : n(_n) {
        adj.resize(n + 1);
    }

    // Agregar arista dirigida
    void add_directed_edge(int from, int to, T weight = 1, int id = -1) {
        adj[from].push_back(edges.size());
        edges.push_back({from, to, weight, id});
    }

    // Agregar arista bidireccional
    void add_undirected_edge(int u, int v, T weight = 1, int id = -1) {
        add_directed_edge(u, v, weight, id);
        add_directed_edge(v, u, weight, id);
    }

    // ---------------------------------------------------------
    // ALGORITMOS INTEGRADOS
    // ---------------------------------------------------------

    // 1. Dijkstra O(E log V)
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

    // 8. Floyd-Warshall O(V^3) - All-Pairs Shortest Path
    // Retorna una matriz de distancias. Útil solo si N <= 400.
    vector<vector<T>> floyd_warshall() {
        vector<vector<T>> dist(n + 1, vector<T>(n + 1, LINF));
        for (int i = 1; i <= n; i++) dist[i][i] = 0;
        for (auto& e : edges) {
            dist[e.from][e.to] = min(dist[e.from][e.to], e.weight);
            // Si el grafo es no dirigido, asegúrate de que 'edges' tenga ambas direcciones
        }

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

    // 2. Ordenamiento Topológico (Kahn's Algorithm)
    vector<int> topo_sort() {
        vector<int> in_degree(n + 1, 0), order;
        for (auto& e : edges) in_degree[e.to]++;

        queue<int> q;
        for (int i = 1; i <= n; i++)
            if (in_degree[i] == 0) q.push(i);

        while (!q.empty()) {
            int u = q.front(); q.pop();
            order.push_back(u);
            for (int id : adj[u]) {
                if (--in_degree[edges[id].to] == 0)
                    q.push(edges[id].to);
            }
        }
        return (order.size() == n) ? order : vector<int>(); // Vacio si hay ciclo
    }

    // 3. Tarjan para Componentes Fuertemente Conexas (SCC)
    // Devuelve un arreglo donde comp[i] es el ID de la componente del nodo i
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

    // 4. Bellman-Ford O(V E)
    // Retorna: <bool (ok), vector<T> (distancias), vector<int> (padres), vector<int> (ciclo negativo)>
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

    // 5. Encontrar Ciclo Dirigido (Round Trip II) - Blindado contra RE
    // Retorna: vector vacio si no hay ciclo. Si hay, retorna los nodos del ciclo en orden [c1, c2, ..., ck, c1].
    vector<int> find_directed_cycle() {
        vector<int> state(n + 1, 0); // 0: Blanco (no visitado), 1: Gris (en pila), 2: Negro (terminado)
        vector<int> parent(n + 1, -1);
        vector<int> cycle;

        auto dfs = [&](auto& self, int u) -> bool {
            state[u] = 1; // Gris
            for (int id : adj[u]) {
                int v = edges[id].to;
                if (state[v] == 0) {
                    parent[v] = u;
                    if (self(self, v)) return true;
                } else if (state[v] == 1) { // Ciclo dirigido encontrado
                    cycle.push_back(v);
                    int curr = u;
                    while (curr != -1 && curr != v) { // Blindado contra parent[-1]
                        cycle.push_back(curr);
                        curr = parent[curr];
                    }
                    cycle.push_back(v);
                    reverse(cycle.begin(), cycle.end());
                    return true;
                }
            }
            state[u] = 2; // Negro
            return false;
        };

        for (int i = 1; i <= n; i++) {
            if (state[i] == 0) {
                if (dfs(dfs, i)) return cycle;
            }
        }
        return cycle;
    }

    // 6. Encontrar Ciclo No Dirigido (Round Trip I) - Blindado contra multi-aristas y RE
    // Retorna: vector vacio si no hay ciclo. Si hay, retorna los nodos del ciclo en orden [c1, c2, ..., ck, c1].
    vector<int> find_undirected_cycle() {
        vector<int> state(n + 1, 0); // 0: Blanco, 1: Gris (en pila), 2: Negro
        vector<int> parent_node(n + 1, -1);
        vector<int> cycle;

        auto dfs = [&](auto& self, int u, int p_edge_idx) -> bool {
            state[u] = 1;
            for (int id : adj[u]) {
                // Si la arista es la inversa de la que usamos para llegar a u, la ignoramos
                if ((id ^ 1) == p_edge_idx) continue;
                int v = edges[id].to;

                if (state[v] == 0) {
                    parent_node[v] = u;
                    if (self(self, v, id)) return true;
                } else if (state[v] == 1) { // Ciclo encontrado (back-edge hacia un ancestro)
                    cycle.push_back(v);
                    int curr = u;
                    while (curr != -1 && curr != v) { // Blindado contra desbordamiento
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

        for (int i = 1; i <= n; i++) {
            if (state[i] == 0) {
                if (dfs(dfs, i, -1)) return cycle;
            }
        }
        return cycle;
    }

    // 7. Algoritmo de Kruskal para Minimum Spanning Tree (MST)
    /**
     * @brief Calcula el Árbol de Expansión Mínima (MST).
     * @return Un par donde el primer elemento es el peso total del MST,
     *         y el segundo elemento es un vector con los índices de las aristas del MST.
     */
    pair<T, vector<int>> kruskal() {
        vector<int> edge_indices(edges.size());
        iota(edge_indices.begin(), edge_indices.end(), 0);

        // Ordenar índices de las aristas por peso
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

    // 8. Circuito Euleriano (Hierholzer's Algorithm) O(V + E)
    // IMPORTANTE: Para grafos no dirigidos, debes pasar un 'id' unico (0, 1, 2...) en add_undirected_edge.
    // Retorna vector vacio si no existe circuito.
    vector<int> eulerian_circuit(int start_node = 1, bool undirected = true) {
        vector<int> in_deg(n + 1, 0), out_deg(n + 1, 0);
        int max_id = -1;
        for (auto& e : edges) {
            out_deg[e.from]++;
            in_deg[e.to]++;
            max_id = max(max_id, e.id);
        }

        // 1. Verificación de grados (Condición de Euler)
        for (int i = 1; i <= n; i++) {
            if (undirected) {
                if (out_deg[i] % 2 != 0) return {};
            } else {
                if (in_deg[i] != out_deg[i]) return {};
            }
        }

        vector<bool> used_edge(max_id + 1, false);
        vector<int> circuit;
        vector<int> head(n + 1, 0); // Puntero O(E) para iterar adj

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
                    // Si es dirigido y no pasaron ID, el puntero head es suficiente
                    self(self, e.to);
                }
            }
            circuit.push_back(u);
        };

        dfs(dfs, start_node);
        reverse(circuit.begin(), circuit.end());

        // 2. Verificación de conectividad (deben haberse usado todas las aristas)
        // Para grafos no dirigidos, circuit.size() debe ser igual a M + 1 (donde M = edges.size() / 2)
        int required_size = undirected ? (edges.size() / 2 + 1) : (edges.size() + 1);
        if (circuit.size() != required_size && edges.size() > 0) return {};

        return circuit;
    }
};

// GRAFOS FUNCIONALES (Successor Graphs)
// ==========================================
// Cada nodo u tiene un único sucesor succ[u] (out-degree = 1).
// Incluye:
// - Descomposición de TODOS los ciclos en O(N) sin recursión profunda.
// - Binary Lifting para saltar K pasos en O(log K).
// - Identificación de nodos en ciclo, distancias a ciclo y tamaños de ciclo.
struct FunctionalGraph {
    int n, log_k;
    vector<vector<int>> up;
    vector<int> succ;

    // Información de ciclos
    vector<vector<int>> cycles;    // Lista con todos los ciclos [ [c1, c2, c3], ... ]
    vector<int> in_cycle;          // 1 si el nodo u pertenece a un ciclo, 0 si no
    vector<int> cycle_id;          // ID del ciclo al que pertenece o al que llega
    vector<int> cycle_pos;         // Posición (0-indexed) del nodo dentro de su ciclo
    vector<int> cycle_size;        // Tamaño del ciclo al que pertenece o al que llega
    vector<int> dist_to_cycle;     // Distancia en pasos desde u hasta entrar a su ciclo

    FunctionalGraph(int _n, int max_k_bits = 60) : n(_n), log_k(max_k_bits) {
        up.assign(n + 1, vector<int>(log_k, 0));
        succ.assign(n + 1, 0);
        in_cycle.assign(n + 1, 0);
        cycle_id.assign(n + 1, -1);
        cycle_pos.assign(n + 1, -1);
        cycle_size.assign(n + 1, 0);
        dist_to_cycle.assign(n + 1, 0);
    }

    // 1. Descomposición completa de ciclos en O(N) lineal e iterativo (cero riesgo de Stack Overflow)
    void decompose_cycles(const vector<int>& _succ) {
        succ = _succ;
        vector<int> state(n + 1, 0); // 0: No visitado, 1: En camino actual, 2: Terminado
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

            // Si chocamos con un nodo en la ruta actual, encontramos un ciclo nuevo
            if (curr >= 1 && curr <= n && state[curr] == 1) {
                vector<int> cyc;
                bool found_start = false;
                for (int u : path) {
                    if (u == curr) found_start = true;
                    if (found_start) cyc.push_back(u);
                }

                int cid = cycles.size();
                int c_sz = cyc.size();
                for (int pos = 0; pos < c_sz; pos++) {
                    int u = cyc[pos];
                    in_cycle[u] = 1;
                    cycle_id[u] = cid;
                    cycle_pos[u] = pos;
                    cycle_size[u] = c_sz;
                    dist_to_cycle[u] = 0;
                }
                cycles.push_back(cyc);
            }

            for (int u : path) state[u] = 2;
        }

        // Calcular distancias y ciclo destino para los nodos que no están en ciclo
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
                d++;
                int u = path[j];
                cycle_id[u] = target_cid;
                cycle_size[u] = target_csz;
                dist_to_cycle[u] = d;
            }
        }
    }

    // 2. Binary Lifting para consultas de saltos en O(log K)
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

    // Salto de K pasos en O(log K)
    int get_kth_successor(int u, ll k) {
        for (int j = 0; j < log_k; j++) {
            if (k & (1LL << j)) {
                u = up[u][j];
                if (u < 1 || u > n) return 0;
            }
        }
        return u;
    }

    // Verifica si existe al menos un ciclo de longitud exacta K en O(Número de Ciclos) <= O(N)
    bool has_cycle_of_length(int target_len) {
        for (auto& cyc : cycles) {
            if ((int)cyc.size() == target_len) return true;
        }
        return false;
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

    // adj = lista de adyacencia bidireccional del árbol
    void dfs(int u, int p, const vector<vector<int>>& adj) {
        up[u][0] = p;
        for (int j = 1; j < log_n; j++) {
            up[u][j] = up[up[u][j - 1]][j - 1];
        }
        for (int v : adj[u]) {
            if (v != p) {
                depth[v] = depth[u] + 1;
                dfs(v, u, adj);
            }
        }
    }

    void build(int root, const vector<vector<int>>& adj) {
        dfs(root, root, adj);
    }

    int get_lca(int u, int v) {
        if (depth[u] < depth[v]) swap(u, v);
        int k = depth[u] - depth[v];
        for (int j = 0; j < log_n; j++) {
            if (k & (1 << j)) u = up[u][j];
        }
        if (u == v) return u;
        for (int j = log_n - 1; j >= 0; j--) {
            if (up[u][j] != up[v][j]) {
                u = up[u][j];
                v = up[v][j];
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

    // n es el número de variables booleanas
    TwoSat(int _n) : n(_n), G(2 * _n) {
        assignment.assign(n + 1, false);
    }

    // Retorna el nodo de la variable u. Si is_true es false, retorna su negación.
    int get_node(int u, bool is_true) {
        return is_true ? u : u + n;
    }

    // Agrega la cláusula (u OR v)
    // Ejemplo: Si quiero "X_2 o NO X_3", llamo add_clause(2, true, 3, false)
    void add_clause(int u, bool is_u_true, int v, bool is_v_true) {
        int not_u = get_node(u, !is_u_true);
        int node_v = get_node(v, is_v_true);
        int not_v = get_node(v, !is_v_true);
        int node_u = get_node(u, is_u_true);

        G.add_directed_edge(not_u, node_v); // !u -> v
        G.add_directed_edge(not_v, node_u); // !v -> u
    }

    // Forzar que una variable sea obligatoriamente Verdadera o Falsa
    void force_value(int u, bool is_true) {
        add_clause(u, is_true, u, is_true);
    }

    // Agregar implicación: Si U pasa, entonces V tiene que pasar obligatoriamente (U => V)
    // Es lógicamente equivalente a (!U OR V)
    void add_implication(int u, bool is_u_true, int v, bool is_v_true) {
        add_clause(u, !is_u_true, v, is_v_true);
    }

    // Agregar XOR: u y v deben tener valores DISTINTOS (u != v)
    // Equivale a (u OR v) AND (!u OR !v)
    void add_xor(int u, bool is_u_true, int v, bool is_v_true) {
        add_clause(u, is_u_true, v, is_v_true);
        add_clause(u, !is_u_true, v, !is_v_true);
    }

    // Agregar XNOR / Equivalencia: u y v deben tener el MISMO valor (u == v)
    // Equivale a (!u OR v) AND (u OR !v)
    void add_equivalence(int u, bool is_u_true, int v, bool is_v_true) {
        add_implication(u, is_u_true, v, is_v_true);
        add_implication(v, is_v_true, u, is_u_true);
    }

    // Intenta resolver el 2-SAT. Retorna true si es posible.
    // Los resultados quedan en el arreglo booleano 'assignment'
    bool solve() {
        vector<int> comp = G.get_scc(); // Tarjan O(V+E)

        for (int i = 1; i <= n; i++) {
            if (comp[i] == comp[i + n]) {
                return false; // Contradicción: u y !u están en el mismo ciclo
            }
            // Magia de Tarjan: Los componentes terminados primero (menor ID) son sumideros.
            // Siempre asignamos True a los sumideros para no forzar errores hacia atrás.
            assignment[i] = comp[i] < comp[i + n];
        }
        return true;
    }
};

// ==========================================
// EMPAREJAMIENTO BIPARTITO (Kuhn's Algorithm) O(V * E)
// ==========================================
// Útil para Maximum Bipartite Matching. Para N, M <= 1000.
struct BipartiteMatcher {
    int n, m;
    vector<vector<int>> adj;
    vector<int> match; // match[v] = nodo de 'U' emparejado con 'v' (de 'V')
    vector<bool> vis;

    // n = tamaño del conjunto U (izquierdo), m = tamaño del conjunto V (derecho)
    BipartiteMatcher(int _n, int _m) : n(_n), m(_m) {
        adj.resize(n + 1);
        match.assign(m + 1, -1);
    }

    void add_edge(int u, int v) {
        adj[u].push_back(v);
    }

    bool dfs(int u) {
        for (int v : adj[u]) {
            if (vis[v]) continue;
            vis[v] = true;
            if (match[v] < 0 || dfs(match[v])) {
                match[v] = u;
                return true;
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
// PUENTES, PUNTOS DE ARTICULACIÓN Y BRIDGE TREE
// ==========================================
struct BridgeTree {
    int n, timer, num_comps;
    vector<vector<pair<int, int>>> adj; // {v, edge_id}
    vector<int> tin, low, comp_id;
    vector<bool> is_bridge, is_articulation;

    // El Bridge Tree (Árbol condensado de componentes biconexas)
    vector<vector<int>> tree_adj;

    BridgeTree(int _n, int num_edges) : n(_n) {
        adj.resize(n + 1);
        tin.assign(n + 1, -1);
        low.assign(n + 1, -1);
        comp_id.assign(n + 1, 0);
        is_bridge.assign(num_edges, false);
        is_articulation.assign(n + 1, false);
        timer = 0;
        num_comps = 0;
    }

    void add_edge(int u, int v, int id) {
        adj[u].push_back({v, id});
        adj[v].push_back({u, id});
    }

    void dfs_tarjan(int u, int p = -1) {
        tin[u] = low[u] = ++timer;
        int children = 0;
        for (auto& edge : adj[u]) {
            int v = edge.first;
            int id = edge.second;
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
            int v = edge.first;
            int id = edge.second;
            if (comp_id[v] == 0) {
                if (is_bridge[id]) {
                    num_comps++;
                    tree_adj.push_back(vector<int>()); // Nuevo nodo en el árbol
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
        for (int i = 1; i <= n; i++) {
            if (tin[i] == -1) dfs_tarjan(i);
        }
        tree_adj.push_back(vector<int>()); // Dummy 0
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
// Muy rápido en la práctica. O(E * sqrt(V)) en grafos bipartitos.
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
        adj.resize(n + 1);
        level.resize(n + 1);
        ptr.resize(n + 1);
    }

    void add_edge(int u, int v, ll cap, bool directed = true) {
        adj[u].push_back(edges.size());
        edges.push_back(FlowEdge(u, v, cap));
        adj[v].push_back(edges.size());
        edges.push_back(FlowEdge(v, u, directed ? 0 : cap));
    }

    bool bfs() {
        fill(level.begin(), level.end(), -1);
        level[s] = 0;
        queue<int> q;
        q.push(s);
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
            while (ll pushed = dfs(s, LINF)) {
                flow += pushed;
            }
        }
        return flow;
    }
};

// ==========================================
// MIN COST MAX FLOW (SPFA-based Successive Shortest Path)
// ==========================================
struct MCMF {
    struct Edge {
        int to;
        ll cap, flow, cost;
        int rev;
    };

    int n;
    vector<vector<Edge>> adj;
    vector<ll> dist;
    vector<int> p_node, p_edge;
    vector<bool> in_queue;

    MCMF(int _n) : n(_n) {
        adj.resize(n + 1);
        dist.resize(n + 1);
        p_node.resize(n + 1);
        p_edge.resize(n + 1);
        in_queue.resize(n + 1);
    }

    void add_edge(int u, int v, ll cap, ll cost) {
        adj[u].push_back({v, cap, 0, cost, (int)adj[v].size()});
        adj[v].push_back({u, 0, 0, -cost, (int)adj[u].size() - 1});
    }

    bool spfa(int s, int t) {
        fill(dist.begin(), dist.end(), LINF);
        fill(in_queue.begin(), in_queue.end(), false);
        queue<int> q;

        dist[s] = 0;
        q.push(s);
        in_queue[s] = true;

        while (!q.empty()) {
            int u = q.front(); q.pop();
            in_queue[u] = false;

            for (int i = 0; i < adj[u].size(); i++) {
                auto& e = adj[u][i];
                if (e.cap - e.flow > 0 && dist[e.to] > dist[u] + e.cost) {
                    dist[e.to] = dist[u] + e.cost;
                    p_node[e.to] = u;
                    p_edge[e.to] = i;
                    if (!in_queue[e.to]) {
                        q.push(e.to);
                        in_queue[e.to] = true;
                    }
                }
            }
        }
        return dist[t] != LINF;
    }

    // Retorna {Max Flow, Min Cost}
    pair<ll, ll> solve(int s, int t) {
        ll flow = 0, cost = 0;
        while (spfa(s, t)) {
            ll push = LINF;
            int curr = t;
            while (curr != s) {
                int p = p_node[curr];
                int idx = p_edge[curr];
                push = min(push, adj[p][idx].cap - adj[p][idx].flow);
                curr = p;
            }

            flow += push;
            curr = t;
            while (curr != s) {
                int p = p_node[curr];
                int idx = p_edge[curr];
                int rev_idx = adj[p][idx].rev;

                adj[p][idx].flow += push;
                adj[curr][rev_idx].flow -= push;
                cost += push * adj[p][idx].cost;
                curr = p;
            }
        }
        return {flow, cost};
    }
};

int main() {
    // ==========================================
    // EJEMPLOS DE USO MÁGICOS
    // ==========================================

    // ---------------------------------------------------------
    // 1. DIJKSTRA & FLOYD-WARSHALL (Caminos Mínimos)
    // ---------------------------------------------------------
    // cout << "--- 1. SHORTEST PATHS ---\n";
    // Instanciamos un grafo g1 de 3 nodos
    // Graph<ll> g1(3);
    // Agregamos aristas dirigidas (origen, destino, peso)
    // g1.add_directed_edge(1, 2, 5);
    // g1.add_directed_edge(2, 3, 10);
    // g1.add_directed_edge(1, 3, 20);

    // Calculamos los caminos mínimos desde el nodo 1 usando Dijkstra
    // vector<ll> dist1 = g1.dijkstra(1);
    // Imprimimos la distancia mínima al nodo 3 (la ruta óptima es 1->2->3 = 15)
    // cout << "Dijkstra (1 -> 3): " << dist1[3] << "\n";

    // Calculamos todas las distancias entre todos los pares (solo si V <= 400)
    // auto fw = g1.floyd_warshall();
    // Verificamos que la distancia en la matriz de 1 a 3 coincida (15)
    // cout << "Floyd-Warshall (1 -> 3): " << fw[1][3] << "\n\n";

    // ---------------------------------------------------------
    // 2. SCC (Tarjan) & ORDENAMIENTO TOPOLÓGICO
    // ---------------------------------------------------------
    // cout << "--- 2. SCC & TOPO SORT ---\n";
    // Grafo g2 de 4 nodos para buscar componentes fuertemente conexas
    // Graph<ll> g2(4);
    // Creamos un ciclo entre 1, 2 y 3
    // g2.add_directed_edge(1, 2);
    // g2.add_directed_edge(2, 3);
    // g2.add_directed_edge(3, 1);
    // Arista de escape hacia el nodo 4 (que no pertenece al ciclo)
    // g2.add_directed_edge(3, 4);

    // Tarjan O(V+E) asigna un ID único a cada componente biconexa
    // vector<int> scc = g2.get_scc();
    // 1, 2 y 3 comparten ID. El 4 tendrá un ID menor por ser sumidero topológico.
    // cout << "IDs SCC -> Nodo 1: " << scc[1] << ", Nodo 3: " << scc[3] << ", Nodo 4: " << scc[4] << "\n";

    // Grafo Acíclico Dirigido (DAG) de 3 nodos para Topo Sort
    // Graph<ll> dag(3);
    // dag.add_directed_edge(1, 2);
    // dag.add_directed_edge(1, 3);
    // dag.add_directed_edge(2, 3);
    // Imprimimos el orden topológico (un nodo solo aparece si sus dependencias ya pasaron)
    // cout << "Topo Sort de un DAG: ";
    // for(int x : dag.topo_sort()) cout << x << " "; // Output: 1 2 3
    // cout << "\n\n";

    // ---------------------------------------------------------
    // 3. KRUSKAL (MST) & DSU
    // ---------------------------------------------------------
    // cout << "--- 3. KRUSKAL (MST) ---\n";
    // Grafo no dirigido g3 de 3 nodos para el Árbol de Expansión Mínima
    // Graph<ll> g3(3);
    // Aristas bidireccionales con pesos
    // g3.add_undirected_edge(1, 2, 10);
    // g3.add_undirected_edge(2, 3, 20);
    // g3.add_undirected_edge(1, 3, 5);

    // Extraemos el peso total y los índices de las aristas usadas en el MST
    // auto [mst_weight, mst_edges] = g3.kruskal();
    // Debe dar 15 (usa las aristas de peso 5 y 10 para conectar los 3 nodos)
    // cout << "Peso del Arbol de Expansion Minima: " << mst_weight << "\n\n";

    // ---------------------------------------------------------
    // 4. 2-SAT (Satisfactibilidad)
    // ---------------------------------------------------------
    // cout << "--- 4. 2-SAT ---\n";
    // Instanciamos 2-SAT para 3 variables booleanas (x1, x2, x3)
    // TwoSat ts(3);
    // Cláusula 1: (x1 OR x2)
    // ts.add_clause(1, true, 2, true);
    // Cláusula 2: (!x1 OR x3)
    // ts.add_clause(1, false, 3, true);
    // Cláusula 3: (!x2 OR !x3)
    // ts.add_clause(2, false, 3, false);

    // resolvemos la ecuación booleana
    // if (ts.solve()) {
        // Imprimimos la asignación válida que encontró Tarjan
       //  cout << "Formula Satisfactible! x1=" << ts.assignment[1]
             << " x2=" << ts.assignment[2] << " x3=" << ts.assignment[3] << "\n\n";
    // } else {
        // cout << "No hay solucion.\n\n"; // Saldría si hubiera una contradicción insalvable
    // }

    // ---------------------------------------------------------
    // 5. FUNCTIONAL GRAPH (Grafos Sucesores)
    // ---------------------------------------------------------
    // cout << "--- 5. FUNCTIONAL GRAPH ---\n";
    // Grafo donde TODO nodo tiene out-degree == 1 exacto
    // FunctionalGraph fg(5);
    // Arreglo de sucesores. El índice 0 es ignorado (1-indexed)
    // 1->2, 2->3, 3->1 (Ciclo 1) | 4->5, 5->4 (Ciclo 2)
    // vector<int> succ = {0, 2, 3, 1, 5, 4};
    // Identifica todos los ciclos en O(N) de forma iterativa (anti Stack-Overflow)
    // fg.decompose_cycles(succ);

    // Imprime la cantidad de ciclos detectados (2 ciclos)
    // cout << "Cantidad de ciclos puros: " << fg.cycles.size() << "\n";
    // Confirma si el nodo 1 está atrapado en un ciclo (true)
    // cout << "Nodo 1 pertenece a un ciclo? " << fg.in_cycle[1] << "\n\n";

    // ---------------------------------------------------------
    // 6. EULERIAN CIRCUIT
    // ---------------------------------------------------------
    // cout << "--- 6. EULERIAN CIRCUIT ---\n";
    // Grafo g6 para buscar un camino que pase por todas las aristas exactamente una vez
    // Graph<ll> g6(3);
    // OBLIGATORIO: En grafos NO dirigidos, debes pasar IDs únicos (0, 1, 2)
    // para evitar cruzar la misma arista de ida y vuelta engañando al algoritmo
    // g6.add_undirected_edge(1, 2, 1, 0);
    // g6.add_undirected_edge(2, 3, 1, 1);
    // g6.add_undirected_edge(3, 1, 1, 2);

    // Genera la secuencia de nodos del circuito empezando en 1
    // cout << "Circuito Euleriano (Triangulo): ";
    // for (int u : g6.eulerian_circuit(1, true)) cout << u << " ";
    // cout << "\n\n"; // Imprime 1 2 3 1

    // ---------------------------------------------------------
    // 7. BIPARTITE MATCHING (Kuhn)
    // ---------------------------------------------------------
    // cout << "--- 7. BIPARTITE MATCHING ---\n";
    // 3 nodos en conjunto A (Trabajadores) y 3 en conjunto B (Tareas)
    // BipartiteMatcher bm(3, 3);
    // Añadimos las aristas (Trabajador -> Tarea que sabe hacer)
    // bm.add_edge(1, 1);
    // bm.add_edge(1, 2); // El T1 es polivalente
    // bm.add_edge(2, 2);
    // bm.add_edge(3, 3);

    // Kuhn calcula el "Maximum Bipartite Matching" en O(VE)
    // cout << "Maximas asignaciones posibles: " << bm.solve() << "\n\n"; // Da 3 tareas asignadas

    // ---------------------------------------------------------
    // 8. BRIDGE TREE (Condensación Biconexa)
    // ---------------------------------------------------------
    // cout << "--- 8. BRIDGE TREE ---\n";
    // 5 nodos, 5 aristas
    // BridgeTree btree(5, 5);
    // Parámetros: (u, v, edge_id_unico). Creamos un ciclo entre 1, 2 y 3.
    // btree.add_edge(1, 2, 0);
    // btree.add_edge(2, 3, 1);
    // btree.add_edge(3, 1, 2);
    // Aristas de escape que actuarán como PUENTES
    // btree.add_edge(3, 4, 3);
    // btree.add_edge(4, 5, 4);
    // Condensamos el grafo. Todo ciclo colapsa en un solo super-nodo.
    // btree.build();

    // El ciclo {1,2,3} es un supernodo, 4 es otro, 5 es otro. Total = 3 nodos en el árbol nuevo.
    // cout << "Cantidad de nodos en el arbol condensado: " << btree.num_comps << "\n\n";

    // ---------------------------------------------------------
    // 9. DINIC'S ALGORITHM (Max Flow)
    // ---------------------------------------------------------
    // cout << "--- 9. MAX FLOW (DINIC) ---\n";
    // Dinic(nodos, Source, Sink)
    // Dinic dinic(4, 1, 4);
    // Aristas: (origen, destino, capacidad)
    // dinic.add_edge(1, 2, 100);
    // dinic.add_edge(1, 3, 50);
    // dinic.add_edge(2, 4, 80);
    // dinic.add_edge(3, 4, 70);
    // Ejecuta BFS por niveles y DFS para encontrar cuellos de botella
    // cout << "Flujo maximo posible: " << dinic.max_flow() << "\n\n"; // Envia 80 por arriba y 50 por abajo = 130

    // ---------------------------------------------------------
    // 10. MIN COST MAX FLOW (MCMF)
    // ---------------------------------------------------------
    // cout << "--- 10. MIN COST MAX FLOW ---\n";
    // Red de flujo de 4 nodos
    // MCMF mcmf(4);
    // Aristas: (origen, destino, capacidad, costo_por_unidad)
    // mcmf.add_edge(1, 2, 10, 5);
    // mcmf.add_edge(1, 3, 10, 1); // Ruta más barata
    // mcmf.add_edge(2, 4, 10, 2);
    // mcmf.add_edge(3, 4, 10, 8);

    // Calcula la forma más barata de enviar el mayor flujo posible de 1 a 4
    // auto [flujo, costo] = mcmf.solve(1, 4);
    // Flujo=20. Costo=(10 cap * (5+2 costo)) + (10 cap * (1+8 costo)) = 70 + 90 = 160.
    // cout << "Max Flujo: " << flujo << ", Al Menor Costo: " << costo << "\n";

    return 0;
}
