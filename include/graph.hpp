/**
 * Uso[1]: dx4, dy4, dx8, dy8, isValid(x, y, R, C);
 * Direcciones de movimiento en grillas 4/8-conexas y validacion de coordenadas dentro de la grilla.
 * Complejidad: $O(1)$.
 */
const int dx4[4] = {1, -1, 0, 0};
const int dy4[4] = {0, 0, 1, -1};
const int dx8[8] = {1, 1, 0, -1, -1, -1, 0, 1};
const int dy8[8] = {0, 1, 1, 1, 0, -1, -1, -1};
#define isValid(x, y, R, C) ((x) >= 0 && (x) < (R) && (y) >= 0 && (y) < (C))

/**
 * Uso: Edge<T> e = {from, to, id, weight};
 * Arista con extremos, identificador opcional y peso generico; usada internamente por Graph.
 * Complejidad: $O(1)$ por arista.
 */
template <typename T = ll> struct Edge { int from, to, id; T weight; };

/**
 * Uso[0]: DSU dsu(n); dsu.find(x); dsu.unite(a, b); dsu.component_size(x); dsu.comps; dsu.get_components();
 * Union-find con union por tamano y acceso a los miembros de cada componente.
 * Complejidad: find y unite casi $O(\alpha(N))$ amortizado.
 */
struct DSU {
  vector<int> p, sz;
  vector<vector<int>> members;
  int comps;
DSU(int n) : p(n), sz(n, 1), members(n), comps(n) {
    iota(all(p), 0);
    for (int i = 0; i < n; i++) members[i].pb(i);
  }
  int find(int x) { return p[x] == x ? x : p[x] = find(p[x]); }
  bool unite(int a, int b) {
    a = find(a); b = find(b);
    if (a == b) return false;
    if (sz[a] < sz[b]) swap(a, b);
    p[b] = a; sz[a] += sz[b];
    members[a].insert(members[a].end(), all(members[b]));
    vector<int>().swap(members[b]);
    comps--;
    return true;
  }
  const vector<int>& component_nodes(int x) { return members[find(x)]; }
  int component_size(int x) { return sz[find(x)]; }
  vector<vector<int>> get_components(int first = 0) const {
    vector<vector<int>> res;
    for (int i = first; i < p.size(); i++) if (p[i] == i && !members[i].empty()) res.pb(members[i]);
    return res;
  }
};

/**
 * Uso[1]: Graph<ll> g(n); g.add_undirected_edge(u, v, w); g.add_directed_edge(u, v, w, id);
 * Grafo con lista de adyacencia de indices a aristas y pesos genericos.
 * Complejidad: agregar arista $O(1)$.
 */
template <typename T = ll> struct Graph {
  int n;
  vector<vector<int>> adj; // indices a edges
  vector<Edge<T>> edges;
  Graph(int n) : n(n) { adj.resize(n + 1); }
  void add_directed_edge(int u, int v, T w = 1, int id = -1) {
    adj[u].pb(edges.size());
    edges.pb({u, v, id, w});
  }
  void add_undirected_edge(int u, int v, T w = 1, int id = -1) {
    add_directed_edge(u, v, w, id);
    add_directed_edge(v, u, w, id);
  }
};

/**
 * Uso[1]: Dijkstra<ll> dih; 
 * auto [ dist, par ] = dih.run(g, src, invalid_edges);
 * auto path = dih.restore_path(target);
 * Caminos minimos desde un origen con pesos no negativos.
 * Complejidad: O(E log V).
 */
template <typename T = ll> struct Dijkstra {
  vector<T> dist;
  vi p;
  pair<vector<T>, vi> run(const Graph<T>& g, int src, const vector<bool>& invalid_edges = {}) {
    dist.assign(g.n + 1, LINF);
    p.assign(g.n + 1, -1);
    priority_queue<pair<T, int>, vector<pair<T, int>>, greater<>> pq;
    dist[src] = 0; pq.push({0, src});
    
    while (!pq.empty()) {
      auto [d, u] = pq.top(); pq.pop();
      if (d > dist[u]) continue;
      
      for (int idx : g.adj[u]) {
        auto& e = g.edges[idx];
        
        if (!invalid_edges.empty() && e.id != -1 && invalid_edges[e.id]) continue;
        
        if (dist[u] + e.weight < dist[e.to]) {
          dist[e.to] = dist[u] + e.weight;
          p[e.to] = u;
          pq.push({dist[e.to], e.to});
        }
      }
    }
    return {dist, p};
  }
  vi restore_path(int target) {
    vi path;
    if (dist.empty() || dist[target] == LINF) return path; 
    for (int v = target; v != -1; v = p[v]) path.pb(v);
    reverse(all(path));
    return path;
  }
};

/**
 * Uso[1]: auto dist = ZeroOneBfs<ll>().run(g, src);
 * Caminos minimos con pesos restringidos a $\{0, 1\}$ mediante un deque.
 * Complejidad: $O(V + E)$.
 */
template <typename T = ll> struct ZeroOneBfs {
  vector<T> run(const Graph<T>& g, int src) {
    vector<T> dist(g.n + 1, LINF);
    deque<int> dq;
    dist[src] = 0; dq.pb(src);
    while (!dq.empty()) {
      int u = dq.front(); dq.pop_front();
      for (int id : g.adj[u]) {
        auto& e = g.edges[id];
        if (dist[u] + e.weight < dist[e.to]) {
          dist[e.to] = dist[u] + e.weight;
          e.weight ? dq.pb(e.to) : dq.push_front(e.to);
        }
      }
    }
    return dist;
  }
};

/**
 * Uso[1]: auto [es_bipartito, colores] = Bipartite<ll>().check(g);
 * Verifica si el grafo es bipartito y colorea cada componente con 0/1.
 * Complejidad: $O(V + E)$.
 */
template <typename T = ll> struct Bipartite {
  pair<bool, vi> check(const Graph<T>& g) {
    vi color(g.n + 1, -1);
    bool ok = true;
    for (int i = 1; i <= g.n && ok; i++) if (color[i] == -1) {
      queue<int> q;
      q.push(i); color[i] = 0;
      while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int id : g.adj[u]) {
          int v = g.edges[id].to;
          if (color[v] == -1) color[v] = color[u] ^ 1, q.push(v);
          else if (color[v] == color[u]) ok = false;
        }
      }
    }
    return {ok, color};
  }
};

/**
 * Uso[1]: auto dist = FloydWarshall<ll>().run(g);
 * Caminos minimos entre todos los pares de nodos.
 * Complejidad: $O(V^3)$ tiempo, $O(V^2)$ espacio.
 */
template <typename T = ll> struct FloydWarshall {
  vector<vector<T>> run(const Graph<T>& g) {
    vector<vector<T>> dist(g.n + 1, vector<T>(g.n + 1, LINF));
    for (int i = 1; i <= g.n; i++) dist[i][i] = 0;
    for (auto& e : g.edges) dist[e.from][e.to] = min(dist[e.from][e.to], e.weight);
    for (int k = 1; k <= g.n; k++)
      for (int i = 1; i <= g.n; i++)
        for (int j = 1; j <= g.n; j++)
          if (dist[i][k] < LINF && dist[k][j] < LINF) dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
    return dist;
  }
};

/**
 * Uso[1]: auto order = TopoSort<ll>().run(g);
 * Orden topologico de un grafo dirigido aciclico; devuelve vacio si existe un ciclo.
 * Complejidad: $O(V + E)$.
 */
template <typename T = ll> struct TopoSort {
  vi run(const Graph<T>& g) {
    vi in(g.n + 1, 0), order;
    for (auto& e : g.edges) in[e.to]++;
    queue<int> q;
    for (int i = 1; i <= g.n; i++) if (!in[i]) q.push(i);
    while (!q.empty()) {
      int u = q.front(); q.pop();
      order.pb(u);
      for (int id : g.adj[u]) if (--in[g.edges[id].to] == 0) q.push(g.edges[id].to);
    }
    return order.size() == g.n ? order : vi();
  }
};

/**
 * Uso[1]: auto comp = TarjanSCC<ll>().run(g);
 * Componentes fuertemente conexas mediante Tarjan; devuelve el id de componente de cada nodo.
 * Complejidad: $O(V + E)$.
 */
template <typename T = ll> struct TarjanSCC {
  vi run(const Graph<T>& g) {
    int n = g.n, timer = 0, ncomps = 0;
    vi val(n + 1, 0), comp(n + 1, -1), stk;
    auto dfs = [&](auto&& self, int u) -> int {
      int low = val[u] = ++timer;
      stk.pb(u);
      for (int id : g.adj[u]) {
        int v = g.edges[id].to;
        if (comp[v] < 0) low = min(low, val[v] ? val[v] : self(self, v));
      }
      if (low == val[u]) {
        int x;
        do { x = stk.back(); stk.pop_back(); comp[x] = ncomps; } while (x != u);
        ncomps++;
      }
      return val[u] = low;
    };
    for (int i = 1; i <= n; i++) if (!val[i]) dfs(dfs, i);
    return comp;
  }
};

/**
 * Uso[1]: auto [sin_ciclo_neg, dist, parent, ciclo] = BellmanFord<ll>().run(g, src);
 * Caminos minimos con pesos negativos; detecta y recupera un ciclo negativo si existe.
 * Complejidad: $O(V \cdot E)$.
 */
template <typename T = ll> struct BellmanFord {
  tuple<bool, vector<T>, vi, vi> run(const Graph<T>& g, int src) {
    vector<T> dist(g.n + 1, LINF);
    vi p(g.n + 1, -1);
    dist[src] = 0;
    int x = -1;
    for (int i = 0; i < g.n; i++) {
      x = -1;
      for (auto& e : g.edges) if (dist[e.from] < LINF && dist[e.to] > dist[e.from] + e.weight) {
        dist[e.to] = max(-LINF, dist[e.from] + e.weight);
        p[e.to] = e.from; x = e.to;
      }
    }
    vi cycle;
    bool ok = true;
    if (x != -1) {
      ok = false;
      int y = x;
      for (int i = 0; i < g.n; i++) y = p[y];
      for (int cur = y;; cur = p[cur]) { cycle.pb(cur); if (cur == y && cycle.size() > 1) break; }
      reverse(all(cycle));
    }
    return {ok, dist, p, cycle};
  }
};

/**
 * Uso[1]: auto ciclo = DirectedCycle<ll>().find(g);
 * Encuentra un ciclo en un grafo dirigido (vacio si no existe).
 * Complejidad: $O(V + E)$.
 */
template <typename T = ll> struct DirectedCycle {
  vi find(const Graph<T>& g) {
    int n = g.n;
    vi state(n + 1, 0), parent(n + 1, -1), cycle;
    auto dfs = [&](auto&& self, int u) -> bool {
      state[u] = 1;
      for (int id : g.adj[u]) {
        int v = g.edges[id].to;
        if (state[v] == 0) {
          parent[v] = u;
          if (self(self, v)) return true;
        } else if (state[v] == 1) {
          cycle.pb(v);
          for (int cur = u; cur != -1 && cur != v; cur = parent[cur]) cycle.pb(cur);
          cycle.pb(v);
          reverse(all(cycle));
          return true;
        }
      }
      state[u] = 2;
      return false;
    };
    for (int i = 1; i <= n; i++) if (!state[i] && dfs(dfs, i)) return cycle;
    return cycle;
  }
};

/**
 * Uso[1]: auto ciclo = UndirectedCycle<ll>().find(g);
 * Encuentra un ciclo en un grafo no dirigido (vacio si no existe).
 * Complejidad: $O(V + E)$.
 */
template <typename T = ll> struct UndirectedCycle {
  vi find(const Graph<T>& g) {
    int n = g.n;
    vi state(n + 1, 0), parent(n + 1, -1), cycle;
    auto dfs = [&](auto&& self, int u, int pe) -> bool {
      state[u] = 1;
      for (int id : g.adj[u]) {
        if ((id ^ 1) == pe) continue;
        int v = g.edges[id].to;
        if (state[v] == 0) {
          parent[v] = u;
          if (self(self, v, id)) return true;
        } else if (state[v] == 1) {
          cycle.pb(v);
          for (int cur = u; cur != -1 && cur != v; cur = parent[cur]) cycle.pb(cur);
          cycle.pb(v);
          reverse(all(cycle));
          return true;
        }
      }
      state[u] = 2;
      return false;
    };
    for (int i = 1; i <= n; i++) if (!state[i] && dfs(dfs, i, -1)) return cycle;
    return cycle;
  }
};

/**
 * Uso[1]: auto [peso, aristas_mst] = Kruskal<ll>().run(g);
 * Arbol de expansion minima (MST) mediante Kruskal; devuelve el peso total y los indices de las aristas.
 * Complejidad: $O(E \log E)$.
 */
template <typename T = ll> struct Kruskal {
  pair<T, vi> run(const Graph<T>& g) {
    vi idx(g.edges.size());
    iota(all(idx), 0);
    sort(all(idx), [&](int a, int b) { return g.edges[a].weight < g.edges[b].weight; });
    DSU dsu(g.n + 1);
    T total = 0;
    vi mst;
    for (int i : idx) {
      auto& e = g.edges[i];
      if (dsu.unite(e.from, e.to)) {
        total += e.weight;
        mst.pb(i);
        if (mst.size() == g.n - 1) break;
      }
    }
    return {total, mst};
  }
};

/**
 * Uso[1]: auto circuito = EulerianCircuit<ll>().find(g, start, undirected);
 * Camino o circuito euleriano: recorre cada arista exactamente una vez (vacio si no existe).
 * Un circuito euleriano existe si el grafo es conexo y, en no dirigido, todos los grados son pares, o en dirigido, $in(u) = out(u)$ en cada nodo.
 * Complejidad: $O(V + E)$.
 */
template <typename T = ll> struct EulerianCircuit {
  vi find(const Graph<T>& g, int start = 1, bool undirected = true) {
    int n = g.n;
    vi in(n + 1, 0), out(n + 1, 0);
    int max_id = -1;
    for (auto& e : g.edges) {
      out[e.from]++; in[e.to]++;
      max_id = max(max_id, e.id);
    }
    for (int i = 1; i <= n; i++) {
      if (undirected ? (out[i] & 1) : (in[i] != out[i])) return {};
    }
    vector<bool> used(max_id + 1, false);
    vi circuit, head(n + 1, 0);
    auto dfs = [&](auto&& self, int u) -> void {
      while (head[u] < g.adj[u].size()) {
        auto& e = g.edges[g.adj[u][head[u]++]];
        if (e.id == -1 || !used[e.id]) {
          if (e.id != -1) used[e.id] = true;
          self(self, e.to);
        }
      }
      circuit.pb(u);
    };
    dfs(dfs, start);
    reverse(all(circuit));
    int need = undirected ? g.edges.size() / 2 + 1 : g.edges.size() + 1;
    if (circuit.size() != need && !g.edges.empty()) return {};
    return circuit;
  }
};

/**
 * Uso[1]: Hamiltonian h(n); h.add_edge(u, v); auto path = h.find_path(); auto cycle = h.find_cycle();
 * Camino hamiltoniano: visita cada vertice exactamente una vez; si existe, find_cycle devuelve un ciclo (el ultimo nodo conecta con el primero). Vacio si no hay.
 * Un camino hamiltoniano es NP-completo en general, por eso se resuelve por backtracking.
 * Indexación: nodos 1-indexados.
 * Complejidad: $O(N!)$ peor caso, practico para $N \leq 15$.
 */
struct Hamiltonian {
  int n;
  vector<vi> adj;
  vector<bool> vis;
  Hamiltonian(int n) : n(n), adj(n + 1), vis(n + 1, false) {}
  void add_edge(int u, int v) { adj[u].pb(v); adj[v].pb(u); }
  bool bt(int u, int cnt, int start, bool cycle, vi& path) {
    if (cnt == n) {
      if (!cycle) return true;
      for (int v : adj[u]) if (v == start) return true;
      return false;
    }
    for (int v : adj[u]) if (!vis[v]) {
      vis[v] = true; path.pb(v);
      if (bt(v, cnt + 1, start, cycle, path)) return true;
      path.pop_back(); vis[v] = false;
    }
    return false;
  }
  vi find(bool cycle) {
    for (int s = 1; s <= n; s++) {
      fill(all(vis), false);
      vi path = {s};
      vis[s] = true;
      if (bt(s, 1, s, cycle, path)) return path;
    }
    return {};
  }
  vi find_path() { return find(false); }  // camino que usa todos los nodos
  vi find_cycle() { return find(true); }  // ciclo que usa todos los nodos
};

/**
 * Uso[1]: FunctionalGraph fg(n); fg.build_lifting(succ); fg.kth_successor(u, k); fg.decompose(succ);
 * Grafo funcional (cada nodo un unico sucesor): k-esimo sucesor con binary lifting y descomposicion en ciclos.
 * Complejidad: preproceso $O(N \log K)$, k-esimo sucesor $O(\log K)$.
 */
struct FunctionalGraph {
  int n, log_k;
  vector<vector<int>> up;
  vi succ, in_cycle, cycle_id, cycle_pos, cycle_size, dist_to_cycle;
  vector<vi> cycles;
  FunctionalGraph(int n, int log_k = 60) : n(n), log_k(log_k),
    up(n + 1, vector<int>(log_k, 0)), succ(n + 1, 0), in_cycle(n + 1, 0),
    cycle_id(n + 1, -1), cycle_pos(n + 1, -1), cycle_size(n + 1, 0), dist_to_cycle(n + 1, 0) {}
  // Descompone el grafo en ciclos y arboles colgantes
  void decompose(const vi& _succ) {
    succ = _succ;
    vi state(n + 1, 0);
    cycles.clear();
    for (int i = 1; i <= n; i++) if (!state[i]) {
      int cur = i;
      vi path;
      while (cur >= 1 && cur <= n && state[cur] == 0) {
        state[cur] = 1;
        path.pb(cur);
        cur = succ[cur];
      }
      if (cur >= 1 && cur <= n && state[cur] == 1) {
        vi cyc;
        bool start = false;
        for (int u : path) { if (u == cur) start = true; if (start) cyc.pb(u); }
        int cid = cycles.size(), csz = cyc.size();
        for (int pos = 0; pos < csz; pos++) {
          int u = cyc[pos];
          in_cycle[u] = 1; cycle_id[u] = cid; cycle_pos[u] = pos; cycle_size[u] = csz;
          dist_to_cycle[u] = 0;
        }
        cycles.pb(cyc);
      }
      for (int u : path) state[u] = 2;
    }
    for (int i = 1; i <= n; i++) if (!in_cycle[i]) {
      int cur = i;
      vi path;
      while (cur >= 1 && cur <= n && !in_cycle[cur] && cycle_id[cur] == -1) path.pb(cur), cur = succ[cur];
      int cid = (cur >= 1 && cur <= n) ? cycle_id[cur] : -1;
      int csz = (cur >= 1 && cur <= n) ? cycle_size[cur] : 0;
      int d = (cur >= 1 && cur <= n) ? dist_to_cycle[cur] : 0;
      for (int j = path.size() - 1; j >= 0; j--) {
        int u = path[j];
        cycle_id[u] = cid; cycle_size[u] = csz; dist_to_cycle[u] = ++d;
      }
    }
  }
  // Binary lifting sobre succ para responder k-esimo sucesor
  void build_lifting(const vi& _succ) {
    succ = _succ;
    for (int i = 1; i <= n; i++) up[i][0] = succ[i];
    for (int j = 1; j < log_k; j++)
      for (int i = 1; i <= n; i++) {
        int p = up[i][j - 1];
        up[i][j] = (p >= 1 && p <= n) ? up[p][j - 1] : 0;
      }
  }
  int kth_successor(int u, ll k) {
    for (int j = 0; j < log_k; j++) if (k & (1LL << j)) {
      u = up[u][j];
      if (u < 1 || u > n) return 0;
    }
    return u;
  }
};

/**
 * Uso[1]: LCA<ll> lca(n); lca.build(root, g); lca.lca(u, v); lca.dist(u, v);
 * Ancestro comun mas bajo y distancia en numero de aristas mediante binary lifting, construido directamente sobre la estructura Graph.
 * Indexación: nodos 1-indexados.
 * Complejidad: preproceso $O(N \log N)$, consulta $O(\log N)$.
 */
template <typename T = ll> struct LCA {
  int n, log_n;
  vector<vi> up;
  vi depth;
  LCA(int n, int log_n = 20) : n(n), log_n(log_n), up(n + 1, vi(log_n, 0)), depth(n + 1, 0) {}
  void dfs(int u, int p, const Graph<T>& g) {
    up[u][0] = p;
    for (int j = 1; j < log_n; j++) up[u][j] = up[up[u][j - 1]][j - 1];
    for (int id : g.adj[u]) {
      int v = g.edges[id].to;
      if (v != p) depth[v] = depth[u] + 1, dfs(v, u, g);
    }
  }
  void build(int root, const Graph<T>& g) { dfs(root, root, g); }
  int lca(int u, int v) {
    if (depth[u] < depth[v]) swap(u, v);
    int k = depth[u] - depth[v];
    for (int j = 0; j < log_n; j++) if (k & (1 << j)) u = up[u][j];
    if (u == v) return u;
    for (int j = log_n - 1; j >= 0; j--) if (up[u][j] != up[v][j]) u = up[u][j], v = up[v][j];
    return up[u][0];
  }
  int dist(int u, int v) { return depth[u] + depth[v] - 2 * depth[lca(u, v)]; }
};

/**
 * Uso[1]: TwoSat s(n); s.add_clause(u, tu, v, tv); s.force_value(u, t); s.solve();
 * Satisfacibilidad booleana 2-SAT con clausulas (or, implicacion, xor, equivalencia) y la asignacion resultante.
 * Complejidad: $O(V + E)$.
 */
struct TwoSat {
  int n;
  Graph<ll> G;
  vector<bool> assignment;
  TwoSat(int n) : n(n), G(2 * n), assignment(n + 1, false) {}
  int node(int u, bool t) { return t ? u : u + n; }
  void add_clause(int u, bool tu, int v, bool tv) {
    G.add_directed_edge(node(u, !tu), node(v, tv));
    G.add_directed_edge(node(v, !tv), node(u, tu));
  }
  void force_value(int u, bool t) { add_clause(u, t, u, t); }
  void add_implication(int u, bool tu, int v, bool tv) { add_clause(u, !tu, v, tv); }
  void add_xor(int u, bool tu, int v, bool tv) {
    add_clause(u, tu, v, tv);
    add_clause(u, !tu, v, !tv);
  }
  void add_equivalence(int u, bool tu, int v, bool tv) {
    add_implication(u, tu, v, tv);
    add_implication(v, tv, u, tu);
  }
  bool solve() {
    vi comp = TarjanSCC<ll>().run(G);
    for (int i = 1; i <= n; i++) {
      if (comp[i] == comp[i + n]) return false;
      assignment[i] = comp[i] < comp[i + n];
    }
    return true;
  }
};

/**
 * Uso[1]: BipartiteMatcher bm(n_izq, n_der); bm.add_edge(u, v); bm.max_matching(); bm.match_of(v);
 * Emparejamiento maximo bipartito con el algoritmo de Kuhn; el resultado queda en match_of(v) = vertice izquierdo que empareja a v (o -1).
 * Indexación: nodos 1-indexados.
 * Complejidad: $O(V \cdot E)$.
 */
struct BipartiteMatcher {
  int n, m;
  vector<vi> adj;
  vi match;
  BipartiteMatcher(int n, int m) : n(n), m(m), adj(n + 1), match(m + 1, -1) {}
  void add_edge(int u, int v) { adj[u].pb(v); }
  bool dfs(int u, vector<bool>& vis) {
    for (int v : adj[u]) if (!vis[v]) {
      vis[v] = true;
      if (match[v] < 0 || dfs(match[v], vis)) { match[v] = u; return true; }
    }
    return false;
  }
  int max_matching() {
    int ans = 0;
    for (int i = 1; i <= n; i++) {
      vector<bool> vis(m + 1, false);
      if (dfs(i, vis)) ans++;
    }
    return ans;
  }
  int match_of(int v) { return match[v]; } // vertice izquierdo emparejado a v, -1 si v esta libre
};

/**
 * Uso[1]: HopcroftKarp hk(n_izq, n_der); hk.add_edge(u, v); hk.max_matching(); hk.match_of(v);
 * Emparejamiento maximo bipartito con Hopcroft-Karp; el resultado queda en match_of(v) = vertice izquierdo que empareja a v (o 0).
 * Indexación: nodos 1-indexados.
 * Complejidad: $O(E \sqrt{V})$.
 */
struct HopcroftKarp {
  int n, m;
  vector<vi> adj;
  vi pairU, pairV, dist;
  HopcroftKarp(int n, int m) : n(n), m(m), adj(n + 1), pairU(n + 1, 0), pairV(m + 1, 0), dist(n + 1) {}
  void add_edge(int u, int v) { adj[u].pb(v); }
  bool bfs() {
    queue<int> q;
    for (int u = 1; u <= n; u++) {
      if (pairU[u]) dist[u] = INF;
      else dist[u] = 0, q.push(u);
    }
    dist[0] = INF;
    while (!q.empty()) {
      int u = q.front(); q.pop();
      if (dist[u] < dist[0])
        for (int v : adj[u]) if (dist[pairV[v]] == INF) dist[pairV[v]] = dist[u] + 1, q.push(pairV[v]);
    }
    return dist[0] != INF;
  }
  bool dfs(int u) {
    if (!u) return true;
    for (int v : adj[u]) if (dist[pairV[v]] == dist[u] + 1 && dfs(pairV[v])) {
      pairV[v] = u; pairU[u] = v; return true;
    }
    dist[u] = INF;
    return false;
  }
  int max_matching() {
    int res = 0;
    while (bfs()) for (int u = 1; u <= n; u++) if (!pairU[u] && dfs(u)) res++;
    return res;
  }
  int match_of(int v) { return pairV[v]; } // vertice izquierdo emparejado a v, 0 si v esta libre
};
/**
 * Uso[0]: MaxClique mc(n); mc.add_edge(u, v); mc.solve();
 * Clique maxima de un grafo con hasta 60 nodos usando bitsets.
 * Complejidad: exponencial en el peor caso, practico para $N \leq 60$.
 */
struct MaxClique {
  int n, max_sz;
  vector<ll> adj;
  MaxClique(int n) : n(n), max_sz(0), adj(n, 0) {}
  void add_edge(int u, int v) { adj[u] |= (1LL << v); adj[v] |= (1LL << u); }
  void dfs(ll R, ll P, ll X) {
    if (!P && !X) { max_sz = max(max_sz, __builtin_popcountll(R)); return; }
    ll iter = P & ~adj[__builtin_ctzll(P | X)];
    while (iter) {
      int v = __builtin_ctzll(iter); ll v_bit = 1LL << v;
      dfs(R | v_bit, P & adj[v], X & adj[v]);
      P &= ~v_bit; X |= v_bit; iter &= ~v_bit;
    }
  }
  int solve() { dfs(0, (1LL << n) - 1, 0); return max_sz; }
};
/**
 * Uso[1]: BridgeTree bt(n, num_aristas); bt.add_edge(u, v, id); bt.build();
 * Arbol de puentes: detecta puentes (aristas cuya eliminacion desconecta el grafo) y puntos de articulacion, y contrae cada componente 2-arista-conexa en un nodo del arbol resultante (tree_adj).
 * Indexación: nodos 1-indexados; aristas identificadas por su id.
 * Complejidad: $O(V + E)$.
 */
struct BridgeTree {
  int n, timer, num_comps;
  vector<vector<pair<int, int>>> adj;
  vi tin, low, comp_id;
  vector<bool> is_bridge, is_articulation;
  vector<vi> tree_adj;
  BridgeTree(int n, int num_edges) : n(n), adj(n + 1), tin(n + 1, -1), low(n + 1, -1),
    comp_id(n + 1, 0), is_bridge(num_edges, false), is_articulation(n + 1, false), timer(0), num_comps(0) {}
  void add_edge(int u, int v, int id) { adj[u].pb({v, id}); adj[v].pb({u, id}); }
  void dfs_tarjan(int u, int p = -1) {
    tin[u] = low[u] = ++timer;
    int children = 0;
    for (auto [v, id] : adj[u]) {
      if (v == p) continue;
      if (tin[v] != -1) low[u] = min(low[u], tin[v]);
      else {
        children++;
        dfs_tarjan(v, u);
        low[u] = min(low[u], low[v]);
        if (low[v] > tin[u]) is_bridge[id] = true;
        if (low[v] >= tin[u] && p != -1) is_articulation[u] = true;
      }
    }
    if (p == -1 && children > 1) is_articulation[u] = true;
  }
  void dfs_comp(int u, int cur) {
    comp_id[u] = cur;
    for (auto [v, id] : adj[u]) if (!comp_id[v]) {
      if (is_bridge[id]) {
        tree_adj.pb(vi());
        tree_adj[cur].pb(++num_comps);
        tree_adj[num_comps].pb(cur);
        dfs_comp(v, num_comps);
      } else dfs_comp(v, cur);
    }
  }
  void build() {
    for (int i = 1; i <= n; i++) if (tin[i] == -1) dfs_tarjan(i);
    tree_adj.pb(vi());
    for (int i = 1; i <= n; i++) if (!comp_id[i]) tree_adj.pb(vi()), dfs_comp(i, ++num_comps);
  }
};

/**
 * Uso[0]: Dinic d(n, s, t); d.add_edge(u, v, cap); d.max_flow(); d.restore_path();
 * Flujo maximo con el algoritmo de Dinic, dirigido o no dirigido; restore_path devuelve un camino $s \to t$ por el que circula flujo positivo.
 * Indexación: nodos 0-indexados.
 * Complejidad: $O(V^2 E)$ en general, $O(E \sqrt{V})$ en redes unitarias.
 */
struct Dinic {
  struct FlowEdge { int u, v; ll cap, flow = 0; };
  int n, s, t;
  vector<FlowEdge> edges;
  vector<vi> adj;
  vi level, ptr;
  Dinic(int n, int s, int t) : n(n), s(s), t(t), adj(n + 1), level(n + 1), ptr(n + 1) {}
  int add_edge(int u, int v, ll cap, bool directed = true) {
    adj[u].pb(edges.size()); edges.pb({u, v, cap});
    adj[v].pb(edges.size()); edges.pb({v, u, directed ? 0 : cap});
    return (int)edges.size() - 2; // indice de la arista forward
  }
  bool bfs() {
    fill(all(level), -1);
    level[s] = 0;
    queue<int> q;
    q.push(s);
    while (!q.empty()) {
      int u = q.front(); q.pop();
      for (int id : adj[u]) {
        auto& e = edges[id];
        if (e.cap - e.flow >= 1 && level[e.v] == -1) level[e.v] = level[u] + 1, q.push(e.v);
      }
    }
    return level[t] != -1;
  }
  ll dfs(int u, ll pushed) {
    if (!pushed) return 0;
    if (u == t) return pushed;
    for (int& cid = ptr[u]; cid < adj[u].size(); cid++) {
      auto& e = edges[adj[u][cid]];
      ll tr = e.cap - e.flow;
      if (level[u] + 1 != level[e.v] || !tr) continue;
      ll push = dfs(e.v, min(pushed, tr));
      if (push) {
        e.flow += push;
        edges[adj[u][cid] ^ 1].flow -= push;
        return push;
      }
    }
    return 0;
  }
  ll max_flow() { return max_flow(s, t); }
  ll max_flow(int src, int snk) {
    s = src; t = snk;
    ll flow = 0;
    while (bfs()) {
      fill(all(ptr), 0);
      while (ll pushed = dfs(s, LINF)) flow += pushed;
    }
    return flow;
  }
  vi restore_path() { // camino s->t por aristas con flujo positivo (vacio si no llega)
    vi path = {s};
    int u = s;
    while (u != t) {
      bool advanced = false;
      for (int id : adj[u]) if (edges[id].flow > 0) {
        u = edges[id].v; path.pb(u); advanced = true; break;
      }
      if (!advanced) return {};
    }
    return path;
  }
};

/**
 * Uso[1]: DinicLB d(n); int id = d.add_edge(u, v, lo, hi); ll flujo = d.max_flow(s, t); ll f = d.edge_flow(id);
 * Flujo maximo con cotas inferiores (lower bounds) en cada arista; devuelve -1 si no existe flujo factible.
 * Complejidad: $O(\text{Dinic})$ sobre el grafo transformado (dos fases de flujo maximo).
 */
struct DinicLB {
  int n, ss, tt;
  Dinic din;
  vector<ll> bal, lo;
  vector<int> edge_idx;
  bool built = false;
  ll total = 0;
  DinicLB(int n) : n(n), ss(0), tt(n + 1), din(n + 2, 0, 0), bal(n + 2, 0) {}
  int add_edge(int u, int v, ll lower, ll upper) {
    bal[u] -= lower; bal[v] += lower;
    lo.pb(lower);
    int id = din.add_edge(u, v, upper - lower);
    edge_idx.pb(id);
    return (int)lo.size() - 1;
  }
  void build_super() {
    if (built) return;
    built = true;
    for (int i = 1; i <= n; i++) {
      if (bal[i] > 0) din.add_edge(ss, i, bal[i]), total += bal[i];
      else if (bal[i] < 0) din.add_edge(i, tt, -bal[i]);
    }
  }
  bool feasible() {
    build_super();
    return din.max_flow(ss, tt) == total;
  }
  ll max_flow(int s, int t) {
    int idx = din.add_edge(t, s, LINF);
    if (!feasible()) return -1;
    ll flow = din.edges[idx].flow; // flujo de la circulacion (valor base s->t)
    // elimina t->s (cap y flow en ambas direcciones, si no queda flujo residual cancelable)
    din.edges[idx].flow = 0; din.edges[idx].cap = 0;
    din.edges[idx ^ 1].flow = 0; din.edges[idx ^ 1].cap = 0;
    return flow + din.max_flow(s, t);
  }
  ll edge_flow(int id) const { return lo[id] + din.edges[edge_idx[id]].flow; }
};

/**
 * Uso[0]: MCMF m(n); m.add_edge(u, v, cap, costo); auto [flujo, costo] = m.solve(s, t);
 * Flujo maximo de costo minimo (min cost max flow) con potenciales para costos negativos.
 * Complejidad: $O(F \cdot E \log V)$ con $F$ el flujo total.
 */
struct MCMF {
  struct Edge { int to, rev; ll cap, flow, cost; };
  int n;
  vector<vector<Edge>> adj;
  vll dist, pot;
  vi p_node, p_edge;
  MCMF(int n) : n(n), adj(n + 1), dist(n + 1), pot(n + 1, 0), p_node(n + 1), p_edge(n + 1) {}
  void add_edge(int u, int v, ll cap, ll cost) {
    adj[u].pb({v, (int)adj[v].size(), cap, 0, cost});
    adj[v].pb({u, (int)adj[u].size() - 1, 0, 0, -cost});
  }
  // Inicializa potenciales si hay costos negativos. Omitir si todos son >= 0.
  void init_potentials(int s) {
    fill(all(pot), LINF);
    pot[s] = 0;
    for (int i = 0; i < n; i++)
      for (int u = 0; u <= n; u++) if (pot[u] != LINF)
        for (auto& e : adj[u]) if (e.cap - e.flow > 0 && pot[e.to] > pot[u] + e.cost) pot[e.to] = pot[u] + e.cost;
  }
  bool dijkstra(int s, int t) {
    fill(all(dist), LINF);
    dist[s] = 0;
    priority_queue<pair<ll, int>, vector<pair<ll, int>>, greater<>> pq;
    pq.push({0, s});
    while (!pq.empty()) {
      auto [d, u] = pq.top(); pq.pop();
      if (d > dist[u]) continue;
      for (int i = 0; i < adj[u].size(); i++) {
        auto& e = adj[u][i];
        ll rc = e.cost + pot[u] - pot[e.to]; // costo reducido
        if (e.cap - e.flow > 0 && dist[e.to] > dist[u] + rc) {
          dist[e.to] = dist[u] + rc;
          p_node[e.to] = u; p_edge[e.to] = i;
          pq.push({dist[e.to], e.to});
        }
      }
    }
    return dist[t] != LINF;
  }
  pair<ll, ll> solve(int s, int t) { // {flujo, costo}
    init_potentials(s);
    ll flow = 0, cost = 0;
    while (dijkstra(s, t)) {
      for (int i = 0; i <= n; i++) if (dist[i] != LINF) pot[i] += dist[i];
      ll push = LINF;
      for (int cur = t; cur != s; cur = p_node[cur])
        push = min(push, adj[p_node[cur]][p_edge[cur]].cap - adj[p_node[cur]][p_edge[cur]].flow);
      for (int cur = t; cur != s; cur = p_node[cur]) {
        auto& e = adj[p_node[cur]][p_edge[cur]];
        e.flow += push;
        adj[cur][e.rev].flow -= push;
        cost += push * e.cost;
      }
      flow += push;
    }
    return {flow, cost};
  }
};

/**
 * Uso[1]: DominatorTree dt(n); dt.add_edge(u, v); dt.build(root);
 * Arbol dominador de un grafo dirigido desde una raiz; dom[u] son los hijos de u en el arbol.
 * Complejidad: $O((V + E) \log V)$.
 */
struct DominatorTree {
  int n, t;
  vector<vi> adj, rev, dom;
  vi dfn, id, p, sdom, idom, dsu, best;
  DominatorTree(int n) : n(n), t(0), adj(n + 1), rev(n + 1), dom(n + 1),
    dfn(n + 1, 0), id(n + 1), p(n + 1), sdom(n + 1), idom(n + 1), dsu(n + 1), best(n + 1) {}
  void add_edge(int u, int v) { adj[u].pb(v); }
  void dfs(int u) {
    dfn[u] = ++t;
    id[t] = u;
    sdom[t] = best[t] = dsu[t] = t;
    for (int v : adj[u]) {
      if (!dfn[v]) dfs(v), p[dfn[v]] = dfn[u];
      rev[dfn[v]].pb(dfn[u]);
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
    vector<vi> bucket(n + 1);
    for (int i = t; i >= 2; i--) {
      for (int u : rev[i]) {
        find(u);
        if (sdom[best[u]] < sdom[i]) sdom[i] = sdom[best[u]];
      }
      bucket[sdom[i]].pb(i);
      dsu[i] = p[i];
      for (int u : bucket[p[i]]) {
        find(u);
        idom[u] = sdom[best[u]] == sdom[u] ? sdom[u] : best[u];
      }
      bucket[p[i]].clear();
    }
    for (int i = 2; i <= t; i++) {
      if (idom[i] != sdom[i]) idom[i] = idom[idom[i]];
      dom[id[idom[i]]].pb(id[i]);
    }
  }
};
