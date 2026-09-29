/**
 * Uso[1]: TreeEdge<T> e = {to, id, weight};
 * Arista de arbol con destino, identificador opcional y peso generico.
 * Complejidad: $O(1)$ por arista.
 */
template <typename T = ll> struct TreeEdge { int to, id; T weight; };

/**
 * Uso[1]: Tree<ll> t(n); t.add_edge(u, v, w); t.init(root); t.get_diameter(); t.get_centers(); t.all_farthest_distances();
 * Arbol con raiz: tiempos de entrada/salida (euler), tamanos de subarbol, diametro con camino y centros.
 * Complejidad: init y diametro $O(N)$.
 */
template <typename T = ll> struct Tree {
  int n, root, timer;
  vector<vector<TreeEdge<T>>> adj;
  vi parent, depth, sz, tin, tout, euler_order;
  vector<T> dist_from_root;
  Tree(int n = 0) : n(n), root(1), timer(0), adj(n + 1), parent(n + 1), depth(n + 1),
    sz(n + 1), tin(n + 1), tout(n + 1), dist_from_root(n + 1) {}
  void add_edge(int u, int v, T w = 1, int id = -1) {
    adj[u].pb({v, id, w});
    adj[v].pb({u, id, w});
  }
  void add_directed_edge(int u, int v, T w = 1, int id = -1) { adj[u].pb({v, id, w}); }
  void dfs_init(int u, int p = 0, int d = 0, T cur = 0) {
    parent[u] = p; depth[u] = d; dist_from_root[u] = cur;
    sz[u] = 1; tin[u] = ++timer;
    euler_order.pb(u);
    for (auto& e : adj[u]) if (e.to != p) {
      dfs_init(e.to, u, d + 1, cur + e.weight);
      sz[u] += sz[e.to];
    }
    tout[u] = timer;
  }
  void init(int root = 1) {
    this->root = root; timer = 0;
    euler_order.clear();
    dfs_init(root, 0, 0, 0);
  }
  struct DiameterResult { T length; int u, v; vi path; };
  DiameterResult get_diameter() {
    auto farthest = [&](int start) {
      vector<T> dist(n + 1, -1);
      queue<int> q;
      q.push(start); dist[start] = 0;
      int far = start;
      while (!q.empty()) {
        int u = q.front(); q.pop();
        if (dist[u] > dist[far]) far = u;
        for (auto& e : adj[u]) if (dist[e.to] == -1) dist[e.to] = dist[u] + e.weight, q.push(e.to);
      }
      return pair<int, T>(far, dist[far]);
    };
    int u = farthest(1).fi;
    auto [v, length] = farthest(u);
    vi prev(n + 1, -1);
    queue<int> q;
    q.push(u); prev[u] = 0;
    while (!q.empty()) {
      int cur = q.front(); q.pop();
      if (cur == v) break;
      for (auto& e : adj[cur]) if (prev[e.to] == -1) prev[e.to] = cur, q.push(e.to);
    }
    vi path;
    for (int cur = v; cur != 0; cur = prev[cur]) path.pb(cur);
    reverse(all(path));
    return {length, u, v, path};
  }
  // 1 o 2 nodos centrales del diametro
  vi get_centers() {
    auto d = get_diameter();
    int m = d.path.size();
    if (m & 1) return {d.path[m / 2]};
    return {d.path[m / 2 - 1], d.path[m / 2]};
  }
  // Distancia maxima desde cada nodo a cualquier otro. O(N)
  vector<T> all_farthest_distances() {
    auto d = get_diameter();
    auto dist_from = [&](int start) {
      vector<T> dist(n + 1, -1);
      queue<int> q;
      q.push(start); dist[start] = 0;
      while (!q.empty()) {
        int cur = q.front(); q.pop();
        for (auto& e : adj[cur]) if (dist[e.to] == -1) dist[e.to] = dist[cur] + e.weight, q.push(e.to);
      }
      return dist;
    };
    auto du = dist_from(d.u), dv = dist_from(d.v);
    vector<T> res(n + 1);
    for (int i = 1; i <= n; i++) res[i] = max(du[i], dv[i]);
    return res;
  }
};

/**
 * Uso[1]: BinaryLifting<ll> bl(t, root); bl.lca(u, v); bl.kth_ancestor(u, k); bl.path_sum(u, v); bl.dist_weighted(u, v);
 * Binary lifting sobre un arbol: lca, k-esimo ancestro, suma/max/min de aristas en el camino y distancias.
 * Complejidad: preproceso $O(N \log N)$, consulta $O(\log N)$.
 */
template <typename T = ll> struct BinaryLifting {
  int n, log_n;
  vi depth;
  vector<vi> up;
  vector<vector<T>> mx, mn, sum;
  vector<T> dist_root;
  BinaryLifting(const Tree<T>& tree, int root = 1) {
    n = tree.n;
    log_n = 32 - __builtin_clz(n) + 1;
    depth = tree.depth;
    dist_root = tree.dist_from_root;
    up.assign(n + 1, vi(log_n, 0));
    mx.assign(n + 1, vector<T>(log_n, -LINF));
    mn.assign(n + 1, vector<T>(log_n, LINF));
    sum.assign(n + 1, vector<T>(log_n, 0));
    for (int u = 1; u <= n; u++) up[u][0] = tree.parent[u] ? tree.parent[u] : u;
    for (int u = 1; u <= n; u++)
      for (auto& e : tree.adj[u]) if (e.to == tree.parent[u]) mx[u][0] = mn[u][0] = sum[u][0] = e.weight;
    for (int j = 1; j < log_n; j++)
      for (int i = 1; i <= n; i++) {
        int p = up[i][j - 1];
        up[i][j] = up[p][j - 1];
        mx[i][j] = max(mx[i][j - 1], mx[p][j - 1]);
        mn[i][j] = min(mn[i][j - 1], mn[p][j - 1]);
        sum[i][j] = sum[i][j - 1] + sum[p][j - 1];
      }
  }
  int kth_ancestor(int u, int k) {
    if (depth[u] < k) return 0;
    for (int j = 0; j < log_n; j++) if (k & (1 << j)) u = up[u][j];
    return u;
  }
  int lca(int u, int v) {
    if (depth[u] < depth[v]) swap(u, v);
    u = kth_ancestor(u, depth[u] - depth[v]);
    if (u == v) return u;
    for (int j = log_n - 1; j >= 0; j--) if (up[u][j] != up[v][j]) u = up[u][j], v = up[v][j];
    return up[u][0];
  }
  int dist_edges(int u, int v) { return depth[u] + depth[v] - 2 * depth[lca(u, v)]; }
  T dist_weighted(int u, int v) { return dist_root[u] + dist_root[v] - 2 * dist_root[lca(u, v)]; }
  // k-esimo nodo en el camino u->v (0-indexed: 0 = u)
  int kth_node_on_path(int u, int v, int k) {
    int anc = lca(u, v);
    int d1 = depth[u] - depth[anc], d2 = depth[v] - depth[anc];
    if (k <= d1) return kth_ancestor(u, k);
    k -= d1;
    if (k <= d2) return kth_ancestor(v, d2 - k);
    return 0;
  }
  T path_max(int u, int v) {
    int anc = lca(u, v);
    T res = -LINF;
    auto lift = [&](int node, int steps) {
      for (int j = 0; j < log_n; j++) if (steps & (1 << j)) res = max(res, mx[node][j]), node = up[node][j];
    };
    lift(u, depth[u] - depth[anc]);
    lift(v, depth[v] - depth[anc]);
    return res;
  }
  T path_min(int u, int v) {
    int anc = lca(u, v);
    T res = LINF;
    auto lift = [&](int node, int steps) {
      for (int j = 0; j < log_n; j++) if (steps & (1 << j)) res = min(res, mn[node][j]), node = up[node][j];
    };
    lift(u, depth[u] - depth[anc]);
    lift(v, depth[v] - depth[anc]);
    return res;
  }
  T path_sum(int u, int v) {
    int anc = lca(u, v);
    T res = 0;
    auto lift = [&](int node, int steps) {
      for (int j = 0; j < log_n; j++) if (steps & (1 << j)) res += sum[node][j], node = up[node][j];
    };
    lift(u, depth[u] - depth[anc]);
    lift(v, depth[v] - depth[anc]);
    return res;
  }
  bool is_ancestor(int u, int v) { return lca(u, v) == u; }
};

/**
 * Uso[1]: EulerTourTree<ll> et(t); et.subtree_range(u); et.subtree_size(u);
 * Mapea el subarbol de cada nodo a un rango contiguo de la secuencia de entrada (tin/tout).
 * Complejidad: $O(1)$ por consulta.
 */
template <typename T = ll> struct EulerTourTree {
  vi in, out;
  EulerTourTree(const Tree<T>& tree) : in(tree.tin), out(tree.tout) {}
  pair<int, int> subtree_range(int u) { return {in[u], out[u]}; }
  int subtree_size(int u) { return out[u] - in[u] + 1; }
};

/**
 * Uso[1]: HLD<ll> hld(n); hld.add_edge(u, v); hld.init(root); hld.update_path(u, v, val); hld.query_path(u, v);
 * Heavy-light decomposition con segment tree perezoso para sumas y updates sobre caminos y subarboles.
 * Complejidad: preproceso $O(N)$, update y query de camino $O(\log^2 N)$.
 */
template <typename T = ll> struct HLD {
  struct SegTree {
    int size;
    vector<T> tree, lazy;
    void init(int n) {
      size = 1;
      while (size < n) size <<= 1;
      tree.assign(2 * size, 0);
      lazy.assign(2 * size, 0);
    }
    void push(int node, int l, int r) {
      if (!lazy[node]) return;
      tree[node] += lazy[node] * (r - l + 1);
      if (l != r) lazy[2 * node] += lazy[node], lazy[2 * node + 1] += lazy[node];
      lazy[node] = 0;
    }
    void update(int node, int l, int r, int ql, int qr, T val) {
      push(node, l, r);
      if (ql > r || qr < l) return;
      if (ql <= l && r <= qr) { lazy[node] += val; push(node, l, r); return; }
      int mid = (l + r) / 2;
      update(2 * node, l, mid, ql, qr, val);
      update(2 * node + 1, mid + 1, r, ql, qr, val);
      tree[node] = tree[2 * node] + tree[2 * node + 1];
    }
    T query(int node, int l, int r, int ql, int qr) {
      push(node, l, r);
      if (ql > r || qr < l) return 0;
      if (ql <= l && r <= qr) return tree[node];
      int mid = (l + r) / 2;
      return query(2 * node, l, mid, ql, qr) + query(2 * node + 1, mid + 1, r, ql, qr);
    }
  } seg;
  int n, cur_pos;
  vector<vi> adj;
  vi parent, depth, heavy, head, pos;
  HLD(int n) : n(n), cur_pos(1), adj(n + 1), parent(n + 1), depth(n + 1),
    heavy(n + 1, -1), head(n + 1), pos(n + 1) { seg.init(n + 1); }
  void add_edge(int u, int v) { adj[u].pb(v); adj[v].pb(u); }
  int dfs_size(int u, int p, int d) {
    int sz = 1, max_c = 0;
    depth[u] = d; parent[u] = p; heavy[u] = -1;
    for (int v : adj[u]) if (v != p) {
      int cs = dfs_size(v, u, d + 1);
      sz += cs;
      if (cs > max_c) max_c = cs, heavy[u] = v;
    }
    return sz;
  }
  void decompose(int u, int h) {
    head[u] = h;
    pos[u] = cur_pos++;
    if (heavy[u] != -1) decompose(heavy[u], h);
    for (int v : adj[u]) if (v != parent[u] && v != heavy[u]) decompose(v, v);
  }
  void init(int root = 1) { cur_pos = 1; dfs_size(root, 0, 0); decompose(root, root); }
  void update_path(int u, int v, T val) {
    while (head[u] != head[v]) {
      if (depth[head[u]] > depth[head[v]]) swap(u, v);
      seg.update(1, 1, n, pos[head[v]], pos[v], val);
      v = parent[head[v]];
    }
    if (depth[u] > depth[v]) swap(u, v);
    seg.update(1, 1, n, pos[u], pos[v], val);
  }
  T query_path(int u, int v) {
    T res = 0;
    while (head[u] != head[v]) {
      if (depth[head[u]] > depth[head[v]]) swap(u, v);
      res += seg.query(1, 1, n, pos[head[v]], pos[v]);
      v = parent[head[v]];
    }
    if (depth[u] > depth[v]) swap(u, v);
    return res + seg.query(1, 1, n, pos[u], pos[v]);
  }
  void update_subtree(int u, int sz_u, T val) { seg.update(1, 1, n, pos[u], pos[u] + sz_u - 1, val); }
  T query_subtree(int u, int sz_u) { return seg.query(1, 1, n, pos[u], pos[u] + sz_u - 1); }
  void update_node(int u, T val) { seg.update(1, 1, n, pos[u], pos[u], val); }
  int lca(int u, int v) {
    while (head[u] != head[v]) {
      if (depth[head[u]] > depth[head[v]]) swap(u, v);
      v = parent[head[v]];
    }
    return depth[u] < depth[v] ? u : v;
  }
};

/**
 * Uso[1]: CentroidDecomposition cd(n); cd.add_edge(u, v, w); cd.K = lim; cd.init(true); cd.total_paths;
 * Descomposicion en centroides: conteo offline de pares con distancia menor o igual a $K$ y consultas online de nodo mas cercano.
 * Complejidad: construccion $O(N \log N)$, update y query online $O(\log N)$.
 */
struct CentroidDecomposition {
  int n;
  vector<vector<pair<int, ll>>> adj;
  vi sz, c_parent;
  vector<bool> removed;
  ll total_paths = 0, K = 0;
  vll local_ans;
  CentroidDecomposition(int n) : n(n), adj(n + 1), sz(n + 1), c_parent(n + 1),
    removed(n + 1, false), local_ans(n + 1, LINF) {}
  void add_edge(int u, int v, ll w = 1) { adj[u].pb({v, w}); adj[v].pb({u, w}); }
  int get_sizes(int u, int p = 0) {
    sz[u] = 1;
    for (auto [v, w] : adj[u]) if (v != p && !removed[v]) sz[u] += get_sizes(v, u);
    return sz[u];
  }
  int get_centroid(int u, int p, int size) {
    for (auto [v, w] : adj[u]) if (v != p && !removed[v] && sz[v] > size / 2) return get_centroid(v, u, size);
    return u;
  }
  void get_paths(int u, int p, ll dist, vll& paths) {
    paths.pb(dist);
    for (auto [v, w] : adj[u]) if (v != p && !removed[v]) get_paths(v, u, dist + w, paths);
  }
  ll count_pairs(vll& paths, ll limit) {
    sort(all(paths));
    ll cnt = 0;
    for (int l = 0, r = paths.size() - 1; l < r;) {
      if (paths[l] + paths[r] <= limit) cnt += r - l, l++;
      else r--;
    }
    return cnt;
  }
  void process_centroid(int c) {
    vll all_paths{0};
    for (auto [v, w] : adj[c]) if (!removed[v]) get_paths(v, c, w, all_paths);
    total_paths += count_pairs(all_paths, K);
    for (auto [v, w] : adj[c]) if (!removed[v]) {
      vll sub;
      get_paths(v, c, w, sub);
      total_paths -= count_pairs(sub, K);
    }
  }
  int build_tree(int u, int p = 0, bool offline_mode = false) {
    int c = get_centroid(u, 0, get_sizes(u, 0));
    if (offline_mode) process_centroid(c);
    removed[c] = true;
    c_parent[c] = p;
    for (auto [v, w] : adj[c]) if (!removed[v]) build_tree(v, c, offline_mode);
    return c;
  }
  void init(bool offline_mode = false) { build_tree(1, 0, offline_mode); }
  // Pinta u de rojo. LCA del arbol original.
  void update(int u, LCA& lca) {
    for (int cur = u; cur; cur = c_parent[cur]) local_ans[cur] = min(local_ans[cur], (ll)lca.dist(u, cur));
  }
  ll query(int u, LCA& lca) {
    ll best = LINF;
    for (int cur = u; cur; cur = c_parent[cur]) best = min(best, (ll)lca.dist(u, cur) + local_ans[cur]);
    return best;
  }
};

/**
 * Uso[1]: DSUOnTree s(n); s.add_edge(u, v); s.set_colors(color_1index); s.solve(root); s.ans[u];
 * Small-to-large (sack) para responder estadisticas de colores en cada subarbol.
 * Complejidad: $O(N \log N)$.
 */
struct DSUOnTree {
  int n, timer, distinct;
  vector<vi> adj;
  vi sz, heavy, tin, tout, euler_node, color, ans, cnt;
  DSUOnTree(int n) : n(n), timer(0), distinct(0), adj(n + 1), sz(n + 1), heavy(n + 1, -1),
    tin(n + 1), tout(n + 1), euler_node(n + 1), color(n + 1), ans(n + 1), cnt(n + 1) {}
  void add_edge(int u, int v) { adj[u].pb(v); adj[v].pb(u); }
  void set_colors(const vi& c) {
    for (int i = 1; i <= n; i++) color[i] = c[i];
  }
  void dfs_size(int u, int p = 0) {
    sz[u] = 1;
    tin[u] = ++timer;
    euler_node[timer] = u;
    int max_sz = 0;
    for (int v : adj[u]) if (v != p) {
      dfs_size(v, u);
      sz[u] += sz[v];
      if (sz[v] > max_sz) max_sz = sz[v], heavy[u] = v;
    }
    tout[u] = timer;
  }
  void add_node(int u) { if (!cnt[color[u]]++) distinct++; }
  void remove_node(int u) { if (!--cnt[color[u]]) distinct--; }
  void dfs_sack(int u, int p, bool keep) {
    for (int v : adj[u]) if (v != p && v != heavy[u]) dfs_sack(v, u, false);
    if (heavy[u] != -1) dfs_sack(heavy[u], u, true);
    add_node(u);
    for (int v : adj[u]) if (v != p && v != heavy[u])
      for (int t = tin[v]; t <= tout[v]; t++) add_node(euler_node[t]);
    ans[u] = distinct;
    if (!keep) for (int t = tin[u]; t <= tout[u]; t++) remove_node(euler_node[t]);
  }
  void solve(int root = 1) { timer = 0; dfs_size(root); dfs_sack(root, 0, false); }
};

/**
 * Uso[1]: VirtualTree vt(t, bl); auto [raiz, nodos] = vt.build(vector_de_nodos);
 * Construye el arbol virtual sobre un subconjunto de nodos, con aristas pesadas por distancia.
 * Complejidad: $O(K \log N)$ con $K$ el numero de nodos solicitados.
 */
struct VirtualTree {
  const Tree<ll>& tree;
  BinaryLifting<ll>& bl;
  vector<vector<pair<int, ll>>> virt_adj;
  VirtualTree(const Tree<ll>& tree, BinaryLifting<ll>& bl) : tree(tree), bl(bl), virt_adj(tree.n + 1) {}
  // Retorna {raiz, nodos del arbol virtual}
  pair<int, vi> build(vi nodes) {
    sort(all(nodes), [&](int a, int b) { return tree.tin[a] < tree.tin[b]; });
    int k = nodes.size();
    for (int i = 0; i + 1 < k; i++) nodes.pb(bl.lca(nodes[i], nodes[i + 1]));
    sort(all(nodes), [&](int a, int b) { return tree.tin[a] < tree.tin[b]; });
    nodes.erase(unique(all(nodes)), nodes.end());
    for (int u : nodes) virt_adj[u].clear();
    vi stk{nodes[0]};
    for (int i = 1; i < nodes.size(); i++) {
      int u = nodes[i];
      while (!stk.empty() && !bl.is_ancestor(stk.back(), u)) stk.pop_back();
      if (!stk.empty()) {
        ll w = bl.dist_weighted(stk.back(), u);
        virt_adj[stk.back()].pb({u, w});
        virt_adj[u].pb({stk.back(), w});
      }
      stk.pb(u);
    }
    return {nodes[0], nodes};
  }
};

/**
 * Uso[1]: TreeHashing::are_isomorphic(n1, adj1, n2, adj2);
 * Determina si dos arboles son isomorfos (sin considerar etiquetas) mediante hashing y centroides.
 * Complejidad: $O(N)$.
 */
struct TreeHashing {
  static uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
    x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
    return x ^ (x >> 31);
  }
  static uint64_t hash_node(uint64_t x) { return splitmix64(x + 0x517cc1b727220a95); }
  static uint64_t rooted_hash(int u, int p, const vector<vi>& adj) {
    vector<uint64_t> child;
    for (int v : adj[u]) if (v != p) child.pb(rooted_hash(v, u, adj));
    sort(all(child));
    uint64_t h = 1;
    for (uint64_t c : child) h += hash_node(c);
    return splitmix64(h);
  }
  static vector<int> get_centers(int n, const vector<vi>& adj) {
    vi deg(n + 1);
    queue<int> leaves;
    for (int i = 1; i <= n; i++) {
      deg[i] = adj[i].size();
      if (deg[i] <= 1) leaves.push(i);
    }
    int remaining = n;
    while (remaining > 2) {
      int sz = leaves.size();
      remaining -= sz;
      for (int i = 0; i < sz; i++) {
        int u = leaves.front(); leaves.pop();
        for (int v : adj[u]) if (--deg[v] == 1) leaves.push(v);
      }
    }
    vi centers;
    while (!leaves.empty()) centers.pb(leaves.front()), leaves.pop();
    return centers;
  }
  static bool are_isomorphic(int n1, const vector<vi>& adj1, int n2, const vector<vi>& adj2) {
    if (n1 != n2) return false;
    auto c1 = get_centers(n1, adj1), c2 = get_centers(n2, adj2);
    for (int r1 : c1)
      for (int r2 : c2)
        if (rooted_hash(r1, 0, adj1) == rooted_hash(r2, 0, adj2)) return true;
    return false;
  }
};

/**
 * Uso[1]: RerootingDP r(n); r.add_edge(u, v); auto ans = r.solve(root);
 * DP con rerooting para calcular la solucion para todas las posibles raices; ajustar merge, extend y finalize.
 * Complejidad: $O(N)$.
 */
struct RerootingDP {
  struct State { ll val; };
  int n;
  vector<vi> adj;
  vector<State> dp_down, ans;
  const State NEUTRAL = {0};
  State merge(State a, State b) { return {max(a.val, b.val)}; }
  State extend(State s, int v, int u) { return {s.val + 1}; } // transicion por arista
  State finalize(State s, int u) { return s; }
  RerootingDP(int n) : n(n), adj(n + 1), dp_down(n + 1, NEUTRAL), ans(n + 1, NEUTRAL) {}
  void add_edge(int u, int v) { adj[u].pb(v); adj[v].pb(u); }
  void dfs_down(int u, int p) {
    State cur = NEUTRAL;
    for (int v : adj[u]) if (v != p) dfs_down(v, u), cur = merge(cur, extend(dp_down[v], v, u));
    dp_down[u] = finalize(cur, u);
  }
  void dfs_reroot(int u, int p, State from_parent) {
    vi children;
    for (int v : adj[u]) if (v != p) children.pb(v);
    int k = children.size();
    vector<State> pref(k + 1, NEUTRAL), suff(k + 1, NEUTRAL);
    for (int i = 0; i < k; i++) pref[i + 1] = merge(pref[i], extend(dp_down[children[i]], children[i], u));
    for (int i = k - 1; i >= 0; i--) suff[i] = merge(suff[i + 1], extend(dp_down[children[i]], children[i], u));
    ans[u] = finalize(merge(pref[k], from_parent), u);
    for (int i = 0; i < k; i++) {
      State outside = merge(merge(pref[i], suff[i + 1]), from_parent);
      dfs_reroot(children[i], u, extend(finalize(outside, u), u, children[i]));
    }
  }
  vector<State> solve(int root = 1) {
    dfs_down(root, 0);
    dfs_reroot(root, 0, NEUTRAL);
    return ans;
  }
};

/**
 * Uso[1]: FastLCA<ll> fl(t); fl.lca(u, v); fl.dist_edges(u, v);
 * LCA por Euler + RMQ con respuesta en tiempo constante por consulta.
 * Complejidad: preproceso $O(N \log N)$, consulta $O(1)$.
 */
template <typename T = ll> struct FastLCA {
  const Tree<T>& tree;
  vi tin, euler, euler_depth, log_table;
  vector<vi> st;
  FastLCA(const Tree<T>& tree) : tree(tree), tin(tree.n + 1) {
    auto dfs = [&](auto&& self, int u, int p, int d) -> void {
      tin[u] = euler.size();
      euler.pb(u);
      euler_depth.pb(d);
      for (auto& e : tree.adj[u]) if (e.to != p) {
        self(self, e.to, u, d + 1);
        euler.pb(u);
        euler_depth.pb(d);
      }
    };
    dfs(dfs, tree.root, 0, 0);
    int m = euler.size();
    log_table.assign(m + 1, 0);
    for (int i = 2; i <= m; i++) log_table[i] = log_table[i / 2] + 1;
    int K = log_table[m] + 1;
    st.assign(K, vi(m));
    for (int i = 0; i < m; i++) st[0][i] = i;
    for (int i = 1; i < K; i++)
      for (int j = 0; j + (1 << i) <= m; j++) {
        int l = st[i - 1][j], r = st[i - 1][j + (1 << (i - 1))];
        st[i][j] = euler_depth[l] < euler_depth[r] ? l : r;
      }
  }
  int lca(int u, int v) {
    int l = tin[u], r = tin[v];
    if (l > r) swap(l, r);
    int i = log_table[r - l + 1];
    int a = st[i][l], b = st[i][r - (1 << i) + 1];
    return euler_depth[a] < euler_depth[b] ? euler[a] : euler[b];
  }
  int dist_edges(int u, int v) { return tree.depth[u] + tree.depth[v] - 2 * tree.depth[lca(u, v)]; }
};

/**
 * Uso[1]: MoOnTrees<ll> mo(n); mo.dfs_euler(root, 0, t); mo.solve(queries, fast_lca);
 * Algoritmo de Mo sobre caminos de un arbol; ajustar toggle para actualizar la respuesta al agregar o quitar un nodo.
 * Complejidad: $O((N + Q)\sqrt{N})$ toggle esperados.
 */
template <typename T = ll> struct MoOnTrees {
  struct Query {
    int l, r, lca, id, block_size;
    bool operator<(const Query& o) const {
      int b1 = l / block_size, b2 = o.l / block_size;
      if (b1 != b2) return b1 < b2;
      return (b1 & 1) ? (r < o.r) : (r > o.r);
    }
  };
  int n, timer;
  vi tin, tout, euler;
  vector<bool> in_path;
  ll current_ans = 0;
  MoOnTrees(int n) : n(n), timer(0), tin(n + 1), tout(n + 1), euler(2 * n + 1), in_path(n + 1, false) {}
  void dfs_euler(int u, int p, const Tree<T>& tree) {
    tin[u] = ++timer;
    euler[timer] = u;
    for (auto& e : tree.adj[u]) if (e.to != p) dfs_euler(e.to, u, tree);
    tout[u] = ++timer;
    euler[timer] = u;
  }
  void toggle(int u) {
    // TODO: actualizar respuesta al agregar/remover el nodo u (p.ej. frecuencias de colores)
    in_path[u] = !in_path[u];
  }
  vll solve(vector<pair<int, int>>& raw_queries, const FastLCA<T>& fast_lca) {
    int q = raw_queries.size();
    int bs = max(1, (int)(2 * n / sqrt(max(1, q))));
    vector<Query> queries;
    for (int i = 0; i < q; i++) {
      auto [u, v] = raw_queries[i];
      if (tin[u] > tin[v]) swap(u, v);
      int anc = fast_lca.lca(u, v);
      queries.pb(anc == u ? Query{tin[u], tin[v], 0, i, bs} : Query{tout[u], tin[v], anc, i, bs});
    }
    sort(all(queries));
    vll ans(q);
    int curr_l = 1, curr_r = 0;
    for (auto& mq : queries) {
      while (curr_l > mq.l) toggle(euler[--curr_l]);
      while (curr_r < mq.r) toggle(euler[++curr_r]);
      while (curr_l < mq.l) toggle(euler[curr_l++]);
      while (curr_r > mq.r) toggle(euler[curr_r--]);
      if (mq.lca) toggle(mq.lca);
      ans[mq.id] = current_ans;
      if (mq.lca) toggle(mq.lca);
    }
    return ans;
  }
};
