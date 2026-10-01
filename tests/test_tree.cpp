#include "test_base.hpp"
#include "graph.hpp"
#include "rq.hpp"
#include "tree.hpp"
using namespace std;
int failures = 0;
#define CHECK(cond) do { if (!(cond)) { cerr << "  FAIL L" << __LINE__ << ": " #cond "\n"; failures++; } } while (0)

int main() {
  // Tree base: camino 1-2-3-4-5
  Tree<ll> t(5);
  for (int i = 1; i < 5; i++) t.add_edge(i, i + 1);
  t.init(1);
  CHECK(t.sz[3] == 3);
  auto d = t.get_diameter();
  CHECK(d.length == 4);
  CHECK(t.get_centers() == (vi{3}));
  auto far = all_farthest_distances(t);
  CHECK(far[3] == 2 && far[1] == 4);
  // arbol de un solo nodo
  Tree<ll> single(1);
  single.init(1);
  CHECK(single.get_diameter().length == 0);

  // BinaryLifting con pesos
  Tree<ll> tw(4);
  tw.add_edge(1, 2, 5);
  tw.add_edge(2, 3, 2);
  tw.add_edge(2, 4, 3);
  tw.init(1);
  BinaryLifting<ll> bl(tw, 1);
  CHECK(bl.lca(3, 4) == 2);
  CHECK(bl.kth_ancestor(3, 1) == 2 && bl.kth_ancestor(3, 2) == 1);
  CHECK(bl.dist_edges(3, 4) == 2);
  CHECK(bl.dist_weighted(3, 4) == 5);
  CHECK(bl.path_sum(3, 4) == 5);
  CHECK(bl.path_max(3, 4) == 3 && bl.path_min(3, 4) == 2);
  CHECK(bl.is_ancestor(2, 3) && !bl.is_ancestor(3, 2));

  // EulerTourTree
  Tree<ll> te(4);
  te.add_edge(1, 2); te.add_edge(2, 3); te.add_edge(2, 4);
  te.init(1);
  EulerTourTree<ll> et(te);
  CHECK(et.subtree_size(2) == 3 && et.subtree_size(1) == 4);

  // HLD (nuevo API): pesos en nodos, caminos y subarboles con LazySegTree
  HLD hld(4);
  for (int i = 1; i < 4; i++) hld.add_edge(i, i + 1);
  hld.build(1);
  vll values(5, 0);
  values[hld.pos[2]] = 5;
  values[hld.pos[4]] = 3;
  LazySegTree st(values);
  auto path_sum = [&](int u, int v) {
    ll ans = 0;
    hld.process_path(u, v, [&](int l, int r) { ans += st.query(l, r); });
    return ans;
  };
  CHECK(path_sum(2, 4) == 8); // nodos 2+3+4 = 5+0+3
  CHECK(path_sum(1, 4) == 8);
  hld.process_path(1, 4, [&](int l, int r) { st.add_range(l, r, 1); });
  CHECK(path_sum(1, 4) == 12); // +1 en cada nodo del camino
  hld.process_subtree(2, [&](int l, int r) { st.add_range(l, r, 10); }); // subarbol 2: nodos 2,3,4
  CHECK(path_sum(2, 4) == 41); // 8 + 3 (add camino) + 30 (add subarbol)
  CHECK(path_sum(1, 4) == 42);

  // HLD con pesos en aristas (values_on_edges = true)
  HLD hld_edges(4, true);
  for (int i = 1; i < 4; i++) hld_edges.add_edge(i, i + 1);
  hld_edges.build(1);
  vll e_values(5, 0);
  e_values[hld_edges.pos[2]] = 2; // arista 1-2
  e_values[hld_edges.pos[3]] = 4; // arista 2-3
  e_values[hld_edges.pos[4]] = 1; // arista 3-4
  LazySegTree est(e_values);
  auto edge_path_sum = [&](int u, int v) {
    ll ans = 0;
    hld_edges.process_path(u, v, [&](int l, int r) { ans += est.query(l, r); });
    return ans;
  };
  CHECK(edge_path_sum(2, 4) == 5); // 4 + 1
  CHECK(edge_path_sum(1, 4) == 7); // 2 + 4 + 1
  CHECK(edge_path_sum(2, 2) == 0); // mismo nodo: sin aristas
  ll edge_sub = 0;
  hld_edges.process_subtree(2, [&](int l, int r) { edge_sub += est.query(l, r); }); // aristas 2-3 y 3-4
  CHECK(edge_sub == 5);

  // CentroidPathProcessing: contar pares a distancia exacta K
  auto brute_dist_pairs = [&](const vector<vector<pair<int, int>>>& ad, int n, int K) {
    int cnt = 0;
    for (int s = 1; s <= n; s++) {
      vector<int> d(n + 1, -1);
      queue<int> q; q.push(s); d[s] = 0;
      while (!q.empty()) {
        int u = q.front(); q.pop();
        for (auto [v, w] : ad[u]) if (d[v] == -1) d[v] = d[u] + 1, q.push(v);
      }
      for (int t = s + 1; t <= n; t++) if (d[t] == K) cnt++;
    }
    return cnt;
  };
  // camino 1-2-3-4-5, K=2 -> pares {1,3},{2,4},{3,5} = 3
  CentroidPathProcessing cpp(5);
  for (int i = 1; i < 5; i++) cpp.add_edge(i, i + 1);
  cpp.K = 2;
  cpp.solve(1);
  CHECK(cpp.global.valid_pairs == 3);
  // arboles aleatorios vs brute force
  {
    mt19937 rnd(123);
    for (int it = 0; it < 30; it++) {
      int n = 2 + rnd() % 12;
      CentroidPathProcessing cpx(n);
      vector<vector<pair<int, int>>> adx(n + 1);
      for (int i = 2; i <= n; i++) {
        int p = 1 + rnd() % (i - 1);
        cpx.add_edge(i, p, 1);
        adx[i].pb({p, 1}); adx[p].pb({i, 1});
      }
      int K = rnd() % n;
      cpx.K = K;
      cpx.solve(1);
      CHECK(cpx.global.valid_pairs == brute_dist_pairs(adx, n, K));
    }
  }

  // CentroidTree: distancia minima al nodo rojo mas cercano (camino 1-2-3-4-5)
  auto path_dist = [](int a, int b) { return (ll)abs(a - b); };
  CentroidTree<ll> ct(5);
  for (int i = 1; i < 5; i++) ct.add_edge(i, i + 1);
  ct.build(1);
  CHECK(ct.query(3, path_dist) == LINF); // ningun rojo
  ct.update(3, path_dist);
  CHECK(ct.query(3, path_dist) == 0);
  CHECK(ct.query(1, path_dist) == 2);
  CHECK(ct.query(5, path_dist) == 2);
  ct.update(1, path_dist);
  CHECK(ct.query(5, path_dist) == 2); // min(dist(5,3), dist(5,1))
  CHECK(ct.query(2, path_dist) == 1); // min(dist(2,3), dist(2,1))

  // CentroidBinarySearch: juez simulado en camino 1-2-3-4-5 con nodo oculto 5
  {
    CentroidBinarySearch cbs(5);
    for (int i = 1; i < 5; i++) cbs.add_edge(i, i + 1);
    // secuencia de queries "? u" -> respuestas dist(u, 5): 3,2,4,4,5,5 -> 2,3,1,1,0,0
    istringstream in("2 3 1 1 0 0");
    ostringstream out;
    auto* old_cin = cin.rdbuf(in.rdbuf());
    auto* old_cout = cout.rdbuf(out.rdbuf());
    int hidden_found = cbs.solve(1);
    cin.rdbuf(old_cin);
    cout.rdbuf(old_cout);
    CHECK(hidden_found == 5);
  }

  // DSU on tree: colores distintos por subarbol
  DSUOnTree sack(5);
  sack.add_edge(1, 2); sack.add_edge(1, 3); sack.add_edge(2, 4); sack.add_edge(2, 5);
  sack.set_colors({0, 1, 2, 1, 2, 1});
  sack.solve(1);
  CHECK(sack.ans[1] == 2 && sack.ans[2] == 2 && sack.ans[4] == 1 && sack.ans[3] == 1);

  // VirtualTree (nuevo API con lambdas): camino 1-2-3-4, query {2,4}
  Tree<ll> tv(4);
  for (int i = 1; i < 4; i++) tv.add_edge(i, i + 1);
  tv.init(1);
  BinaryLifting<ll> blv(tv, 1);
  VirtualTree vt(4);
  auto vt_build = [&](VirtualTree& v, vi qs) {
    return v.build(qs,
      [&](int u, int v) { return blv.lca(u, v); },
      [&](int u, int v) { return blv.is_ancestor(u, v); },
      [&](int u, int v) { return blv.dist_weighted(u, v); },
      [&](int a, int b) { return tv.tin[a] < tv.tin[b]; });
  };
  {
    auto [root, nodes] = vt_build(vt, {2, 4});
    CHECK(root == 2);
    CHECK(nodes.size() == 2);
    CHECK(vt.virt_adj[2].size() == 1 && vt.virt_adj[2][0].se == 2);
  }
  // arbol 1-2, 1-3, 3-4; query {2,4} -> el lca 1 entra al arbol virtual
  {
    Tree<ll> vt_tree(4);
    vt_tree.add_edge(1, 2); vt_tree.add_edge(1, 3); vt_tree.add_edge(3, 4);
    vt_tree.init(1);
    BinaryLifting<ll> bl2(vt_tree, 1);
    VirtualTree vt2(4);
    auto [root, nodes] = vt2.build({2, 4},
      [&](int u, int v) { return bl2.lca(u, v); },
      [&](int u, int v) { return bl2.is_ancestor(u, v); },
      [&](int u, int v) { return bl2.dist_weighted(u, v); },
      [&](int a, int b) { return vt_tree.tin[a] < vt_tree.tin[b]; });
    CHECK(root == 1);
    CHECK(nodes == (vi{1, 2, 4}));
    CHECK(vt2.virt_adj[1].size() == 2); // 1-2 (w=1) y 1-4 (w=2)
    ll w2 = 0, w4 = 0;
    for (auto [v, w] : vt2.virt_adj[1]) { if (v == 2) w2 = w; if (v == 4) w4 = w; }
    CHECK(w2 == 1 && w4 == 2);
  }
  // adaptabilidad: mismo arbol virtual usando FastLCA (O(1) por query)
  {
    Tree<ll> vt_tree(4);
    vt_tree.add_edge(1, 2); vt_tree.add_edge(1, 3); vt_tree.add_edge(3, 4);
    vt_tree.init(1);
    FastLCA<ll> fl(vt_tree);
    VirtualTree vt3(4);
    auto [root, nodes] = vt3.build({2, 4},
      [&](int u, int v) { return fl.lca(u, v); },
      [&](int u, int v) { return fl.lca(u, v) == u; },
      [&](int u, int v) { return (ll)fl.dist_edges(u, v); },
      [&](int a, int b) { return vt_tree.tin[a] < vt_tree.tin[b]; });
    CHECK(root == 1);
    CHECK(nodes == (vi{1, 2, 4}));
    CHECK(vt3.virt_adj[1].size() == 2);
  }

  // TreeHashing: dos estrellas de 3 hojas son isomorfas
  Tree<ll> t1(4), t2(4);
  for (int i = 2; i <= 4; i++) t1.add_edge(1, i), t2.add_edge(1, i);
  CHECK(TreeHashing<ll>::are_isomorphic(t1, t2));
  Tree<ll> t3(3);
  t3.add_edge(1, 2); t3.add_edge(2, 3); // camino de 3
  CHECK(!TreeHashing<ll>::are_isomorphic(t1, t3));

  // RerootingDP: maxima distancia (default) en camino 1-2-3
  RerootingDP rdp(3);
  rdp.add_edge(1, 2); rdp.add_edge(2, 3);
  auto res = rdp.solve(1);
  CHECK(res[2].val == 1 && res[1].val == 2 && res[3].val == 2);

  // FastLCA
  Tree<ll> tf(5);
  for (int i = 1; i < 5; i++) tf.add_edge(i, i + 1);
  tf.init(1);
  FastLCA<ll> fl(tf);
  CHECK(fl.lca(4, 5) == 4 && fl.lca(2, 5) == 2);
  CHECK(fl.dist_edges(1, 5) == 4);

  if (failures) { cout << "test_tree: " << failures << " fallos\n"; return 1; }
  cout << "test_tree: OK\n";
  return 0;
}