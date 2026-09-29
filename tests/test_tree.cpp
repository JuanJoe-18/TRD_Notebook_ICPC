#include "test_base.hpp"
#include "graph.hpp"
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
  auto far = t.all_farthest_distances();
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

  // HLD
  Tree<ll> th(4);
  for (int i = 1; i < 4; i++) th.add_edge(i, i + 1);
  HLD<ll> hld(4);
  for (int i = 1; i < 4; i++) hld.add_edge(i, i + 1);
  hld.init(1);
  hld.update_node(2, 5);
  hld.update_node(4, 3);
  CHECK(hld.query_path(2, 4) == 8);
  CHECK(hld.query_path(1, 4) == 8);
  hld.update_path(1, 4, 1);
  CHECK(hld.query_path(1, 4) == 12);
  CHECK(hld.lca(1, 4) == 1);

  // CentroidDecomposition offline: camino 1-2-3-4-5, contar pares con dist <= 2
  CentroidDecomposition cd(5);
  for (int i = 1; i < 5; i++) cd.add_edge(i, i + 1);
  cd.K = 2;
  cd.init(true);
  CHECK(cd.total_paths == 7);

  // CentroidDecomposition online (requiere LCA)
  Tree<ll> tco(5);
  for (int i = 1; i < 5; i++) tco.add_edge(i, i + 1);
  tco.init(1);
  vector<vi> adj(6);
  for (int u = 1; u <= 5; u++) for (auto& e : tco.adj[u]) adj[u].pb(e.to);
  LCA lca(5);
  lca.build(1, adj);
  CentroidDecomposition cdon(5);
  for (int i = 1; i < 5; i++) cdon.add_edge(i, i + 1);
  cdon.init(false);
  CHECK(cdon.query(5, lca) == LINF); // ningun rojo
  cdon.update(3, lca);
  CHECK(cdon.query(1, lca) == 2);
  CHECK(cdon.query(5, lca) == 2);

  // DSU on tree: colores distintos por subarbol
  DSUOnTree sack(5);
  sack.add_edge(1, 2); sack.add_edge(1, 3); sack.add_edge(2, 4); sack.add_edge(2, 5);
  sack.set_colors({0, 1, 2, 1, 2, 1});
  sack.solve(1);
  CHECK(sack.ans[1] == 2 && sack.ans[2] == 2 && sack.ans[4] == 1 && sack.ans[3] == 1);

  // VirtualTree
  Tree<ll> tv(4);
  for (int i = 1; i < 4; i++) tv.add_edge(i, i + 1);
  tv.init(1);
  BinaryLifting<ll> blv(tv, 1);
  VirtualTree vt(tv, blv);
  auto [root, nodes] = vt.build({2, 4});
  CHECK(root == 2);
  CHECK(nodes.size() == 2);
  CHECK(vt.virt_adj[2].size() == 1 && vt.virt_adj[2][0].se == 2);

  // TreeHashing: dos estrellas de 3 hojas son isomorfas
  vector<vi> g1(5), g2(5);
  for (int i = 2; i <= 4; i++) g1[1].pb(i), g1[i].pb(1);
  for (int i = 2; i <= 4; i++) g2[1].pb(i), g2[i].pb(1);
  CHECK(TreeHashing::are_isomorphic(4, g1, 4, g2));
  vector<vi> g3(4);
  g3[1] = {2}; g3[2] = {1, 3}; g3[3] = {2}; // camino de 3
  CHECK(!TreeHashing::are_isomorphic(4, g1, 3, g3));

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