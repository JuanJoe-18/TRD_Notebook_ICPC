#include "test_base.hpp"
#include "graph.hpp"
using namespace std;
int failures = 0;
#define CHECK(cond) do { if (!(cond)) { cerr << "  FAIL L" << __LINE__ << ": " #cond "\n"; failures++; } } while (0)

int main() {
  // DSU
  DSU dsu(5);
  CHECK(dsu.find(0) == 0);
  dsu.unite(0, 1); dsu.unite(2, 3); dsu.unite(3, 4);
  CHECK(dsu.find(1) == dsu.find(0));
  CHECK(dsu.find(4) == dsu.find(2));
  CHECK(dsu.component_size(0) == 2 && dsu.component_size(4) == 3);
  CHECK(dsu.comps == 2);
  CHECK(!dsu.unite(0, 1)); // ya unidos

  // Dijkstra + camino
  Graph<ll> g(4);
  g.add_directed_edge(1, 2, 10, 0);
  g.add_directed_edge(2, 3, 5, 1);
  g.add_directed_edge(1, 4, 2, 2);
  g.add_directed_edge(4, 3, 20, 3);
  Dijkstra<ll> dih;
  auto [dist, par] = dih.run(g, 1);
  CHECK(dist[3] == 15);
  CHECK(dih.restore_path(3) == (vi{1, 2, 3}));

  // invalid_edges: bloquear la arista id 0 (1->2) -> el camino pasa por 1->4->3
  auto [dist2, par2] = dih.run(g, 1, vector<bool>{true, false, false, false});
  CHECK(dist2[3] == 22);
  CHECK(dist2[2] == LINF);
  CHECK(dih.restore_path(3) == (vi{1, 4, 3}));

  // nodo inalcanzable -> restore_path vacio
  Graph<ll> g2(3);
  g2.add_directed_edge(1, 2, 1, 0); // nodo 3 aislado
  Dijkstra<ll> dih2;
  dih2.run(g2, 1);
  CHECK(dih2.restore_path(3).empty());
  CHECK(dih2.restore_path(2) == (vi{1, 2}));

  // BFS 0-1
  Graph<ll> g01(3);
  g01.add_directed_edge(1, 2, 0);
  g01.add_directed_edge(2, 3, 1);
  g01.add_directed_edge(1, 3, 1);
  auto d01 = ZeroOneBfs<ll>().run(g01, 1);
  CHECK(d01[3] == 1 && d01[2] == 0);

  // bipartito
  Graph<ll> bip(4);
  bip.add_undirected_edge(1, 2); bip.add_undirected_edge(2, 3);
  bip.add_undirected_edge(3, 4); bip.add_undirected_edge(4, 1);
  auto [ok, col] = Bipartite<ll>().check(bip);
  CHECK(ok);
  CHECK(col[1] == col[3] && col[1] != col[2]);
  Graph<ll> tri(3);
  tri.add_undirected_edge(1, 2); tri.add_undirected_edge(2, 3); tri.add_undirected_edge(3, 1);
  CHECK(!Bipartite<ll>().check(tri).fi);

  // Floyd-Warshall
  auto allp = FloydWarshall<ll>().run(g);
  CHECK(allp[1][3] == 15 && allp[1][1] == 0);

  // topo
  Graph<ll> dag(4);
  dag.add_directed_edge(1, 2); dag.add_directed_edge(2, 3); dag.add_directed_edge(1, 4);
  auto topo = TopoSort<ll>().run(dag);
  CHECK(topo.size() == 4);
  CHECK(find(all(topo), 1) < find(all(topo), 2) && find(all(topo), 2) < find(all(topo), 3));
  Graph<ll> cyc(2);
  cyc.add_directed_edge(1, 2); cyc.add_directed_edge(2, 1);
  CHECK(TopoSort<ll>().run(cyc).empty());

  // SCC
  Graph<ll> scc(3);
  scc.add_directed_edge(1, 2); scc.add_directed_edge(2, 1); scc.add_directed_edge(3, 2);
  auto comp = TarjanSCC<ll>().run(scc);
  CHECK(comp[1] == comp[2]);
  CHECK(comp[3] != comp[1]);

  // Bellman-Ford: sin ciclo negativo
  Graph<ll> bf(3);
  bf.add_directed_edge(1, 2, 4); bf.add_directed_edge(2, 3, -2); bf.add_directed_edge(1, 3, 5);
  auto [okBf, distBf, pBf, cycleBf] = BellmanFord<ll>().run(bf, 1);
  CHECK(okBf && distBf[3] == 2);
  // con ciclo negativo
  Graph<ll> neg(2);
  neg.add_directed_edge(1, 2, -1); neg.add_directed_edge(2, 1, -1);
  auto [okNeg, dNeg, pNeg, cyNeg] = BellmanFord<ll>().run(neg, 1);
  CHECK(!okNeg && !cyNeg.empty());

  // ciclos
  Graph<ll> dc(3);
  dc.add_directed_edge(1, 2); dc.add_directed_edge(2, 3); dc.add_directed_edge(3, 1);
  auto dcycle = DirectedCycle<ll>().find(dc);
  CHECK(dcycle.size() == 4 && dcycle[0] == dcycle.back());
  Graph<ll> uc(3);
  uc.add_undirected_edge(1, 2); uc.add_undirected_edge(2, 3); uc.add_undirected_edge(3, 1);
  CHECK(!UndirectedCycle<ll>().find(uc).empty());

  // Kruskal MST
  Graph<ll> mst(3);
  mst.add_undirected_edge(1, 2, 1);
  mst.add_undirected_edge(2, 3, 2);
  mst.add_undirected_edge(1, 3, 3);
  auto [wt, edges_mst] = Kruskal<ll>().run(mst);
  CHECK(wt == 3 && edges_mst.size() == 2);

  // Euleriano no dirigido: triangulo
  Graph<ll> eul(3);
  eul.add_undirected_edge(1, 2, 1, 0);
  eul.add_undirected_edge(2, 3, 1, 1);
  eul.add_undirected_edge(3, 1, 1, 2);
  auto ec = EulerianCircuit<ll>().find(eul, 1, true);
  CHECK(ec.size() == 4 && ec[0] == ec.back());

  // FunctionalGraph: 1->2, 2->3, 3->2
  FunctionalGraph fg(3);
  fg.build_lifting({0, 2, 3, 2});
  CHECK(fg.kth_successor(1, 1) == 2 && fg.kth_successor(1, 2) == 3 && fg.kth_successor(1, 3) == 2);
  fg.decompose({0, 2, 3, 2});
  CHECK(fg.in_cycle[2] && fg.in_cycle[3] && !fg.in_cycle[1]);
  CHECK(fg.cycle_size[2] == 2 && fg.dist_to_cycle[1] == 1);

  // LCA en arbol 1-2,1-3,3-4
  vector<vi> tree_adj(5);
  tree_adj[1] = {2, 3}; tree_adj[2] = {1}; tree_adj[3] = {1, 4}; tree_adj[4] = {3};
  LCA lca(4);
  lca.build(1, tree_adj);
  CHECK(lca.lca(2, 4) == 1 && lca.lca(3, 4) == 3);
  CHECK(lca.dist(2, 4) == 3);

  // 2-SAT
  TwoSat sat(2);
  sat.add_clause(1, true, 2, true);
  sat.add_clause(1, false, 2, true);
  CHECK(sat.solve());
  TwoSat unsat(1);
  unsat.force_value(1, true);
  unsat.force_value(1, false);
  CHECK(!unsat.solve());

  // Kuhn
  BipartiteMatcher bm(2, 2);
  bm.add_edge(1, 1); bm.add_edge(1, 2); bm.add_edge(2, 2);
  CHECK(bm.max_matching() == 2);
  // Hopcroft-Karp
  HopcroftKarp hk(2, 2);
  hk.add_edge(1, 1); hk.add_edge(1, 2); hk.add_edge(2, 2);
  CHECK(hk.max_matching() == 2);

  // BridgeTree: triangulo + puente 3-4
  BridgeTree bt(4, 4);
  bt.add_edge(1, 2, 0); bt.add_edge(2, 3, 1); bt.add_edge(3, 1, 2);
  bt.add_edge(3, 4, 3);
  bt.build();
  CHECK(bt.is_bridge[3] && !bt.is_bridge[0]);

  // Dinic: max flow
  Dinic din(4, 1, 4);
  din.add_edge(1, 2, 3); din.add_edge(1, 3, 2);
  din.add_edge(2, 3, 1); din.add_edge(2, 4, 2); din.add_edge(3, 4, 3);
  CHECK(din.max_flow() == 5);

  // MCMF
  MCMF mcmf(3);
  mcmf.add_edge(1, 2, 2, 3);
  mcmf.add_edge(2, 3, 2, 1);
  mcmf.add_edge(1, 3, 2, 5);
  auto [fl, cs] = mcmf.solve(1, 3);
  CHECK(fl == 4 && cs == 18);

  // DominatorTree: 1->2, 1->3, 2->4, 3->4
  DominatorTree dom(4);
  dom.add_edge(1, 2); dom.add_edge(1, 3); dom.add_edge(2, 4); dom.add_edge(3, 4);
  dom.build(1);
  CHECK(dom.dom[1].size() == 3); // 2, 3, 4

  // MaxClique (0-based): K4 -> 4; cuadrado con diagonal -> 3; triangulo + aislado -> 3
  {
    MaxClique k4(4);
    for (int i = 0; i < 4; i++)
      for (int j = i + 1; j < 4; j++) k4.add_edge(i, j);
    CHECK(k4.solve() == 4);
  }
  {
    MaxClique sq(4);
    sq.add_edge(0, 1); sq.add_edge(1, 2); sq.add_edge(2, 3); sq.add_edge(3, 0);
    sq.add_edge(0, 2); // diagonal
    CHECK(sq.solve() == 3);
  }
  {
    MaxClique tri(4);
    tri.add_edge(0, 1); tri.add_edge(1, 2); tri.add_edge(2, 0); // nodo 3 aislado
    CHECK(tri.solve() == 3);
  }

  // DinicLB: flujo con lower bounds
  {
    // 1->2 [1,2], 2->3 [1,2], 1->3 [0,10]; max flow 1->3 = 12
    DinicLB dlb(3);
    int e0 = dlb.add_edge(1, 2, 1, 2);
    int e1 = dlb.add_edge(2, 3, 1, 2);
    int e2 = dlb.add_edge(1, 3, 0, 10);
    CHECK(dlb.max_flow(1, 3) == 12);
    CHECK(dlb.edge_flow(e0) == 2);
    CHECK(dlb.edge_flow(e1) == 2);
    CHECK(dlb.edge_flow(e2) == 10);
  }
  {
    // 1->2 [3,5]; max flow 1->2 = 5
    DinicLB dlb(2);
    int e0 = dlb.add_edge(1, 2, 3, 5);
    CHECK(dlb.max_flow(1, 2) == 5);
    CHECK(dlb.edge_flow(e0) == 5);
  }
  {
    // Ciclo con demandas inconsistentes -> infeasible
    // 1->2 [3,3], 2->3 [1,1], 3->1 [5,5]
    DinicLB dlb(3);
    dlb.add_edge(1, 2, 3, 3);
    dlb.add_edge(2, 3, 1, 1);
    dlb.add_edge(3, 1, 5, 5);
    CHECK(dlb.max_flow(1, 3) == -1);
  }
  {
    // circulacion factible sin s/t: 1->2 [2,4], 2->3 [1,3], 3->1 [1,2]
    DinicLB dlb(3);
    dlb.add_edge(1, 2, 2, 4);
    dlb.add_edge(2, 3, 1, 3);
    dlb.add_edge(3, 1, 1, 2);
    CHECK(dlb.feasible());
    // circulacion infactible: 1->2 [5,5], 2->1 [3,3] (desbalance neto)
    DinicLB bad(2);
    bad.add_edge(1, 2, 5, 5);
    bad.add_edge(2, 1, 3, 3);
    CHECK(!bad.feasible());
  }

  if (failures) { cout << "test_graph: " << failures << " fallos\n"; return 1; }
  cout << "test_graph: OK\n";
  return 0;
}