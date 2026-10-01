#include "test_base.hpp"
#include "rq.hpp"
using namespace std;
int failures = 0;
#define CHECK(cond) do { if (!(cond)) { cerr << "  FAIL L" << __LINE__ << ": " #cond "\n"; failures++; } } while (0)

int main() {
  vll a = {5, 2, 8, 1, 9, 3, 7, 4};

  // Sparse table (min)
  SparseTable st(a);
  CHECK(st.query(1, 4) == 1);
  CHECK(st.query(0, 0) == 5);
  CHECK(st.query(0, 7) == 1);
  CHECK(st.query(5, 7) == 3);

  // Fenwick
  FenwickTree bit(a.size());
  for (int i = 0; i < a.size(); i++) bit.add(i + 1, a[i]);
  CHECK(bit.query(2, 5) == 20);  // 2+8+1+9
  CHECK(bit.query(1, 8) == 39);  // suma total
  CHECK(bit.sum(1) == 5);

  // Iterative segtree (suma)
  IterativeSegTree ist(a);
  CHECK(ist.query(0, 3) == 16);
  ist.update(0, 100);
  CHECK(ist.query(0, 3) == 111);

  // Recursive segtree (max) + find_first
  RecursiveSegTree rst(a);
  CHECK(rst.query(0, 7) == 9);
  CHECK(rst.query(1, 4) == 9);
  rst.update(2, 20);
  CHECK(rst.query(0, 7) == 20);
  CHECK(rst.find_first(9) == 2); // primer indice con valor >= 9
  CHECK(rst.find_first(25) == -1);

  // Lazy segtree (suma, add y set)
  LazySegTree lst(a);
  CHECK(lst.query(1, 3) == 11);
  lst.add_range(1, 3, 10);
  CHECK(lst.query(1, 3) == 41);
  lst.set_range(1, 2, 5);
  CHECK(lst.query(1, 3) == 21);
  CHECK(lst.get(1) == 5 && lst.get(3) == 11);
  lst.set_range(0, 7, 1);
  CHECK(lst.query(0, 7) == 8);

  // Dynamic segtree
  DynamicSegTree dst(1e9);
  dst.update(5, 10);
  dst.update(100000000, 7);
  dst.update(5, 3);
  CHECK(dst.query(1, 1e9) == 20);
  CHECK(dst.query(6, 99999999) == 0);
  CHECK(dst.query(6, 100000000) == 7);
  CHECK(dst.query(100000001, 1e9) == 0);

  // Persistent segtree
  PersistentSegTree pst(a);
  int v1 = pst.update(0, 2, 100); // suma[2] += 100
  int v2 = pst.update(0, 2, 50);
  CHECK(pst.query(0, 0, 7) == 39);   // version original intacta
  CHECK(pst.query(v1, 0, 7) == 139);
  CHECK(pst.query(v2, 0, 7) == 89);

  // Fenwick 2D
  FenwickTree2D bit2(3, 3);
  bit2.add(2, 2, 10);
  bit2.add(3, 1, 5);
  CHECK(bit2.query(1, 1, 3, 3) == 15);
  CHECK(bit2.query(1, 1, 2, 2) == 10);

  // Merge sort tree
  MergeSortTree mst(a);
  CHECK(mst.count_less_than(0, 7, 5) == 4); // {2,1,3,4}
  CHECK(mst.count_less_equal(0, 7, 5) == 5); // + el 5
  CHECK(mst.count_less_than(2, 4, 9) == 2); // {8,1}
  mst.update(3, 100); // a[3]: 1 -> 100
  CHECK(mst.count_less_than(0, 7, 5) == 3); // {2,3,4}
  CHECK(mst.count_less_equal(0, 7, 5) == 4); // {2,3,4,5}

  // Sqrt decomposition
  SqrtDecomposition sq(a);
  CHECK(sq.query(0, 7) == 39);
  CHECK(sq.query(2, 5) == 21); // 8+1+9+3
  sq.update(0, 100);
  CHECK(sq.query(0, 7) == 134);

  // Treap
  Treap tp;
  tp.insert(5); tp.insert(2); tp.insert(8); tp.insert(2); tp.insert(1);
  CHECK(tp.kth_element(0) == 1);
  CHECK(tp.kth_element(1) == 2);
  CHECK(tp.kth_element(4) == 8);
  CHECK(tp.count_less_than(5) == 3); // 1,2,2
  tp.erase(2);
  CHECK(tp.count_less_than(5) == 2);
  CHECK(tp.kth_element(4) == LINF);

  // Implicit treap
  ImplicitTreap it;
  for (int i = 0; i < a.size(); i++) it.insert(i, a[i]);
  CHECK(it.query_range(0, 2) == 15);
  it.reverse_range(0, 7);
  CHECK(it.query_range(0, 2) == 14); // {4,7,3}
  it.erase(0);
  CHECK(it.query_range(0, 0) == 7);

  // Mo clasico: contar distintos en [l, r]
  {
    vi arr = {1, 2, 1, 3, 2, 3, 1};
    vector<pii> queries = {{0, 6}, {1, 3}, {2, 5}, {0, 2}, {3, 6}, {2, 2}, {4, 6}};
    auto brute = [&](int l, int r) {
      set<int> st;
      for (int i = l; i <= r; i++) st.insert(arr[i]);
      return (ll)st.size();
    };
    Mo mo(7);
    for (auto& [l, r] : queries) mo.add_query(l, r);
    vi cnt(4, 0);
    int distinct = 0;
    auto add = [&](int idx) { if (!cnt[arr[idx]]++) distinct++; };
    auto rem = [&](int idx) { if (!--cnt[arr[idx]]) distinct--; };
    auto get = [&]() { return (ll)distinct; };
    auto ans = mo.solve(add, rem, get);
    for (int i = 0; i < queries.size(); i++)
      CHECK(ans[i] == brute(queries[i].fi, queries[i].se));
  }

  // Mo con updates: distintos en [l, r] sobre un arreglo que cambia
  {
    vi arr = {1, 2, 1, 3, 2, 3, 1};
    vi base = arr; // copia original para el brute force
    MoWithUpdates mo(7);
    // updates: (pos, old, new); las queries se etiquetan con el tiempo actual
    mo.add_update(1, 2, 5); // arr[1]: 2 -> 5
    mo.add_query(0, 4);
    mo.add_update(4, 2, 1); // arr[4]: 2 -> 1
    mo.add_update(0, 1, 3); // arr[0]: 1 -> 3
    mo.add_query(0, 6);
    mo.add_query(2, 5);
    // tiempos de las queries: 1, 3, 3
    struct U { int pos, oldv, newv; };
    vector<U> ups = {{1, 2, 5}, {4, 2, 1}, {0, 1, 3}};
    int qt[3] = {1, 3, 3};
    pii qr[3] = {{0, 4}, {0, 6}, {2, 5}};
    vi cnt(6, 0);
    int distinct = 0;
    auto add = [&](int idx) { if (!cnt[arr[idx]]++) distinct++; };
    auto rem = [&](int idx) { if (!--cnt[arr[idx]]) distinct--; };
    auto apply = [&](const MoWithUpdates::Update& u, int L, int R) {
      if (u.pos >= L && u.pos <= R) {
        if (!--cnt[arr[u.pos]]) distinct--;
      }
      arr[u.pos] = u.new_val;
      if (u.pos >= L && u.pos <= R) {
        if (!cnt[arr[u.pos]]++) distinct++;
      }
    };
    auto revert = [&](const MoWithUpdates::Update& u, int L, int R) {
      if (u.pos >= L && u.pos <= R) {
        if (!--cnt[arr[u.pos]]) distinct--;
      }
      arr[u.pos] = u.old_val;
      if (u.pos >= L && u.pos <= R) {
        if (!cnt[arr[u.pos]]++) distinct++;
      }
    };
    auto get = [&]() { return (ll)distinct; };
    auto ans = mo.solve(add, rem, apply, revert, get);
    for (int i = 0; i < 3; i++) {
      vi cur = base;
      for (int j = 0; j < qt[i]; j++) cur[ups[j].pos] = ups[j].newv;
      set<int> st;
      for (int j = qr[i].fi; j <= qr[i].se; j++) st.insert(cur[j]);
      CHECK(ans[i] == (ll)st.size());
    }
  }

  // Rollback Mo: distintos en [l, r] solo con add (sin remove)
  {
    vi arr = {1, 2, 1, 3, 2, 3, 1};
    vector<pii> queries = {{0, 6}, {1, 3}, {2, 5}, {0, 2}, {3, 6}, {2, 2}, {4, 6}, {0, 4}};
    auto brute = [&](int l, int r) {
      set<int> st;
      for (int i = l; i <= r; i++) st.insert(arr[i]);
      return (ll)st.size();
    };
    RollbackMo rmo(7, 3); // bloque pequeno a proposito
    for (auto& [l, r] : queries) rmo.add_query(l, r);
    vi cnt(4, 0);
    int distinct = 0;
    vector<int> changes;
    auto add = [&](int idx) { changes.pb(idx); if (!cnt[arr[idx]]++) distinct++; };
    auto current = [&]() { return (ll)distinct; };
    auto snapshot = [&]() { changes.clear(); };
    auto rollback = [&]() {
      for (int idx : changes) if (--cnt[arr[idx]] == 0) distinct--;
    };
    auto reset = [&]() { fill(all(cnt), 0); distinct = 0; changes.clear(); };
    auto ans = rmo.solve(add, current, snapshot, rollback, reset);
    for (int i = 0; i < queries.size(); i++)
      CHECK(ans[i] == brute(queries[i].fi, queries[i].se));
  }

  // Stress aleatorio: Mo clasico vs brute force
  {
    mt19937 rnd(777);
    for (int it = 0; it < 20; it++) {
      int n = 1 + rnd() % 30, V = 1 + rnd() % 8;
      vi arr(n);
      for (int& x : arr) x = rnd() % V;
      int Q = 1 + rnd() % 40;
      vector<pii> queries(Q);
      for (auto& [l, r] : queries) {
        l = rnd() % n, r = rnd() % n;
        if (l > r) swap(l, r);
      }
      auto brute = [&](int l, int r) {
        set<int> st;
        for (int i = l; i <= r; i++) st.insert(arr[i]);
        return (ll)st.size();
      };
      Mo mo(n);
      for (auto& [l, r] : queries) mo.add_query(l, r);
      vi cnt(V, 0);
      int distinct = 0;
      auto add = [&](int idx) { if (!cnt[arr[idx]]++) distinct++; };
      auto rem = [&](int idx) { if (!--cnt[arr[idx]]) distinct--; };
      auto get = [&]() { return (ll)distinct; };
      auto ans = mo.solve(add, rem, get);
      for (int i = 0; i < Q; i++) CHECK(ans[i] == brute(queries[i].fi, queries[i].se));
      // Rollback Mo en el mismo caso
      RollbackMo rmo(n);
      for (auto& [l, r] : queries) rmo.add_query(l, r);
      vi cnt2(V, 0);
      int d2 = 0;
      vector<int> ch;
      auto add2 = [&](int idx) { ch.pb(idx); if (!cnt2[arr[idx]]++) d2++; };
      auto cur2 = [&]() { return (ll)d2; };
      auto snap2 = [&]() { ch.clear(); };
      auto rb2 = [&]() { for (int idx : ch) if (--cnt2[arr[idx]] == 0) d2--; };
      auto rst2 = [&]() { fill(all(cnt2), 0); d2 = 0; ch.clear(); };
      auto ans2 = rmo.solve(add2, cur2, snap2, rb2, rst2);
      for (int i = 0; i < Q; i++) CHECK(ans2[i] == brute(queries[i].fi, queries[i].se));
    }
  }

  // Stress aleatorio: Mo con updates vs brute force
  {
    mt19937 rnd(778);
    for (int it = 0; it < 10; it++) {
      int n = 1 + rnd() % 20, V = 1 + rnd() % 7;
      vi arr(n), base(n);
      for (int i = 0; i < n; i++) arr[i] = base[i] = rnd() % V;
      MoWithUpdates mo(n);
      vector<tuple<int, int, int>> ups; // pos, old, new
      vector<pair<int, int>> qrs;       // l, r
      vector<int> qtimes;               // numero de updates al agregar la query
      int steps = 1 + rnd() % 15;
      for (int s = 0; s < steps; s++) {
        if (rnd() % 2 == 0 || ups.empty()) { // update
          int pos = rnd() % n, oldv = arr[pos], newv = rnd() % V;
          arr[pos] = newv;
          ups.pb({pos, oldv, newv});
          mo.add_update(pos, oldv, newv);
        } else { // query
          int l = rnd() % n, r = rnd() % n;
          if (l > r) swap(l, r);
          qrs.pb({l, r});
          qtimes.pb(ups.size());
          mo.add_query(l, r);
        }
      }
      auto state_at = [&](int t) {
        vi cur = base;
        for (int i = 0; i < t; i++) cur[get<0>(ups[i])] = get<2>(ups[i]);
        return cur;
      };
      arr = base; // apply/revert presumen el arreglo inicial
      vi cnt(V, 0);
      int distinct = 0;
      auto add = [&](int idx) { if (!cnt[arr[idx]]++) distinct++; };
      auto rem = [&](int idx) { if (!--cnt[arr[idx]]) distinct--; };
      auto apply = [&](const MoWithUpdates::Update& u, int L, int R) {
        if (u.pos >= L && u.pos <= R) { if (!--cnt[arr[u.pos]]) distinct--; }
        arr[u.pos] = u.new_val;
        if (u.pos >= L && u.pos <= R) { if (!cnt[arr[u.pos]]++) distinct++; }
      };
      auto revert = [&](const MoWithUpdates::Update& u, int L, int R) {
        if (u.pos >= L && u.pos <= R) { if (!--cnt[arr[u.pos]]) distinct--; }
        arr[u.pos] = u.old_val;
        if (u.pos >= L && u.pos <= R) { if (!cnt[arr[u.pos]]++) distinct++; }
      };
      auto get = [&]() { return (ll)distinct; };
      auto ans = mo.solve(add, rem, apply, revert, get);
      for (int i = 0; i < qrs.size(); i++) {
        vi cur = state_at(qtimes[i]);
        set<int> st;
        for (int j = qrs[i].fi; j <= qrs[i].se; j++) st.insert(cur[j]);
        CHECK(ans[i] == (ll)st.size());
      }
    }
  }

  if (failures) { cout << "test_rq: " << failures << " fallos\n"; return 1; }
  cout << "test_rq: OK\n";
  return 0;
}