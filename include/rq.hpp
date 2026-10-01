/**
 * Uso[0]: SparseTable st(a); st.query(L, R);
 * Consultas estaticas en rangos para operaciones idempotentes (min, max, gcd).
 * Complejidad: construccion $O(N \log N)$, consulta $O(1)$.
 */
struct SparseTable {
  int n;
  vector<vll> st;
  SparseTable(const vll& a) : n(a.size()) {
    int log_n = log2(n) + 1;
    st.assign(n, vll(log_n));
    for (int i = 0; i < n; i++) st[i][0] = a[i];
    for (int j = 1; j < log_n; j++)
      for (int i = 0; i + (1 << j) <= n; i++)
        st[i][j] = min(st[i][j - 1], st[i + (1 << (j - 1))][j - 1]); // min/max/gcd
  }
  ll query(int L, int R) {
    int j = log2(R - L + 1);
    return min(st[L][j], st[R - (1 << j) + 1][j]);
  }
};

/**
 * Uso[1]: FenwickTree ft(n); ft.add(i, x); ll s = ft.query(2, 3);
 * Arbol de Fenwick para sumas de prefijos y rangos con actualizaciones puntuales.
 * Complejidad: actualizacion y consulta $O(\log N)$.
 */
struct FenwickTree {
  int n;
  vll bit;
  FenwickTree(int n) : n(n), bit(n + 1, 0) {}
  void add(int idx, ll val) { for (; idx <= n; idx += idx & -idx) bit[idx] += val; }
  ll sum(int idx) {
    ll res = 0;
    for (; idx > 0; idx -= idx & -idx) res += bit[idx];
    return res;
  }
  ll query(int L, int R) { return sum(R) - sum(L - 1); }
};

/**
 * Uso[0]: IterativeSegTree st(a); st.update(pos, val); st.query(L, R);
 * Arbol de segmentos iterativo para suma (o min/max) con actualizacion puntual.
 * Complejidad: actualizacion y consulta $O(\log N)$.
 */
struct IterativeSegTree {
  int n;
  vll st;
  IterativeSegTree(const vll& a) : n(a.size()), st(2 * a.size(), 0) {
    for (int i = 0; i < n; i++) st[n + i] = a[i];
    for (int i = n - 1; i > 0; i--) st[i] = st[i << 1] + st[i << 1 | 1]; // + / min / max
  }
  void update(int pos, ll val) {
    for (st[pos += n] = val; pos > 1; pos >>= 1) st[pos >> 1] = st[pos] + st[pos ^ 1];
  }
  ll query(int L, int R) {
    ll resL = 0, resR = 0;
    for (L += n, R += n + 1; L < R; L >>= 1, R >>= 1) {
      if (L & 1) resL += st[L++];
      if (R & 1) resR += st[--R];
    }
    return resL + resR;
  }
};

/**
 * Uso[0]: RecursiveSegTree st(a); st.update(pos, val); st.query(L, R); st.find_first(val);
 * Arbol de segmentos recursivo con estructura de nodo y metodo merge para hacer variaciones faciles (max, +, min).
 * Para cambiar la operacion, edita merge (mezcla de dos nodos) y neutral (identidad de la operacion).
 * Indexación: 0-indexado.
 * Complejidad: actualizacion y consulta $O(\log N)$.
 */
struct RecursiveSegTree {
  struct Node { ll val; };
  int n;
  vector<Node> st;
  Node merge(const Node& a, const Node& b) { return {max(a.val, b.val)}; } // max / + / min
  Node neutral() { return {-LINF}; } // identidad de merge: -LINF (max), 0 (suma), LINF (min)
  RecursiveSegTree(const vll& a) : n(a.size()), st(4 * a.size()) { build(1, 0, n - 1, a); }
  void build(int p, int L, int R, const vll& a) {
    if (L == R) { st[p] = {a[L]}; return; }
    int mid = (L + R) / 2;
    build(p << 1, L, mid, a);
    build(p << 1 | 1, mid + 1, R, a);
    st[p] = merge(st[p << 1], st[p << 1 | 1]);
  }
  void update(int p, int L, int R, int pos, ll val) {
    if (L == R) { st[p] = {val}; return; } // o st[p].val += val para sumar
    int mid = (L + R) / 2;
    pos <= mid ? update(p << 1, L, mid, pos, val) : update(p << 1 | 1, mid + 1, R, pos, val);
    st[p] = merge(st[p << 1], st[p << 1 | 1]);
  }
  Node query(int p, int L, int R, int qL, int qR) {
    if (qL > R || qR < L) return neutral();
    if (qL <= L && R <= qR) return st[p];
    int mid = (L + R) / 2;
    return merge(query(p << 1, L, mid, qL, qR), query(p << 1 | 1, mid + 1, R, qL, qR));
  }
  // Primer indice con valor >= val (asume arbol de max). -1 si no existe
  int find_first(int p, int L, int R, ll val) {
    if (st[p].val < val) return -1;
    if (L == R) return L;
    int mid = (L + R) / 2;
    return st[p << 1].val >= val ? find_first(p << 1, L, mid, val) : find_first(p << 1 | 1, mid + 1, R, val);
  }
  void update(int pos, ll val) { update(1, 0, n - 1, pos, val); }
  ll query(int L, int R) { return query(1, 0, n - 1, L, R).val; }
  int find_first(ll val) { return find_first(1, 0, n - 1, val); }
};

/**
 * Uso[0]: LazySegTree st(a); st.add_range(L, R, v); st.set_range(L, R, v); st.query(L, R); st.get(pos);
 * Arbol de segmentos con propagacion perezosa y estructura de nodo + metodo merge para suma (o min/max) con add y set por rango.
 * Para cambiar la operacion, edita merge y la manera en que apply_add/apply_set combinan el valor de cada nodo.
 * Indexación: 0-indexado.
 * Complejidad: actualizacion y consulta de rango $O(\log N)$.
 */
struct LazySegTree {
  struct Node { ll sum = 0; };
  int n;
  vector<Node> st;
  vll lazy_add, lazy_set;
  vector<bool> marked;
  Node merge(const Node& a, const Node& b) { return {a.sum + b.sum}; } // + / min / max
  LazySegTree(const vll& a) : n(a.size()), st(4 * a.size()), lazy_add(4 * a.size()),
    lazy_set(4 * a.size()), marked(4 * a.size(), false) { build(1, 0, n - 1, a); }
  void apply_set(int p, ll val, int L, int R) {
    st[p].sum = val * (R - L + 1); // para min/max: st[p].sum = val
    lazy_set[p] = val;
    lazy_add[p] = 0;
    marked[p] = true;
  }
  void apply_add(int p, ll val, int L, int R) {
    st[p].sum += val * (R - L + 1); // para min/max: st[p].sum += val
    if (marked[p]) lazy_set[p] += val;
    else lazy_add[p] += val;
  }
  void push(int p, int L, int R) {
    int mid = (L + R) / 2;
    if (marked[p]) {
      apply_set(p << 1, lazy_set[p], L, mid);
      apply_set(p << 1 | 1, lazy_set[p], mid + 1, R);
      marked[p] = false;
    }
    if (lazy_add[p]) {
      apply_add(p << 1, lazy_add[p], L, mid);
      apply_add(p << 1 | 1, lazy_add[p], mid + 1, R);
      lazy_add[p] = 0;
    }
  }
  void build(int p, int L, int R, const vll& a) {
    if (L == R) { st[p] = {a[L]}; return; }
    int mid = (L + R) / 2;
    build(p << 1, L, mid, a);
    build(p << 1 | 1, mid + 1, R, a);
    st[p] = merge(st[p << 1], st[p << 1 | 1]);
  }
  void update(int p, int L, int R, int qL, int qR, ll val, int type) { // 1 = add, 2 = set
    if (qL > R || qR < L) return;
    if (qL <= L && R <= qR) { type == 1 ? apply_add(p, val, L, R) : apply_set(p, val, L, R); return; }
    push(p, L, R);
    int mid = (L + R) / 2;
    update(p << 1, L, mid, qL, qR, val, type);
    update(p << 1 | 1, mid + 1, R, qL, qR, val, type);
    st[p] = merge(st[p << 1], st[p << 1 | 1]);
  }
  Node query(int p, int L, int R, int qL, int qR) {
    if (qL > R || qR < L) return {0};
    if (qL <= L && R <= qR) return st[p];
    push(p, L, R);
    int mid = (L + R) / 2;
    return merge(query(p << 1, L, mid, qL, qR), query(p << 1 | 1, mid + 1, R, qL, qR));
  }
  Node get_point(int p, int L, int R, int pos) {
    if (L == R) return st[p];
    push(p, L, R);
    int mid = (L + R) / 2;
    return pos <= mid ? get_point(p << 1, L, mid, pos) : get_point(p << 1 | 1, mid + 1, R, pos);
  }
  void add_range(int L, int R, ll val) { update(1, 0, n - 1, L, R, val, 1); }
  void set_range(int L, int R, ll val) { update(1, 0, n - 1, L, R, val, 2); }
  ll query(int L, int R) { return query(1, 0, n - 1, L, R).sum; }
  ll get(int pos) { return get_point(1, 0, n - 1, pos).sum; }
};

/**
 * Uso[1]: DynamicSegTree st(max_rango); st.update(pos, val); st.query(L, R);
 * Arbol de segmentos dinamico (sparse) para indices grandes sin compresion previa.
 * Complejidad: actualizacion y consulta $O(\log \text{max\_rango})$.
 */
struct DynamicSegTree {
  struct Node { ll sum = 0; int left = 0, right = 0; };
  vector<Node> st;
  ll max_range;
  int root;
  DynamicSegTree(ll max_range = 1e9) : max_range(max_range), root(0) { st.emplace_back(); }
  int new_node() { st.emplace_back(); return st.size() - 1; }
  int update(int node, ll L, ll R, ll pos, ll val) {
    if (!node) node = new_node();
    st[node].sum += val;
    if (L == R) return node;
    ll mid = (L + R) / 2;
    if (pos <= mid) st[node].left = update(st[node].left, L, mid, pos, val);
    else st[node].right = update(st[node].right, mid + 1, R, pos, val);
    return node;
  }
  ll query(int node, ll L, ll R, ll qL, ll qR) {
    if (!node || qL > R || qR < L) return 0;
    if (qL <= L && R <= qR) return st[node].sum;
    ll mid = (L + R) / 2;
    return query(st[node].left, L, mid, qL, qR) + query(st[node].right, mid + 1, R, qL, qR);
  }
  void update(ll pos, ll val) { root = update(root, 1, max_range, pos, val); }
  ll query(ll L, ll R) { return query(root, 1, max_range, L, R); }
};

/**
 * Uso[0]: PersistentSegTree pst(a); int v = pst.update(version, pos, val); pst.query(v, L, R);
 * Arbol de segmentos persistente con versiones historicas y consultas sobre cualquier version.
 * Complejidad: actualizacion y consulta $O(\log N)$; espacio $O(N \log N)$.
 */
struct PersistentSegTree {
  struct Node { ll sum; int left, right; Node(ll s = 0, int l = 0, int r = 0) : sum(s), left(l), right(r) {} };
  vector<Node> st;
  vi roots;
  int n;
  int new_node(ll sum, int left, int right) { st.emplace_back(sum, left, right); return st.size() - 1; }
  int build(int L, int R, const vll& a) {
    if (L == R) return new_node(a[L], 0, 0);
    int mid = (L + R) / 2;
    int l = build(L, mid, a), r = build(mid + 1, R, a);
    return new_node(st[l].sum + st[r].sum, l, r);
  }
  int update(int old, int L, int R, int pos, ll val) {
    if (L == R) return new_node(st[old].sum + val, 0, 0);
    int mid = (L + R) / 2;
    int l = st[old].left, r = st[old].right;
    if (pos <= mid) l = update(l, L, mid, pos, val);
    else r = update(r, mid + 1, R, pos, val);
    return new_node(st[l].sum + st[r].sum, l, r);
  }
  ll query(int node, int L, int R, int qL, int qR) {
    if (!node || qL > R || qR < L) return 0;
    if (qL <= L && R <= qR) return st[node].sum;
    int mid = (L + R) / 2;
    return query(st[node].left, L, mid, qL, qR) + query(st[node].right, mid + 1, R, qL, qR);
  }
  PersistentSegTree(const vll& a) : n(a.size()) {
    st.emplace_back();
    roots.pb(build(0, n - 1, a));
  }
  int update(int prev_version, int pos, ll val) { // crea nueva version
    roots.pb(update(roots[prev_version], 0, n - 1, pos, val));
    return roots.size() - 1;
  }
  ll query(int version, int L, int R) { return query(roots[version], 0, n - 1, L, R); }
};

/**
 * Uso[1]: FenwickTree2D ft(n, m); ft.add(r, c, v); ft.query(r1, c1, r2, c2);
 * Arbol de Fenwick bidimensional para sumas de submatrices con actualizaciones puntuales.
 * Complejidad: actualizacion y consulta $O(\log n \cdot \log m)$.
 */
struct FenwickTree2D {
  int n, m;
  vector<vll> bit;
  FenwickTree2D(int n, int m) : n(n), m(m), bit(n + 1, vll(m + 1, 0)) {}
  void add(int r, int c, ll val) {
    for (int i = r; i <= n; i += i & -i)
      for (int j = c; j <= m; j += j & -j) bit[i][j] += val;
  }
  ll sum(int r, int c) {
    ll res = 0;
    for (int i = r; i > 0; i -= i & -i)
      for (int j = c; j > 0; j -= j & -j) res += bit[i][j];
    return res;
  }
  ll query(int r1, int c1, int r2, int c2) {
    return sum(r2, c2) - sum(r1 - 1, c2) - sum(r2, c1 - 1) + sum(r1 - 1, c1 - 1);
  }
};

/**
 * Uso[0]: MergeSortTree st(a); st.count_less_than(L, R, val); st.count_less_equal(L, R, val); st.update(pos, val);
 * Arbol de segmentos con listas ordenadas para contar elementos menores o iguales a un valor en un rango; soporta actualizaciones puntuales.
 * Indexación: 0-indexado.
 * Complejidad: construccion $O(N \log N)$, consulta $O(\log^2 N)$, actualizacion $O(N)$ (reconstruye el camino por mezcla).
 */
struct MergeSortTree {
  int n;
  vector<vll> st;
  MergeSortTree(const vll& a) : n(a.size()), st(4 * a.size()) { build(1, 0, n - 1, a); }
  void build(int p, int L, int R, const vll& a) {
    if (L == R) { st[p] = {a[L]}; return; }
    int mid = (L + R) / 2;
    build(p << 1, L, mid, a);
    build(p << 1 | 1, mid + 1, R, a);
    merge_children(p);
  }
  void merge_children(int p) {
    st[p].resize(st[p << 1].size() + st[p << 1 | 1].size());
    merge(all(st[p << 1]), all(st[p << 1 | 1]), st[p].begin());
  }
  void update(int p, int L, int R, int pos, ll val) {
    if (L == R) { st[p] = {val}; return; }
    int mid = (L + R) / 2;
    pos <= mid ? update(p << 1, L, mid, pos, val) : update(p << 1 | 1, mid + 1, R, pos, val);
    merge_children(p);
  }
  int less_than(int p, int L, int R, int qL, int qR, ll val) {
    if (qL > R || qR < L) return 0;
    if (qL <= L && R <= qR) return lower_bound(all(st[p]), val) - st[p].begin();
    int mid = (L + R) / 2;
    return less_than(p << 1, L, mid, qL, qR, val) + less_than(p << 1 | 1, mid + 1, R, qL, qR, val);
  }
  int less_equal(int p, int L, int R, int qL, int qR, ll val) {
    if (qL > R || qR < L) return 0;
    if (qL <= L && R <= qR) return upper_bound(all(st[p]), val) - st[p].begin();
    int mid = (L + R) / 2;
    return less_equal(p << 1, L, mid, qL, qR, val) + less_equal(p << 1 | 1, mid + 1, R, qL, qR, val);
  }
  void update(int pos, ll val) { update(1, 0, n - 1, pos, val); }
  int count_less_than(int L, int R, ll val) { return less_than(1, 0, n - 1, L, R, val); }
  int count_less_equal(int L, int R, ll val) { return less_equal(1, 0, n - 1, L, R, val); }
};

/**
 * Uso[0]: SqrtDecomposition sq(a); sq.update(pos, val); sq.query(L, R);
 * Descomposicion en bloques para sumas de rangos con actualizaciones puntuales.
 * Complejidad: actualizacion y consulta $O(\sqrt{N})$.
 */
struct SqrtDecomposition {
  int n, block_size;
  vll a, blocks;
  SqrtDecomposition(const vll& arr) : n(arr.size()), block_size(max(1, (int)sqrt(n))), a(arr) {
    blocks.assign((n + block_size - 1) / block_size, 0);
    for (int i = 0; i < n; i++) blocks[i / block_size] += a[i];
  }
  void update(int pos, ll val) {
    int b = pos / block_size;
    a[pos] = val;
    blocks[b] = 0;
    int end = min(n, (b + 1) * block_size);
    for (int i = b * block_size; i < end; i++) blocks[b] += a[i];
  }
  ll query(int L, int R) {
    int bL = L / block_size, bR = R / block_size;
    ll res = 0;
    if (bL == bR) {
      for (int i = L; i <= R; i++) res += a[i];
      return res;
    }
    for (int i = L; i < (bL + 1) * block_size; i++) res += a[i];
    for (int b = bL + 1; b < bR; b++) res += blocks[b];
    for (int i = bR * block_size; i <= R; i++) res += a[i];
    return res;
  }
};

/**
 * Uso[0]: Treap t; t.insert(x); t.erase(x); t.kth_element(k); t.count_less_than(x);
 * Arbol binario de busqueda balanceado aleatoriamente (treap) con consultas por orden.
 * Complejidad: cada operacion $O(\log N)$ esperado.
 */

mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
struct TreapNode {
  ll key, prior;
  int sz;
  TreapNode *l, *r;
  TreapNode(ll key) : key(key), prior(rng()), sz(1), l(nullptr), r(nullptr) {}
};
typedef TreapNode* pNode;
struct Treap {
  pNode root = nullptr;
  int sz(pNode t) { return t ? t->sz : 0; }
  void upd(pNode t) { if (t) t->sz = 1 + sz(t->l) + sz(t->r); }
  void split(pNode t, ll key, pNode& l, pNode& r) { // l: keys <= key
    if (!t) { l = r = nullptr; return; }
    if (key >= t->key) split(t->r, key, t->r, r), l = t;
    else split(t->l, key, l, t->l), r = t;
    upd(t);
  }
  void merge(pNode& t, pNode l, pNode r) {
    if (!l || !r) { t = l ? l : r; return; }
    if (l->prior > r->prior) merge(l->r, l->r, r), t = l;
    else merge(r->l, l, r->l), t = r;
    upd(t);
  }
  void insert(ll key) {
    pNode l, r;
    split(root, key, l, r);
    pNode node = new TreapNode(key);
    merge(l, l, node);
    merge(root, l, r);
  }
  void erase(ll key) {
    pNode l, m, r;
    split(root, key - 1, l, r);
    split(r, key, m, r);
    if (m) {
      pNode tmp = m;
      merge(m, m->l, m->r);
      delete tmp;
    }
    merge(root, l, m);
    merge(root, root, r);
  }
  ll kth_element(int k) { // k-esimo menor (0-indexed). LINF si fuera de rango
    pNode cur = root;
    while (cur) {
      int ls = sz(cur->l);
      if (ls == k) return cur->key;
      if (ls > k) cur = cur->l;
      else k -= ls + 1, cur = cur->r;
    }
    return LINF;
  }
  int count_less_than(ll key) {
    pNode l, r;
    split(root, key - 1, l, r);
    int res = sz(l);
    merge(root, l, r);
    return res;
  }
};

/**
 * Uso[0]: ImplicitTreap it; it.insert(pos, val); it.erase(pos); it.query_range(L, R); it.reverse_range(L, R);
 * Treap implicito para secuencias dinamicas con insercion/borrado por posicion, suma de rangos y reversiones perezosas.
 * Complejidad: cada operacion $O(\log N)$ esperado.
 */
struct ImplicitNode {
  ll val, sum, prior;
  int sz;
  bool lazy_rev;
  ImplicitNode *l, *r;
  ImplicitNode(ll val) : val(val), sum(val), prior(rng()), sz(1), lazy_rev(false), l(nullptr), r(nullptr) {}
};
typedef ImplicitNode* pINode;
struct ImplicitTreap {
  pINode root = nullptr;
  int sz(pINode t) { return t ? t->sz : 0; }
  ll sum(pINode t) { return t ? t->sum : 0; }
  void push(pINode t) {
    if (t && t->lazy_rev) {
      t->lazy_rev = false;
      swap(t->l, t->r);
      if (t->l) t->l->lazy_rev ^= true;
      if (t->r) t->r->lazy_rev ^= true;
    }
  }
  void upd(pINode t) {
    if (t) t->sz = 1 + sz(t->l) + sz(t->r), t->sum = t->val + sum(t->l) + sum(t->r);
  }
  void split(pINode t, int k, pINode& l, pINode& r) { // l: primeros k elementos
    if (!t) { l = r = nullptr; return; }
    push(t);
    int implicit = sz(t->l) + 1;
    if (k >= implicit) split(t->r, k - implicit, t->r, r), l = t;
    else split(t->l, k, l, t->l), r = t;
    upd(t);
  }
  void merge(pINode& t, pINode l, pINode r) {
    push(l); push(r);
    if (!l || !r) { t = l ? l : r; return; }
    if (l->prior > r->prior) merge(l->r, l->r, r), t = l;
    else merge(r->l, l, r->l), t = r;
    upd(t);
  }
  void insert(int pos, ll val) {
    pINode l, r;
    split(root, pos, l, r);
    pINode node = new ImplicitNode(val);
    merge(l, l, node);
    merge(root, l, r);
  }
  void erase(int pos) {
    pINode l, mid, r;
    split(root, pos, l, r);
    split(r, 1, mid, r);
    delete mid;
    merge(root, l, r);
  }
  void reverse_range(int L, int R) {
    pINode l, mid, r;
    split(root, L, l, mid);
    split(mid, R - L + 1, mid, r);
    mid->lazy_rev ^= true;
    merge(root, l, mid);
    merge(root, root, r);
  }
  ll query_range(int L, int R) {
    pINode l, mid, r;
    split(root, L, l, mid);
    split(mid, R - L + 1, mid, r);
    ll ans = sum(mid);
    merge(root, l, mid);
    merge(root, root, r);
    return ans;
  }
};

/**
 * Uso[0]: Mo mo(n); mo.add_query(l, r); auto ans = mo.solve(add, remove, get_ans);
 * Algoritmo de Mo 1D para responder consultas de rango offline cuando add y remove son $O(1)$.
 * Ejemplo de los callbacks (contar distintos en [l, r]): int cnt[MAXV] = {}, d = 0;
 *   auto add = [&](int i){ if (!cnt[a[i]]++) d++; };
 *   auto remove = [&](int i){ if (!--cnt[a[i]]) d--; };
 *   auto get = [&](){ return (ll)d; };
 * Complejidad: $O((N + Q)\sqrt{N})$ movimientos.
 */
struct Mo {
  struct Query { int l, r, id, block; };
  int n, bs;
  vector<Query> qs;
  Mo(int n, int bs = 0) : n(n), bs(bs ? bs : max(1, (int)sqrt(n))) {}
  void add_query(int l, int r) { qs.pb({l, r, (int)qs.size(), l / bs}); }
  template <class Add, class Remove, class GetAns>
  vector<ll> solve(Add add, Remove remove, GetAns get_ans) {
    vector<ll> res(qs.size());
    sort(all(qs), [](const Query& a, const Query& b) {
      if (a.block != b.block) return a.block < b.block;
      return (a.block & 1) ? a.r > b.r : a.r < b.r; // odd-even para acelerar
    });
    int cur_l = 0, cur_r = -1;
    for (const Query& q : qs) {
      while (cur_l > q.l) add(--cur_l);
      while (cur_r < q.r) add(++cur_r);
      while (cur_l < q.l) remove(cur_l++);
      while (cur_r > q.r) remove(cur_r--);
      res[q.id] = get_ans();
    }
    return res;
  }
};

/**
 * Uso[0]: MoWithUpdates mo(n); mo.add_update(pos, antes, despues); mo.add_query(l, r); auto ans = mo.solve(add, remove, apply, revert, get_ans);
 * Algoritmo de Mo con actualizaciones puntuales (Mo 3D); responde consultas sobre el arreglo en un instante dado.
 * apply(u, L, R) aplica el update u al rango actual [L,R] y revert lo deshace; L y R son los parametros que recibe el callback.
 * Ejemplo (distintos): si u.pos cae en [L,R] se retira el valor viejo y se agrega el nuevo; luego arr[u.pos] = u.new_val (y al reves en revert).
 * Complejidad: $O((N + Q)^{5/3})$ movimientos esperados.
 */
struct MoWithUpdates {
  struct Query { int l, r, t, id; };
  struct Update { int pos; ll old_val, new_val; };
  int n, bs;
  vector<Query> qs;
  vector<Update> ups;
  MoWithUpdates(int n, int bs = 0) : n(n), bs(bs ? bs : max(1, (int)pow(n, 2.0 / 3.0))) {}
  void add_update(int pos, ll old_val, ll new_val) { ups.pb({pos, old_val, new_val}); }
  void add_query(int l, int r) { qs.pb({l, r, (int)ups.size(), (int)qs.size()}); }
  template <class Add, class Remove, class Apply, class Revert, class GetAns>
  vector<ll> solve(Add add, Remove remove, Apply apply, Revert revert, GetAns get_ans) {
    vector<ll> res(qs.size());
    sort(all(qs), [&](const Query& a, const Query& b) {
      int ba = a.l / bs, bb = b.l / bs;
      if (ba != bb) return ba < bb;
      int ca = a.r / bs, cb = b.r / bs;
      if (ca != cb) return ca < cb;
      return a.t < b.t;
    });
    int cur_l = 0, cur_r = -1, cur_t = 0;
    for (const Query& q : qs) {
      while (cur_t < q.t) { apply(ups[cur_t], cur_l, cur_r); cur_t++; }
      while (cur_t > q.t) { cur_t--; revert(ups[cur_t], cur_l, cur_r); }
      while (cur_l > q.l) add(--cur_l);
      while (cur_r < q.r) add(++cur_r);
      while (cur_l < q.l) remove(cur_l++);
      while (cur_r > q.r) remove(cur_r--);
      res[q.id] = get_ans();
    }
    return res;
  }
};

/**
 * Uso[0]: RollbackMo mo(n); mo.add_query(l, r); mo.solve(add, current, snapshot, rollback, reset);
 * Variante de Mo con rollback para estructuras donde agregar es $O(1)$ pero eliminar es costoso o imposible.
 * Ejemplo (distintos): add(i) inserta a[i] y guarda i en un stack; snapshot() limpia el stack; rollback() deshace los adds del stack; reset() vacia toda la estructura.
 * Complejidad: $O((N + Q)\sqrt{N})$ operaciones de add.
 */
struct RollbackMo {
  struct Query { int l, r, id; };
  int n, bs;
  vector<Query> qs;
  RollbackMo(int n, int bs = 0) : n(n), bs(bs ? bs : max(1, (int)sqrt(n))) {}
  void add_query(int l, int r) { qs.pb({l, r, (int)qs.size()}); }
  template <class Add, class Current, class Snapshot, class Rollback, class Reset>
  vector<ll> solve(Add add, Current current, Snapshot snapshot, Rollback rollback, Reset reset) {
    vector<ll> res(qs.size());
    sort(all(qs), [&](const Query& a, const Query& b) {
      if (a.l / bs != b.l / bs) return a.l / bs < b.l / bs;
      return a.r < b.r;
    });
    int block = -1, cur_r = -1;
    for (const Query& q : qs) {
      int b = q.l / bs;
      int block_end = min(n - 1, (b + 1) * bs - 1);
      if (b != block) {
        block = b;
        cur_r = block_end; // rango persistente vacio
        reset();
      }
      if (q.r <= block_end) {
        // query contenida en el bloque: brute force sobre estructura vacia
        snapshot();
        for (int i = q.l; i <= q.r; i++) add(i);
        res[q.id] = current();
        rollback();
      } else {
        while (cur_r < q.r) add(++cur_r);
        snapshot();
        for (int i = block_end; i >= q.l; i--) add(i);
        res[q.id] = current();
        rollback();
      }
    }
    return res;
  }
};
