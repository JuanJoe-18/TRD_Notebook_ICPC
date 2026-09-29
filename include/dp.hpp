/**
 * Uso: LiChaoTree cht; cht.add(m, c); ll v = cht.query(x);
 * Arbol de Li Chao para consultar el minimo de $m \cdot x + c$ entre las rectas insertadas (cambiar < por > para max).
 * Complejidad: insercion y consulta $O(\log R)$ con $R$ el tamano del rango de $x$.
 */
struct LiChaoTree {
  struct Line {
    ll m, c;
    Line(ll m = 0, ll c = LINF) : m(m), c(c) {}
    ll eval(ll x) const { return m * x + c; }
  };
  struct Node { Line line; int lc = -1, rc = -1; };
  vector<Node> st;
  ll min_x, max_x;
  LiChaoTree(ll mn = -1e9, ll mx = 1e9) : min_x(mn), max_x(mx) { st.emplace_back(); }
  void add_line(Line nw, int v, ll l, ll r) {
    ll mid = l + (r - l) / 2;
    bool left = nw.eval(l) < st[v].line.eval(l);
    bool midb = nw.eval(mid) < st[v].line.eval(mid);
    if (midb) swap(st[v].line, nw);
    if (l == r) return;
    if (left != midb) {
      if (st[v].lc == -1) st[v].lc = st.size(), st.emplace_back();
      add_line(nw, st[v].lc, l, mid);
    } else {
      if (st[v].rc == -1) st[v].rc = st.size(), st.emplace_back();
      add_line(nw, st[v].rc, mid + 1, r);
    }
  }
  void add(ll m, ll c) { add_line(Line(m, c), 0, min_x, max_x); }
  ll query(ll x, int v, ll l, ll r) {
    if (v == -1) return LINF;
    ll res = st[v].line.eval(x);
    if (l == r) return res;
    ll mid = l + (r - l) / 2;
    if (x <= mid) return min(res, query(x, st[v].lc, l, mid));
    return min(res, query(x, st[v].rc, mid + 1, r));
  }
  ll query(ll x) { return query(x, 0, min_x, max_x); }
};

/**
 * Uso: dnc_dp(1, n, 0, n, dp_prev, dp_curr, cost);
 * Divide & conquer DP: resuelve $dp[i][j] = \min_{k<j}(dp[i-1][k] + cost(k+1,j))$ cuando la decision optima es monotona.
 * Complejidad: $O(K N \log N)$ con $cost$ evaluada en $O(1)$.
 */
template<class CostFunc>
void dnc_dp(int l, int r, int optl, int optr, const vll& prev, vll& curr, CostFunc& cost) {
  if (l > r) return;
  int mid = (l + r) / 2;
  pair<ll, int> best = {LINF, -1};
  for (int k = optl; k <= min(mid - 1, optr); k++) {
    ll cur = prev[k] == LINF ? LINF : prev[k] + cost(k + 1, mid);
    if (cur < best.fi) best = {cur, k};
  }
  curr[mid] = best.fi;
  int opt = best.se == -1 ? optl : best.se;
  dnc_dp(l, mid - 1, optl, opt, prev, curr, cost);
  dnc_dp(mid + 1, r, opt, optr, prev, curr, cost);
}

/**
 * Uso: auto dp = knuth_dp(n, cost);
 * Optimizacion de Knuth para DP de intervalos $dp[i][j] = \min_{k}(dp[i][k] + dp[k][j]) + cost(i,j)$ con decision monotona.
 * Complejidad: $O(N^2)$.
 */
template<class CostFunc>
vector<vll> knuth_dp(int n, CostFunc& cost) {
  vector<vll> dp(n + 1, vll(n + 1, 0));
  vector<vi> opt(n + 1, vi(n + 1, 0));
  for (int i = 0; i < n; i++) opt[i][i + 1] = i;
  for (int len = 2; len <= n; len++)
    for (int i = 0; i + len <= n; i++) {
      int j = i + len;
      dp[i][j] = LINF;
      for (int k = opt[i][j - 1]; k <= min(j - 1, opt[i + 1][j]); k++) {
        ll cur = dp[i][k] + dp[k][j] + cost(i, j);
        if (cur < dp[i][j]) dp[i][j] = cur, opt[i][j] = k;
      }
    }
  return dp;
}

/**
 * Uso: DigitDP::solve(R) - DigitDP::solve(L-1);
 * DP sobre digitos para contar numeros en $[0, N]$ que cumplen condiciones; ajustar memo y transiciones segun el problema.
 * Complejidad: $O(\text{digitos} \times \text{estados})$ por llamada.
 */
namespace DigitDP {
  string s;
  ll memo[20][2][2]; // agregar estados segun problema
  ll dp(int idx, bool tight, bool lz) {
    if (idx == s.size()) return 1; // 1 si el numero es valido
    if (memo[idx][tight][lz] != -1) return memo[idx][tight][lz];
    int lim = tight ? s[idx] - '0' : 9;
    ll ans = 0;
    for (int d = 0; d <= lim; d++)
      ans += dp(idx + 1, tight && d == lim, lz && d == 0);
    return memo[idx][tight][lz] = ans;
  }
  ll solve(ll n) {
    if (n < 0) return 0;
    s = to_string(n);
    memset(memo, -1, sizeof(memo));
    return dp(0, true, true);
  }
}

/**
 * Uso: sos_dp(dp, n_bits);
 * Suma sobre subconjuntos: $dp[mask]$ acumula la suma de todos sus subconjuntos (sum over subsets).
 * Complejidad: $O(N \cdot 2^N)$.
 */
void sos_dp(vll& dp, int n) {
  for (int i = 0; i < n; i++)
    for (int mask = 0; mask < (1 << n); mask++)
      if (mask & (1 << i)) dp[mask] += dp[mask ^ (1 << i)];
}

/**
 * Uso: auto [len, seq] = lis(a);
 * Subsecuencia creciente mas larga estricta, con reconstruccion de la secuencia.
 * Complejidad: $O(N \log N)$.
 */
pair<int, vi> lis(const vi& a) {
  int n = a.size();
  vi d, pos, parent(n, -1);
  for (int i = 0; i < n; i++) {
    auto it = lower_bound(all(d), a[i]);
    int idx = it - d.begin();
    if (it == d.end()) d.pb(a[i]), pos.pb(i);
    else *it = a[i], pos[idx] = i;
    if (idx > 0) parent[i] = pos[idx - 1];
  }
  vi seq;
  for (int cur = pos.empty() ? -1 : pos.back(); cur != -1; cur = parent[cur]) seq.pb(a[cur]);
  reverse(all(seq));
  return {d.size(), seq};
}

/**
 * Uso: bitset_knapsack(pesos, target);
 * Determina si un subconjunto de pesos suma exactamente un valor objetivo mediante bitsets.
 * Complejidad: $O(N \cdot S / 64)$ con $S$ la suma maxima ($10^5$ fija en MAX_SUM).
 */
const int MAX_SUM = 1e5 + 5;
bool bitset_knapsack(const vi& weights, int target) {
  bitset<MAX_SUM> dp;
  dp[0] = 1;
  for (int w : weights) dp |= dp << w;
  return dp[target];
}

/**
 * Uso: ll ans = WQS::solve(check, K);
 * Busqueda binaria de Lagrange (Alien's trick) para optimizar un DP con la restriccion de usar exactamente $K$ elementos o particiones.
 * Complejidad: $O(\log R \cdot f)$ con $f$ el costo del DP por evaluacion y $R$ el rango de la penalidad.
 */
namespace WQS {
  template <class F>
  ll solve(F check, int K) {
    ll lo = -LINF, hi = LINF;
    while (lo + 1 < hi) {
      ll mid = lo + (hi - lo) / 2;
      if (check(mid).se >= K) lo = mid;
      else hi = mid;
    }
    auto [c, k] = check(lo);
    return c - lo * (ll)K; // quita la penalidad de los K elementos
  }
}

/**
 * Uso: CHTDeque cht; cht.add(m, c); cht.query(x); // y cht.query_bs(x) para x arbitraria
 * Convex hull trick con deque para minimos $m \cdot x + c$ con pendientes monotonas y consultas con $x$ creciente.
 * Complejidad: insercion y consulta amortizadas $O(1)$; consulta binaria $O(\log N)$.
 */
struct CHTDeque {
  struct Line { ll m, c; };
  deque<Line> dq;
  int slope_inc = 0; // +1 pendientes crecientes, -1 decrecientes (auto-detectado)
  static ll eval(const Line& l, ll x) { return l.m * x + l.c; }
  // true si l2 queda obsoleta entre l1 y l3
  bool is_bad(const Line& l1, const Line& l2, const Line& l3) {
    if (slope_inc > 0) return (__int128)(l1.c - l2.c) * (l3.m - l2.m) <= (__int128)(l2.c - l3.c) * (l2.m - l1.m);
    return (__int128)(l3.c - l1.c) * (l1.m - l2.m) <= (__int128)(l2.c - l1.c) * (l1.m - l3.m);
  }
  void add(ll m, ll c) {
    if (!dq.empty() && dq.back().m == m) {
      if (dq.back().c <= c) return; // la nueva no aporta (min)
      dq.pop_back();
    }
    if (dq.size() >= 2 && slope_inc == 0) slope_inc = dq.back().m < m ? 1 : -1;
    Line nw{m, c};
    while (dq.size() >= 2 && is_bad(dq[dq.size() - 2], dq.back(), nw)) dq.pop_back();
    dq.push_back(nw);
  }
  ll query(ll x) { // x debe ser no decreciente
    if (slope_inc > 0) { // pendientes crecientes: optimo va del fondo al frente
      while (dq.size() >= 2 && eval(dq.back(), x) >= eval(dq[dq.size() - 2], x)) dq.pop_back();
      return eval(dq.back(), x);
    }
    while (dq.size() >= 2 && eval(dq.front(), x) >= eval(dq[1], x)) dq.pop_front();
    return eval(dq.front(), x);
  }
  ll query_bs(ll x) { // x arbitraria
    int lo = 0, hi = (int)dq.size() - 1;
    while (lo < hi) {
      int mid = (lo + hi) / 2;
      if (eval(dq[mid], x) <= eval(dq[mid + 1], x)) hi = mid;
      else lo = mid + 1;
    }
    return eval(dq[lo], x);
  }
};