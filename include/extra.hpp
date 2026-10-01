#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

/**
 * Uso[0]: ordered_set<T> s; s.insert(x); s.order_of_key(x); *s.find_by_order(k);
 * Conjunto ordenado de pbds con consultas por rango (k-esimo menor y posicion).
 * Complejidad: $O(\log N)$ por operacion.
 */
template <typename T>
using ordered_set = __gnu_pbds::tree<T, __gnu_pbds::null_type, less<T>, __gnu_pbds::rb_tree_tag, __gnu_pbds::tree_order_statistics_node_update>;
/**
 * Uso[0]: ordered_multiset<T> s; s.insert({x, id});
 * Conjunto ordenado de pbds que admite duplicados emparejando cada valor con un id unico.
 * Complejidad: $O(\log N)$ por operacion.
 */
template <typename T>
using ordered_multiset = __gnu_pbds::tree<pair<T, int>, __gnu_pbds::null_type, less<pair<T, int>>, __gnu_pbds::rb_tree_tag, __gnu_pbds::tree_order_statistics_node_update>;

/**
 * Uso: unordered_map<ll, int, custom_hash> m;
 * Functor de hash anti-colision (splitmix64 con semilla aleatoria) para tablas hash de claves numericas.
 * Complejidad: $O(1)$ por operacion de hash.
 */
struct custom_hash {
  static uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
    x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
    return x ^ (x >> 31);
  }
  size_t operator()(uint64_t x) const {
    static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
    return splitmix64(x + FIXED_RANDOM);
  }
};

/**
 * Uso: MeetInTheMiddle::subset_sum(a, target);
 * Resuelve subset sum exacto con la tecnica meet in the middle para arreglos de hasta 40 elementos.
 * Complejidad: $O(N \cdot 2^{N/2})$ tiempo, $O(2^{N/2})$ espacio.
 */
struct MeetInTheMiddle {
  static bool subset_sum(const vll& a, ll target) {
    int n = a.size(), mid = n / 2;
    vll left, right;
    for (int i = 0; i < (1 << mid); i++) {
      ll s = 0;
      for (int j = 0; j < mid; j++) if (i >> j & 1) s += a[j];
      left.pb(s);
    }
    int n2 = n - mid;
    for (int i = 0; i < (1 << n2); i++) {
      ll s = 0;
      for (int j = 0; j < n2; j++) if (i >> j & 1) s += a[mid + j];
      right.pb(s);
    }
    sort(all(right));
    for (ll s : left) if (binary_search(all(right), target - s)) return true;
    return false;
  }
};

/**
 * Uso[0]: ZobristHash zh; auto pref = zh.build_pref(a); zh.query_range(pref, l, r);
 * Hash de prefijos xor para comparar multiconjuntos de rangos en $O(1)$.
 * Complejidad: construccion $O(N)$, consulta $O(1)$.
 */
struct ZobristHash {
  using u64 = uint64_t;
  mt19937_64 rng;
  unordered_map<int, u64> val_hash;
  ZobristHash(u64 seed = 1337) : rng(seed) { val_hash[0] = 0; }
  u64 get_hash(int x) {
    if (!val_hash.count(x)) val_hash[x] = rng();
    return val_hash[x];
  }
  vector<u64> build_pref(const vi& a) { // prefijos XOR
    vector<u64> pref(a.size() + 1, 0);
    for (int i = 0; i < a.size(); i++) pref[i + 1] = pref[i] ^ get_hash(a[i]);
    return pref;
  }
  u64 query_range(const vector<u64>& pref, int l, int r) { return pref[r + 1] ^ pref[l]; }
};

/**
 * Uso[1]: BitsetReachability br(n); br.add_edge(u, v); br.build(); br.can_reach(u, v);
 * Cierre transitivo (alcanzabilidad todos-pares) de un grafo dirigido usando bitsets de palabras de 64 bits; el tamano se fija en el constructor.
 * Indexación: nodos 1-indexados.
 * Complejidad: construccion $O(N^3 / 64)$, consulta $O(1)$.
 */
struct BitsetReachability {
  int n, W;
  vector<vector<unsigned long long>> reach;
  BitsetReachability(int n) : n(n), W((n + 63) / 64), reach(n + 1, vector<unsigned long long>(W, 0)) {}
  void add_edge(int u, int v) { reach[u][v >> 6] |= 1ULL << (v & 63); }
  void build() {
    for (int i = 1; i <= n; i++) reach[i][i >> 6] |= 1ULL << (i & 63);
    for (int k = 1; k <= n; k++)
      for (int i = 1; i <= n; i++)
        if ((reach[i][k >> 6] >> (k & 63)) & 1ULL)
          for (int w = 0; w < W; w++) reach[i][w] |= reach[k][w];
  }
  bool can_reach(int u, int v) { return (reach[u][v >> 6] >> (v & 63)) & 1ULL; }
};

/**
 * Uso[0]: auto [vals, ids] = compress_coords(a);
 * Compresion de coordenadas: ordena y elimina duplicados (sort + unique) y mapea cada valor a su rango consecutivo $0..K-1$.
 * Indexación: ids 0-indexados en orden creciente de valor.
 * Complejidad: $O(N \log N)$ tiempo, $O(N)$ espacio.
 */
template<typename T>
pair<vector<T>, vector<int>> compress_coords(const vector<T>& a) {
  vector<T> vals = a;
  sort(all(vals));
  vals.erase(unique(all(vals)), vals.end());
  vector<int> ids(a.size());
  for (int i = 0; i < (int)a.size(); i++)
    ids[i] = lower_bound(all(vals), a[i]) - vals.begin();
  return {vals, ids};
}