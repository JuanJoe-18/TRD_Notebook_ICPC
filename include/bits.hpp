/**
 * Uso: check_bit(m, i), set_bit(m, i), clear_bit(m, i), toggle_bit(m, i), lsb(m). i 0-based.
 * Manipulacion basica de bits de una mascara: consultar, activar, limpiar, alternar y aislar el bit menos significativo.
 * Complejidad: $O(1)$ por operacion.
 */
inline bool check_bit(ll mask, int i) { return (mask >> i) & 1; }
inline ll set_bit(ll mask, int i) { return mask | (1LL << i); }
inline ll clear_bit(ll mask, int i) { return mask & ~(1LL << i); }
inline ll toggle_bit(ll mask, int i) { return mask ^ (1LL << i); }
inline ll lsb(ll mask) { return mask & -mask; }
/**
 * Uso: count_bits(x), count_bits_ll(x), highest_bit(x) = bit mas alto (-1 si 0), ctz(x) = ceros a la derecha.
 * Conteo de bits activos, posicion del bit mas alto y numero de ceros finales con instrucciones nativas.
 * Complejidad: $O(1)$ por operacion.
 */
inline int count_bits(int mask) { return __builtin_popcount(mask); }
inline int count_bits_ll(ll mask) { return __builtin_popcountll(mask); }
inline int highest_bit(ll mask) { return mask == 0 ? -1 : 63 - __builtin_clzll(mask); }
inline int ctz(ll mask) { return mask == 0 ? 64 : __builtin_ctzll(mask); }

/**
 * Uso: iterate_submasks(m);
 * Itera todas las sub-mascaras de una mascara dada en orden decreciente (sin incluir la vacia).
 * Complejidad: $O(3^N)$ en total para una mascara de $N$ bits.
 */
inline void iterate_submasks(ll mask) {
  for (ll sub = mask; sub; sub = (sub - 1) & mask) { /* procesar sub */ }
}

/**
 * Uso: auto masks = gospers_hack(k, n);
 * Genera todas las mascaras de tamano $n$ con exactamente $k$ bits activos (gospers hack).
 * Complejidad: $O(\binom{n}{k})$ mascaras generadas, $O(1)$ amortizado por mascara.
 */
inline vector<ll> gospers_hack(int k, int n) {
  vector<ll> res;
  if (!k) return {0};
  ll mask = (1LL << k) - 1, limit = 1LL << n;
  while (mask < limit) {
    res.pb(mask);
    ll c = mask & -mask, r = mask + c;
    mask = (((r ^ mask) >> 2) / c) | r;
  }
  return res;
}

/**
 * Uso: XorBasis b; b.insert(x); b.can_form(k); b.get_max(init); b.get_min();
 * Base lineal de xor sobre enteros de 64 bits: maximo xor de un subconjunto, minimo y test de representabilidad.
 * Complejidad: insercion y consultas $O(60)$.
 */
struct XorBasis {
  static const int BITS = 60;
  ll basis[BITS] = {};
  int sz = 0;
  // Inserta; true si fue linealmente independiente
  bool insert(ll x) {
    for (int i = BITS - 1; i >= 0; i--) if (check_bit(x, i)) {
      if (!basis[i]) return basis[i] = x, sz++, true;
      x ^= basis[i];
    }
    return false;
  }
  bool can_form(ll k) {
    for (int i = BITS - 1; i >= 0; i--) if (check_bit(k, i)) {
      if (!basis[i]) return false;
      k ^= basis[i];
    }
    return true;
  }
  // Maximo de (initial ^ subset)
  ll get_max(ll initial = 0) {
    for (int i = BITS - 1; i >= 0; i--) if ((initial ^ basis[i]) > initial) initial ^= basis[i];
    return initial;
  }
  ll get_min() {
    for (int i = 0; i < BITS; i++) if (basis[i]) return basis[i];
    return 0;
  }
};

/**
 * Uso: BitTrie t; t.insert(x); t.max_xor(x); t.insert(x, -1) borra.
 * Trie binario para consultar el maximo xor de un valor contra el conjunto insertado, con soporte de borrado.
 * Complejidad: insercion y consulta $O(B)$ con $B = 30$ bits (60 si los valores llegan a $10^{18}$).
 */
struct BitTrie {
  static const int BITS = 30; // 60 si los numeros llegan a 1e18
  struct Node { int next[2] = {-1, -1}; int cnt = 0; };
  vector<Node> t;
  BitTrie() { t.emplace_back(); }
  void insert(ll x, int val = 1) { // val=-1 para borrar
    int u = 0;
    for (int i = BITS - 1; i >= 0; i--) {
      int b = check_bit(x, i);
      if (t[u].next[b] == -1) t[u].next[b] = t.size(), t.emplace_back();
      u = t[u].next[b];
      t[u].cnt += val;
    }
  }
  ll max_xor(ll x) {
    int u = 0; ll ans = 0;
    for (int i = BITS - 1; i >= 0; i--) {
      int b = check_bit(x, i), opp = b ^ 1;
      int v = t[u].next[opp];
      if (v != -1 && t[v].cnt > 0) ans |= 1LL << i, u = v;
      else u = t[u].next[b];
    }
    return ans;
  }
};