/**
 * Uso: to_lower(s), to_upper(s), trim(s), replace_all(s, "a", "b"), split(s, ','), join(v, " ");
 * Utilidades de procesamiento de strings: mayusculas/minusculas, recorte, reemplazo, division y union.
 * Complejidad: $O(N)$ por operacion.
 */
void to_lower(string& s) { transform(all(s), s.begin(), ::tolower); }
void to_upper(string& s) { transform(all(s), s.begin(), ::toupper); }
void trim(string& s) {
  s.erase(s.begin(), find_if(all(s), [](unsigned char c) { return !isspace(c); }));
  s.erase(find_if(rall(s), [](unsigned char c) { return !isspace(c); }).base(), s.end());
}
void replace_all(string& s, const string& old, const string& nw) {
  size_t pos = 0;
  while ((pos = s.find(old, pos)) != string::npos) s.replace(pos, old.size(), nw), pos += nw.size();
}
vector<string> split(const string& s, char delim) {
  vector<string> res;
  string item;
  istringstream iss(s);
  while (getline(iss, item, delim)) res.pb(item);
  return res;
}
string join(const vector<string>& v, const string& delim) {
  string res;
  for (int i = 0; i < v.size(); i++) {
    res += v[i];
    if (i + 1 < v.size()) res += delim;
  }
  return res;
}

/**
 * Uso: auto pi = prefix_function(s);
 * Funcion prefijo de KMP: $pi[i]$ es el prefijo propio mas largo que tambien es sufijo de $s[0..i]$.
 * Complejidad: $O(N)$.
 */
vi prefix_function(const string& s) {
  int n = s.size();
  vi pi(n, 0);
  for (int i = 1; i < n; i++) {
    int j = pi[i - 1];
    while (j > 0 && s[i] != s[j]) j = pi[j - 1];
    if (s[i] == s[j]) j++;
    pi[i] = j;
  }
  return pi;
}

/**
 * Uso: auto z = z_function(s);
 * Funcion Z: $z[i]$ es la longitud del prefijo de $s$ que coincide con el sufijo que inicia en $i$.
 * Complejidad: $O(N)$.
 */
vi z_function(const string& s) {
  int n = s.size(), l = 0, r = 0;
  vi z(n, 0);
  for (int i = 1; i < n; i++) {
    if (i <= r) z[i] = min(r - i + 1, z[i - l]);
    while (i + z[i] < n && s[z[i]] == s[i + z[i]]) z[i]++;
    if (i + z[i] - 1 > r) l = i, r = i + z[i] - 1;
  }
  return z;
}

/**
 * Uso: string res = remove_all_occurrences(s, p);
 * Elimina todas las apariciones de un patron $p$ en $s$, incluidas las encadenadas.
 * Complejidad: $O(N + M)$ con $N = |s|$ y $M = |p|$.
 */
string remove_all_occurrences(const string& s, const string& p) {
  vi pi = prefix_function(p);
  string res;
  vi state(s.size() + 1, 0);
  for (char c : s) {
    res += c;
    int j = state[res.size() - 1];
    while (j > 0 && p[j] != c) j = pi[j - 1];
    if (p[j] == c) j++;
    state[res.size()] = j;
    if (j == p.size()) res.resize(res.size() - p.size());
  }
  return res;
}

/**
 * Uso[0]: auto occ = rabin_karp(text, pattern);
 * Encuentra las posiciones donde comienza un patron usando hashing (rolling hash).
 * Complejidad: $O(T + S)$ con $T = |text|$ y $S = |pattern|$.
 */
vi rabin_karp(const string& text, const string& pattern) {
  const ll P = 31, M = 1e9 + 9;
  int S = pattern.size(), T = text.size();
  if (!S || S > T) return {};
  vll p_pow(max(S, T));
  p_pow[0] = 1;
  for (int i = 1; i < p_pow.size(); i++) p_pow[i] = p_pow[i - 1] * P % M;
  vll h(T + 1, 0);
  for (int i = 0; i < T; i++) h[i + 1] = (h[i] + (text[i] - 'a' + 1) * p_pow[i]) % M;
  ll hs = 0;
  for (int i = 0; i < S; i++) hs = (hs + (pattern[i] - 'a' + 1) * p_pow[i]) % M;
  vi occ;
  for (int i = 0; i + S - 1 < T; i++) {
    ll cur = (h[i + S] + M - h[i]) % M;
    if (cur == hs * p_pow[i] % M) occ.pb(i);
  }
  return occ;
}

/**
 * Uso: auto p = manacher(s);
 * Radios de los palindromos en cada posicion (incluyendo centros pares); el maximo da el palindromo mas largo.
 * Complejidad: $O(N)$.
 */
vi manacher(const string& s) {
  string t = "^#";
  for (char c : s) t += c, t += '#';
  t += '$';
  int n = t.size(), c = 0, r = 0;
  vi p(n, 0);
  for (int i = 1; i < n - 1; i++) {
    int mir = 2 * c - i;
    p[i] = r > i ? min(r - i, p[mir]) : 0;
    while (t[i + 1 + p[i]] == t[i - 1 - p[i]]) p[i]++;
    if (i + p[i] > r) c = i, r = i + p[i];
  }
  return vi(p.begin() + 1, p.end() - 1);
}

/**
 * Uso[0]: StringHash sh(s); sh.get_hash(l, r); StringHash sh2(s, true); // bases aleatorias
 * Hash doble de substrings para compararlos en O(1). Con random_base = true las
 * bases se eligen al azar en tiempo de ejecucion (inmune a colisiones forzadas
 * contra bases conocidas en jueces adversariales).
 * Complejidad: preproceso O(N), consulta O(1).
 */
struct StringHash {
  static constexpr ll MOD1 = 1000000007LL; // 1e9 + 7
  static constexpr ll MOD2 = 1000000009LL; // 1e9 + 9
  string s;
  int n;
  ll B1, B2;
  vll h1, h2, p1, p2;
  StringHash(const string& s, bool random_base = false)
    : s(s), n(s.size()), B1(31), B2(37), h1(n + 1, 0), h2(n + 1, 0), p1(n + 1, 1), p2(n + 1, 1) {
    if (random_base) {
      static mt19937_64 rrng(chrono::steady_clock::now().time_since_epoch().count());
      uniform_int_distribution<ll> d1(300, MOD1 - 2), d2(300, MOD2 - 2);
      B1 = d1(rrng); B2 = d2(rrng);
    }
    for (int i = 0; i < n; i++) {
      h1[i + 1] = (h1[i] * B1 + (s[i] - 'a' + 1)) % MOD1;
      h2[i + 1] = (h2[i] * B2 + (s[i] - 'a' + 1)) % MOD2;
      p1[i + 1] = p1[i] * B1 % MOD1;
      p2[i + 1] = p2[i] * B2 % MOD2;
    }
  }
  uint64_t get_hash(int l, int r) {
    ll x1 = (h1[r + 1] - h1[l] * p1[r - l + 1] % MOD1 + MOD1) % MOD1;
    ll x2 = (h2[r + 1] - h2[l] * p2[r - l + 1] % MOD2 + MOD2) % MOD2;
    return (uint64_t)x1 << 32 | (uint32_t)x2;
  }
};

/**
 * Uso: Trie t; t.insert(s); t.search(s);
 * Trie (arbol de prefijos) sobre el alfabeto de minusculas con insercion y busqueda de palabras.
 * Complejidad: insercion y busqueda $O(L)$ con $L$ la longitud de la cadena.
 */
struct Trie {
  struct Node { int next[26] = {}; bool is_end = false; int count = 0; };
  vector<Node> t;
  Trie() { t.emplace_back(); }
  void insert(const string& s) {
    int v = 0;
    for (char ch : s) {
      int c = ch - 'a';
      if (!t[v].next[c]) t[v].next[c] = t.size(), t.emplace_back();
      v = t[v].next[c];
      t[v].count++;
    }
    t[v].is_end = true;
  }
  bool search(const string& s) {
    int v = 0;
    for (char ch : s) {
      int c = ch - 'a';
      if (!t[v].next[c]) return false;
      v = t[v].next[c];
    }
    return t[v].is_end;
  }
};

/**
 * Uso[0]: AhoCorasick ac; ac.insert(p); ac.build(); ac.count_occurrences(text); ac.first_occurrences(text);
 * Automata Aho-Corasick para buscar multiples patrones en un texto: conteo y primera aparicion de cada uno.
 * Complejidad: construccion $O(\sum |p_i| \cdot 26)$, procesamiento del texto $O(T)$.
 */
struct AhoCorasick {
  struct Node {
    int next[26] = {};
    int link = 0;
    vi word_ids;
  };
  vector<Node> t;
  vi bfs_order, pat_len;
  int num_patterns = 0;
  AhoCorasick() { t.emplace_back(); }
  int insert(const string& s) {
    int v = 0, id = num_patterns++;
    pat_len.pb(s.size());
    for (char ch : s) {
      int c = ch - 'a';
      if (!t[v].next[c]) t[v].next[c] = t.size(), t.emplace_back();
      v = t[v].next[c];
    }
    t[v].word_ids.pb(id);
    return id;
  }
  void build() {
    queue<int> q;
    for (int c = 0; c < 26; c++) {
      if (t[0].next[c]) q.push(t[0].next[c]);
      else t[0].next[c] = 0;
    }
    while (!q.empty()) {
      int u = q.front(); q.pop();
      bfs_order.pb(u);
      for (int c = 0; c < 26; c++) {
        int v = t[u].next[c];
        if (v) t[v].link = t[t[u].link].next[c], q.push(v);
        else t[u].next[c] = t[t[u].link].next[c];
      }
    }
  }
  // Numero de ocurrencias de cada patron en text
  vi count_occurrences(const string& text) {
    vi freq(t.size(), 0);
    int u = 0;
    for (char ch : text) u = t[u].next[ch - 'a'], freq[u]++;
    for (int i = bfs_order.size() - 1; i >= 0; i--) {
      int cur = bfs_order[i];
      freq[t[cur].link] += freq[cur];
    }
    vi ans(num_patterns, 0);
    for (int i = 1; i < t.size(); i++) if (freq[i]) for (int id : t[i].word_ids) ans[id] += freq[i];
    return ans;
  }
  // Posicion de la primera ocurrencia de cada patron (o -1)
  vi first_occurrences(const string& text) {
    vi first(t.size(), INF);
    int u = 0;
    for (int i = 0; i < text.size(); i++) {
      u = t[u].next[text[i] - 'a'];
      first[u] = min(first[u], i);
    }
    for (int i = bfs_order.size() - 1; i >= 0; i--) {
      int cur = bfs_order[i];
      first[t[cur].link] = min(first[t[cur].link], first[cur]);
    }
    vi ans(num_patterns, -1);
    for (int i = 1; i < t.size(); i++) if (first[i] != INF)
      for (int id : t[i].word_ids) {
        int st = first[i] - pat_len[id] + 1;
        ans[id] = ans[id] == -1 ? st : min(ans[id], st);
      }
    return ans;
  }
};

/**
 * Uso[0]: SuffixAutomaton sam(s); sam.distinct_substrings(); sam.count_occurrences(p); sam.first_occurrence(p);
 * Automata de sufijos: substrings distintos, ocurrencias y primera aparicion de un patron.
 * Complejidad: construccion $O(N \cdot 26)$, consultas $O(|p|)$.
 */
struct SuffixAutomaton {
  struct Node {
    int len = 0, link = -1, next[26] = {};
    int firstpos = -1;
    ll cnt = 0;
  };
  vector<Node> t;
  int last;
  bool occ_computed = false;
  SuffixAutomaton(const string& s = "") {
    t.reserve(s.size() * 2 + 1);
    t.emplace_back();
    last = 0;
    for (char c : s) extend(c - 'a');
  }
  void extend(int c) {
    int cur = t.size();
    t.emplace_back();
    t[cur].len = t[last].len + 1;
    t[cur].firstpos = t[cur].len - 1;
    t[cur].cnt = 1;
    int p = last;
    while (p != -1 && !t[p].next[c]) t[p].next[c] = cur, p = t[p].link;
    if (p == -1) t[cur].link = 0;
    else {
      int q = t[p].next[c];
      if (t[p].len + 1 == t[q].len) t[cur].link = q;
      else {
        int clone = t.size();
        t.push_back(t[q]);
        t[clone].len = t[p].len + 1;
        t[clone].cnt = 0;
        while (p != -1 && t[p].next[c] == q) t[p].next[c] = clone, p = t[p].link;
        t[q].link = t[cur].link = clone;
      }
    }
    last = cur;
  }
  void compute_occurrences() {
    if (occ_computed) return;
    occ_computed = true;
    int n = t.size();
    vi c(n + 1, 0), order(n);
    for (int i = 0; i < n; i++) c[t[i].len]++;
    for (int i = 1; i <= n; i++) c[i] += c[i - 1];
    for (int i = 0; i < n; i++) order[--c[t[i].len]] = i;
    for (int i = n - 1; i > 0; i--) {
      int u = order[i];
      if (t[u].link != -1) t[t[u].link].cnt += t[u].cnt;
    }
  }
  ll distinct_substrings() {
    ll ans = 0;
    for (int i = 1; i < t.size(); i++) ans += t[i].len - t[t[i].link].len;
    return ans;
  }
  int first_occurrence(const string& p) {
    int u = 0;
    for (char ch : p) {
      int c = ch - 'a';
      if (!t[u].next[c]) return -1;
      u = t[u].next[c];
    }
    return t[u].firstpos - p.size() + 1;
  }
  ll count_occurrences(const string& p) {
    compute_occurrences();
    int u = 0;
    for (char ch : p) {
      int c = ch - 'a';
      if (!t[u].next[c]) return 0;
      u = t[u].next[c];
    }
    return t[u].cnt;
  }
};

/**
 * Uso[0]: SuffixArray sa(s); sa.p (sufijos ordenados); sa.get_lcp(i, j);
 * Arreglo de sufijos con tabla RMQ de lcp para consultar el prefijo comun mas largo de dos sufijos en $O(1)$.
 * Complejidad: construccion $O(N \log N)$, consulta lcp $O(1)$.
 */
struct SuffixArray {
  string s;
  int n, log_n;
  vi p, c, lcp, rank;
  vector<vi> st;
  SuffixArray(string s) : s(s + "$"), n(s.size()) {
    build_sa(); build_lcp(); build_rmq();
  }
  void build_sa() {
    const int alphabet = 256;
    p.assign(n, 0); c.assign(n, 0);
    vi cnt(max(alphabet, n), 0), p_new(n), c_new(n);
    for (int i = 0; i < n; i++) cnt[s[i]]++;
    for (int i = 1; i < alphabet; i++) cnt[i] += cnt[i - 1];
    for (int i = 0; i < n; i++) p[--cnt[s[i]]] = i;
    c[p[0]] = 0;
    int classes = 1;
    for (int i = 1; i < n; i++) {
      if (s[p[i]] != s[p[i - 1]]) classes++;
      c[p[i]] = classes - 1;
    }
    for (int k = 0; (1 << k) < n; k++) {
      for (int i = 0; i < n; i++) {
        p_new[i] = p[i] - (1 << k);
        if (p_new[i] < 0) p_new[i] += n;
      }
      fill(cnt.begin(), cnt.begin() + classes, 0);
      for (int i = 0; i < n; i++) cnt[c[p_new[i]]]++;
      for (int i = 1; i < classes; i++) cnt[i] += cnt[i - 1];
      for (int i = n - 1; i >= 0; i--) p[--cnt[c[p_new[i]]]] = p_new[i];
      c_new[p[0]] = 0;
      classes = 1;
      for (int i = 1; i < n; i++) {
        pair<int, int> cur = {c[p[i]], c[(p[i] + (1 << k)) % n]};
        pair<int, int> prev = {c[p[i - 1]], c[(p[i - 1] + (1 << k)) % n]};
        if (cur != prev) classes++;
        c_new[p[i]] = classes - 1;
      }
      c.swap(c_new);
    }
  }
  void build_lcp() {
    lcp.assign(n, 0);
    rank.assign(n, 0);
    for (int i = 0; i < n; i++) rank[p[i]] = i;
    int k = 0;
    for (int i = 0; i < n - 1; i++) {
      int pi = rank[i], j = p[pi - 1];
      while (s[i + k] == s[j + k]) k++;
      lcp[pi] = k;
      if (k) k--;
    }
  }
  void build_rmq() {
    log_n = log2(n) + 1;
    st.assign(n, vi(log_n));
    for (int i = 0; i < n; i++) st[i][0] = lcp[i];
    for (int j = 1; j < log_n; j++)
      for (int i = 0; i + (1 << j) <= n; i++)
        st[i][j] = min(st[i][j - 1], st[i + (1 << (j - 1))][j - 1]);
  }
  // LCP de los sufijos que empiezan en i y j
  int get_lcp(int i, int j) {
    if (i == j) return n - 1 - i;
    int u = rank[i], v = rank[j];
    if (u > v) swap(u, v);
    u++;
    int k = log2(v - u + 1);
    return min(st[u][k], st[v - (1 << k) + 1][k]);
  }
};

/**
 * Uso[0]: PalindromicTree pt(s); pt.count_palindromes(); pt.compute_occurrences();
 * Arbol de palindromos (eertree): todos los palindromos distintos, sus ocurrencias y su primera posicion.
 * Complejidad: construccion $O(N \cdot 26)$, ocurrencias $O(N)$.
 */
struct PalindromicTree {
  struct Node {
    int next[26] = {};
    int len = 0, link = 0;
    ll occ = 0;
    int first_pos = -1; // posicion donde aparece por primera vez
  };
  vector<Node> t;
  string s; // s[0] = sentinel, nunca coincide con una letra
  int last;
  PalindromicTree(const string& str = "") {
    t.emplace_back(); t.emplace_back();
    t[0].len = -1; t[0].link = 0; // raiz imaginaria (len -1)
    t[1].len = 0;  t[1].link = 0; // raiz del vacio (len 0)
    last = 1;
    s = "$";
    for (char c : str) add(c);
  }
  int get_link(int v, int pos) {
    while (true) {
      int idx = pos - t[v].len - 1;
      if (idx >= 0 && s[idx] == s[pos]) break;
      v = t[v].link;
    }
    return v;
  }
  void add(char c) {
    s += c;
    int pos = s.size() - 1;
    int cur = get_link(last, pos);
    int idx = c - 'a';
    if (!t[cur].next[idx]) {
      int now = t.size();
      t.emplace_back();
      t[now].len = t[cur].len + 2;
      t[now].first_pos = pos - t[now].len; // 0-indexado
      t[now].link = t[now].len == 1 ? 1 : t[get_link(t[cur].link, pos)].next[idx];
      t[cur].next[idx] = now;
    }
    last = t[cur].next[idx];
    t[last].occ++;
  }
  ll count_palindromes() { return t.size() - 2; } // palindromos distintos
  void compute_occurrences() { // occ[u] = veces que aparece el palindromo u
    for (int i = t.size() - 1; i > 1; i--) t[t[i].link].occ += t[i].occ;
  }
};
