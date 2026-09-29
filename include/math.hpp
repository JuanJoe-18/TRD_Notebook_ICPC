namespace MathAlgo {
  /**
   * Uso: modpow(a, b, mod), modinv(a, mod), ext_gcd(a, b, x, y);
   * Potencia modular, inverso modular y euclides extendido ($ax + by = \gcd(a,b)$).
   * Complejidad: $O(\log b)$ por operacion.
   */
  ll modpow(ll a, ll b, ll m = MOD) {
    ll r = 1; a %= m;
    while (b) { if (b & 1) r = r * a % m; a = a * a % m; b >>= 1; }
    return r;
  }
  // Euclides extendido: resuelve ax + by = gcd(a,b)
  ll ext_gcd(ll a, ll b, ll& x, ll& y) {
    if (!b) { x = 1; y = 0; return a; }
    ll x1, y1, d = ext_gcd(b, a % b, x1, y1);
    x = y1; y = x1 - y1 * (a / b); return d;
  }
  /**
   * Uso: auto [x, lcm] = modinv(a, mod);
   * Inverso modular de $a$ modulo $m$ (exige $\gcd(a, m) = 1$).
   * Complejidad: $O(\log m)$.
   */
  ll modinv(ll a, ll m = MOD) {
    ll x, y; assert(ext_gcd(a, m, x, y) == 1);
    return (x % m + m) % m;
  }
  /**
   * Uso: auto [x, lcm] = crt(a1, m1, a2, m2);
   * Teorema chino del resto para dos congruencias; devuelve $\{-1,-1\}$ si no tienen solucion.
   * Complejidad: $O(\log \min(m_1, m_2))$.
   */
  pair<ll, ll> crt(ll a1, ll m1, ll a2, ll m2) {
    ll p, q, g = ext_gcd(m1, m2, p, q);
    if (a1 % g != a2 % g) return {-1, -1};
    ll lcm = (m1 / g) * m2;
    ll res = a1 + (ll)((__int128_t)(a2 - a1) / g * p % (m2 / g)) * m1;
    return {(res % lcm + lcm) % lcm, lcm};
  }

  /**
   * Uso: auto [x, lcm] = solve_crt_array(A, M);
   * Teorema chino del resto para un arreglo de congruencias; devuelve $\{-1,-1\}$ si alguna es inconsistente.
   * Complejidad: $O(K \log m)$ con $K$ el numero de congruencias.
   */
  pair<ll, ll> solve_crt_array(const vector<ll>& A, const vector<ll>& M) {
      ll ans = A[0];
      ll lcm = M[0];
      for (int i = 1; i < A.size(); i++) {
          auto res = crt(ans, lcm, A[i], M[i]);
          if (res.first == -1) return {-1, -1}; // No hay solución
          ans = res.first;
          lcm = res.second;
      }
      return {ans, lcm};
  }

  /**
   * Uso: Sieve sv(n); sv.primes; sv.spf; sv.mu; sv.phi; sv.factorize(x);
   * Criba lineal con lista de primos, menor factor primo, funcion de Mobius y totient de Euler.
   * Complejidad: criba $O(N)$, factorizacion de un valor $O(\log x)$.
   */
  struct Sieve {
    vector<int> primes, spf, mu, phi;
    Sieve(int n) {
      spf.assign(n + 1, 0); mu.assign(n + 1, 0); phi.assign(n + 1, 0);
      mu[1] = phi[1] = 1;
      for (int i = 2; i <= n; i++) {
        if (!spf[i]) spf[i] = i, primes.pb(i), mu[i] = -1, phi[i] = i - 1;
        for (int p : primes) {
          if (p > spf[i] || 1LL * i * p > n) break;
          spf[i * p] = p;
          if (i % p == 0) { mu[i * p] = 0; phi[i * p] = phi[i] * p; break; }
          mu[i * p] = -mu[i]; phi[i * p] = phi[i] * (p - 1);
        }
      }
    }
    vector<int> factorize(int x) { // factores primos con multiplicidad
      vector<int> r;
      while (x > 1) r.pb(spf[x]), x /= spf[x];
      return r;
    }
  };

  /**
   * Uso: Combinatorics c(n, mod); c.nCr(n, r); c.catalan(n);
   * Factoriales e inversos para combinar y numeros de Catalan modulo primo.
   * Complejidad: preproceso $O(N)$, consulta $O(1)$.
   */
  struct Combinatorics {
    vector<ll> fact, invfact;
    Combinatorics(int n, ll mod = MOD) {
      fact.assign(n + 1, 1); invfact.assign(n + 1, 1);
      for (int i = 1; i <= n; i++) fact[i] = fact[i - 1] * i % mod;
      invfact[n] = modinv(fact[n], mod);
      for (int i = n - 1; i >= 0; i--) invfact[i] = invfact[i + 1] * (i + 1) % mod;
      this->mod = mod;
    }
    ll mod;
    ll nCr(int n, int r) {
      if (r < 0 || r > n) return 0;
      return fact[n] * invfact[r] % mod * invfact[n - r] % mod;
    }
    ll catalan(int n) { return nCr(2 * n, n) * modinv(n + 1, mod) % mod; }
  };

  /**
   * Uso: gauss(matriz_extendida, soluciones);
   * Eliminacion gaussiana con pivoteo para sistemas lineales; devuelve 1, INF o 0 segun el numero de soluciones.
   * Complejidad: $O(n \cdot m^2)$.
   */
  int gauss(vector<vector<double>> a, vector<double>& ans) {
    int n = a.size(), m = a[0].size() - 1;
    vector<int> where(m, -1);
    for (int col = 0, row = 0; col < m && row < n; col++) {
      int sel = row;
      for (int i = row; i < n; i++) if (fabs(a[i][col]) > fabs(a[sel][col])) sel = i;
      if (fabs(a[sel][col]) < EPS) continue;
      for (int i = col; i <= m; i++) swap(a[sel][i], a[row][i]);
      where[col] = row;
      for (int i = 0; i < n; i++) if (i != row) {
        double c = a[i][col] / a[row][col];
        for (int j = col; j <= m; j++) a[i][j] -= a[row][j] * c;
      }
      row++;
    }
    ans.assign(m, 0);
    for (int i = 0; i < m; i++) if (where[i] != -1) ans[i] = a[where[i]][m] / a[where[i]][i];
    for (int i = 0; i < n; i++) {
      double s = 0;
      for (int j = 0; j < m; j++) s += ans[j] * a[i][j];
      if (fabs(s - a[i][m]) > EPS) return 0;
    }
    for (int i = 0; i < m; i++) if (where[i] == -1) return INF;
    return 1;
  }

  /**
   * Uso: fft(vec_complejos, invert); auto c = multiply(a, b);
   * Transformada rapida de Fourier iterativa y multiplicacion exacta de polinomios con coeficientes enteros.
   * Complejidad: $O(N \log N)$ con $N$ potencia de dos.
   */
  using cd = complex<double>;
  // FFT iterativa, invert = true para transformada inversa
  void fft(vector<cd>& a, bool invert) {
    int n = a.size();
    for (int i = 1, j = 0; i < n; i++) {
      int bit = n >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      if (i < j) swap(a[i], a[j]);
    }
    for (int len = 2; len <= n; len <<= 1) {
      double ang = 2 * acos(-1.0) / len * (invert ? -1 : 1);
      cd wlen(cos(ang), sin(ang));
      for (int i = 0; i < n; i += len) {
        cd w(1);
        for (int j = 0; j < len / 2; j++) {
          cd u = a[i + j], v = a[i + j + len / 2] * w;
          a[i + j] = u + v; a[i + j + len / 2] = u - v; w *= wlen;
        }
      }
    }
    if (invert) for (cd& x : a) x /= n;
  }
  // Multiplicacion de polinomios (coeficientes enteros) O(n log n)
  vector<int> multiply(vector<int> const& a, vector<int> const& b) {
    vector<cd> fa(all(a)), fb(all(b));
    int n = 1;
    while (n < a.size() + b.size()) n <<= 1;
    fa.resize(n); fb.resize(n);
    fft(fa, false); fft(fb, false);
    for (int i = 0; i < n; i++) fa[i] *= fb[i];
    fft(fa, true);
    vector<int> res(n);
    for (int i = 0; i < n; i++) res[i] = round(fa[i].real());
    while (res.size() > 1 && !res.back()) res.pop_back();
    return res;
  }

  /**
   * Uso: ntt(vec, invert); auto c = multiply_mod(a, b);
   * Transformada numerica rapida (NTT) y multiplicacion de polinomios con coeficientes modulares.
   * Complejidad: $O(N \log N)$ con $N$ potencia de dos.
   */
  void ntt(vector<ll>& a, bool invert) {
    int n = a.size();
    for (int i = 1, j = 0; i < n; i++) {
      int bit = n >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      if (i < j) swap(a[i], a[j]);
    }
    for (int len = 2; len <= n; len <<= 1) {
      ll wlen = modpow(3, (MOD - 1) / len);
      if (invert) wlen = modinv(wlen);
      for (int i = 0; i < n; i += len) {
        ll w = 1;
        for (int j = 0; j < len / 2; j++) {
          ll u = a[i + j], v = a[i + j + len / 2] * w % MOD;
          a[i + j] = (u + v) % MOD;
          a[i + j + len / 2] = (u - v + MOD) % MOD;
          w = w * wlen % MOD;
        }
      }
    }
    if (invert) {
      ll ninv = modinv(n);
      for (ll& x : a) x = x * ninv % MOD;
    }
  }
  vector<ll> multiply_mod(vector<ll> const& a, vector<ll> const& b) {
    vector<ll> fa(all(a)), fb(all(b));
    int n = 1;
    while (n < a.size() + b.size()) n <<= 1;
    fa.resize(n); fb.resize(n);
    ntt(fa, false); ntt(fb, false);
    for (int i = 0; i < n; i++) fa[i] = fa[i] * fb[i] % MOD;
    ntt(fa, true);
    while (fa.size() > 1 && !fa.back()) fa.pop_back();
    return fa;
  }

  /**
   * Uso: nim_game(pilas), sum_n(n), sum_squares(n), geom_sum(a, r, n);
   * Decide el ganador del juego de Nim y calcula sumatorias cerradas (aritmetica, cuadrados y geometrica) modulo $M$.
   * Complejidad: $O(N)$ para Nim, $O(\log M)$ para las sumatorias.
   */
  bool nim_game(const vi& piles) {
    int x = 0;
    for (int p : piles) x ^= p;
    return x != 0;
  }
  ll sum_n(ll n, ll m = MOD) { n %= m; return n * (n + 1) % m * modinv(2, m) % m; }
  ll sum_squares(ll n, ll m = MOD) {
    n %= m;
    return n * (n + 1) % m * (2 * n + 1) % m * modinv(6, m) % m;
  }
  ll geom_sum(ll a, ll r, ll n, ll m = MOD) {
    if (r == 1) return (n % m) * (a % m) % m;
    ll num = (modpow(r, n, m) - 1 + m) % m;
    return (a % m) * num % m * modinv((r - 1 + m) % m, m) % m;
  }

  /**
   * Uso: Matrix m(r, c); m.mat; m * o; m.power(p); Matrix::identity(n);
   * Matriz con multiplicacion modular y exponenciacion para recurrencias lineales y conteo de caminos.
   * Complejidad: producto $O(r \cdot c^2)$, potencia $O(r^3 \log p)$.
   */
  struct Matrix {
    int r, c;
    vector<vector<ll>> mat;
    Matrix(int r, int c) : r(r), c(c), mat(r, vector<ll>(c, 0)) {}
    static Matrix identity(int n) {
      Matrix I(n, n);
      for (int i = 0; i < n; i++) I.mat[i][i] = 1;
      return I;
    }
    Matrix operator*(const Matrix& o) const {
      Matrix res(r, o.c);
      for (int i = 0; i < r; i++)
        for (int k = 0; k < c; k++) if (mat[i][k])
          for (int j = 0; j < o.c; j++)
            res.mat[i][j] = (res.mat[i][j] + mat[i][k] * o.mat[k][j]) % MOD;
      return res;
    }
    Matrix power(ll p) const {
      Matrix res = identity(r), base = *this;
      while (p) { if (p & 1) res = res * base; base = base * base; p >>= 1; }
      return res;
    }
  };

  /**
   * Uso: auto C = berlekamp_massey(secuencia);
   * Encuentra la recurrencia lineal minima de una secuencia; $C[0] \cdot S[n] + C[1] \cdot S[n-1] + \dots = 0$.
   * Complejidad: $O(N^2)$.
   */
  vector<ll> berlekamp_massey(const vector<ll>& s) {
    vector<ll> C{1}, B{1};
    int L = 0; ll m = 1, b = 1;
    for (int i = 0; i < s.size(); i++) {
      ll d = 0;
      for (int j = 0; j <= L; j++) d = (d + C[j] * s[i - j]) % MOD;
      if (d) {
        vector<ll> T = C;
        ll coef = MOD - d * modpow(b, MOD - 2) % MOD;
        while (C.size() <= B.size() + m) C.pb(0);
        for (int j = 0; j < B.size(); j++) C[j + m] = (C[j + m] + coef * B[j]) % MOD;
        if (2 * L <= i) L = i + 1 - L, B = T, b = d, m = 1;
        else m++;
      } else m++;
    }
    C.resize(L + 1);
    return C;
  }

  /**
   * Uso: PollardRho::is_prime(n); PollardRho::factorize(n);
   * Test de primalidad Miller-Rabin y factorizacion de enteros de hasta $2^{64}$ con Pollard-Rho.
   * Complejidad: $O(n^{1/4})$ esperado por factor encontrado.
   */
  struct PollardRho {
    using u64 = uint64_t; using u128 = __uint128_t;
    static u64 modpow(u64 b, u64 e, u64 m) {
      u64 r = 1; b %= m;
      while (e) { if (e & 1) r = (u128)r * b % m; b = (u128)b * b % m; e >>= 1; }
      return r;
    }
    static bool is_prime(u64 n) {
      if (n < 2) return false;
      for (u64 p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
        if (n == p) return true;
        if (n % p == 0) return false;
      }
      u64 d = n - 1; int s = 0;
      while (!(d & 1)) d >>= 1, s++;
      for (u64 a : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
        if (a >= n) continue;
        u64 x = modpow(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool comp = true;
        for (int r = 1; r < s; r++) {
          x = (u128)x * x % n;
          if (x == n - 1) { comp = false; break; }
        }
        if (comp) return false;
      }
      return true;
    }
    static u64 get_factor(u64 n) {
      if (!(n & 1)) return 2;
      if (is_prime(n)) return n;
      u64 x = 2, y = 2, d = 1, c = 1;
      auto f = [&](u64 x) { return (u64)(((u128)x * x % n + c) % n); };
      while (d == 1) {
        x = f(x); y = f(f(y));
        d = gcd(x > y ? x - y : y - x, n);
        if (d == n) { x = rand() % (n - 2) + 2; y = x; c = rand() % (n - 1) + 1; d = 1; }
      }
      return d;
    }
    static void factor(u64 n, map<u64, int>& f) {
      if (n == 1) return;
      if (is_prime(n)) { f[n]++; return; }
      u64 d = get_factor(n);
      factor(d, f); factor(n / d, f);
    }
    static map<u64, int> factorize(u64 n) { // {primo, exponente}
      map<u64, int> f;
      factor(n, f);
      return f;
    }
  };

  /**
   * Uso: Lucas l(p); l.solve(n, r);
   * Combinaciones modulo un primo pequeno con $n$ y $r$ potencialmente enormes (teorema de Lucas).
   * Complejidad: $O(p)$ de preproceso, $O(\log_p n)$ por consulta.
   */
  struct Lucas {
    Combinatorics comb;
    ll P;
    Lucas(ll p) : P(p), comb(p - 1, p) {}
    ll solve(ll n, ll r) {
      if (r < 0 || r > n) return 0;
      if (!r) return 1;
      if (n % P < r % P) return 0;
      return solve(n / P, r / P) * comb.nCr(n % P, r % P) % P;
    }
  };

  /**
   * Uso: fwht_xor(v, invert); auto c = xor_convolution(a, b); // igual para or_ y and_
   * Transformada de Walsh-Hadamard y convoluciones sobre operaciones de bits (xor, and, or) modulo $M$.
   * Complejidad: $O(N \log N)$.
   */
  void fwht_xor(vector<ll>& a, bool invert) {
    int n = a.size();
    for (int len = 2; len <= n; len <<= 1) {
      for (int i = 0; i < n; i += len) {
        for (int j = 0; j < len / 2; j++) {
          ll u = a[i + j], v = a[i + j + len / 2];
          a[i + j] = (u + v) % MOD;
          a[i + j + len / 2] = (u - v + MOD) % MOD;
        }
      }
    }
    if (invert) {
      ll inv = modinv(n);
      for (ll& x : a) x = x * inv % MOD;
    }
  }
  void fwht_or(vector<ll>& a, bool invert) {
    int n = a.size();
    for (int len = 2; len <= n; len <<= 1) {
      for (int i = 0; i < n; i += len) {
        for (int j = 0; j < len / 2; j++) {
          ll u = a[i + j], v = a[i + j + len / 2];
          if (invert) a[i + j + len / 2] = (v - u + MOD) % MOD;
          else a[i + j + len / 2] = (u + v) % MOD;
        }
      }
    }
  }
  void fwht_and(vector<ll>& a, bool invert) {
    int n = a.size();
    for (int len = 2; len <= n; len <<= 1) {
      for (int i = 0; i < n; i += len) {
        for (int j = 0; j < len / 2; j++) {
          ll u = a[i + j], v = a[i + j + len / 2];
          if (invert) a[i + j] = (u - v + MOD) % MOD;
          else a[i + j] = (u + v) % MOD;
        }
      }
    }
  }
  vector<ll> xor_convolution(vector<ll> a, vector<ll> b) {
    int n = 1;
    while (n < max(a.size(), b.size())) n <<= 1;
    a.resize(n); b.resize(n);
    fwht_xor(a, false); fwht_xor(b, false);
    for (int i = 0; i < n; i++) a[i] = a[i] * b[i] % MOD;
    fwht_xor(a, true);
    return a;
  }
  vector<ll> or_convolution(vector<ll> a, vector<ll> b) {
    int n = 1;
    while (n < max(a.size(), b.size())) n <<= 1;
    a.resize(n); b.resize(n);
    fwht_or(a, false); fwht_or(b, false);
    for (int i = 0; i < n; i++) a[i] = a[i] * b[i] % MOD;
    fwht_or(a, true);
    return a;
  }
  vector<ll> and_convolution(vector<ll> a, vector<ll> b) {
    int n = 1;
    while (n < max(a.size(), b.size())) n <<= 1;
    a.resize(n); b.resize(n);
    fwht_and(a, false); fwht_and(b, false);
    for (int i = 0; i < n; i++) a[i] = a[i] * b[i] % MOD;
    fwht_and(a, true);
    return a;
  }

  /**
   * Uso: floor_sum(n, m, a, b);
   * Calcula $\sum_{i=0}^{n-1} \lfloor (a \cdot i + b) / m \rfloor$.
   * Complejidad: $O(\log m)$.
   */
  ll floor_sum(ll n, ll m, ll a, ll b) {
    ll ans = 0;
    while (true) {
      if (a >= m) { ans += (n - 1) * n * (a / m) / 2; a %= m; }
      if (b >= m) { ans += n * (b / m); b %= m; }
      ll y_max = a * n + b;
      if (y_max < m) break;
      n = y_max / m;
      b = y_max % m;
      swap(m, a);
    }
    return ans;
  }
}