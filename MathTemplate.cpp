/**
 * @file MathTemplateICPC.cpp
 * @brief Plantilla de matemáticas para programación competitiva.
 * @details Teoría de Números, Combinatoria, FFT, Eliminación Gaussiana, Teoría de Juegos y Fórmulas.
 */
//  __  __    _    _____ _   _
// |  \/  |  / \  |_   _| | | |
// | |\/| | / _ \   | | | |_| |
// | |  | |/ ___ \  | | |  _  |
// |_|  |_/_/   \_\ |_| |_| |_|
//
//           MATH TEMPLATE

#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
using u64 = uint64_t;
using u128 = __uint128_t;

const ll MOD = 1e9 + 7;
const double PI = acos(-1.0);
const ll INF = 1e15;

// ==========================================
// ESTRUCTURA MAESTRA DE MATEMÁTICAS (ICPC)
// ==========================================

namespace MathAlgo {

    /**
     * 1. Exponenciación Binaria (O(log B))
     */
    ll modpow(ll a, ll b, ll m = MOD) {
        ll res = 1;
        a %= m;
        while (b > 0) {
            if (b & 1) res = res * a % m;
            a = a * a % m;
            b >>= 1;
        }
        return res;
    }

    /**
     * 2. Algoritmo de Euclides Extendido
     */
    ll extGCD(ll a, ll b, ll &x, ll &y) {
        if (b == 0) { x = 1; y = 0; return a; }
        ll x1, y1;
        ll d = extGCD(b, a % b, x1, y1);
        x = y1; y = x1 - y1 * (a / b);
        return d;
    }

    /**
     * 3. Inverso Modular (O(log M))
     */
    ll modinv(ll a, ll m = MOD) {
        ll x, y;
        ll g = extGCD(a, m, x, y);
        assert(g == 1);
        return (x % m + m) % m;
    }

    /**
     * 4. Teorema Chino del Resto (CRT)
     */
    pair<ll, ll> crt(ll a1, ll m1, ll a2, ll m2) {
        ll p, q;
        ll g = extGCD(m1, m2, p, q);
        if (a1 % g != a2 % g) return {-1, -1};
        ll lcm = (m1 / g) * m2;
        ll res = a1 + (ll)((__int128_t)(a2 - a1) / g * p % (m2 / g)) * m1;
        return {(res % lcm + lcm) % lcm, lcm};
    }

    /**
     * 5. Criba Factorizadora (O(N) precalculo, O(log N) factorización)
     */
/**
     * 5. Criba Extendida (Factorización, Möbius, Totiente de Euler)
     */
    struct Sieve {
        vector<int> primes, spf, mu, phi;
        Sieve(int n) {
            spf.assign(n + 1, 0);
            mu.assign(n + 1, 0);
            phi.assign(n + 1, 0);
            mu[1] = 1; phi[1] = 1;
            
            for (int i = 2; i <= n; i++) {
                if (spf[i] == 0) { 
                    spf[i] = i; 
                    primes.push_back(i);
                    mu[i] = -1;
                    phi[i] = i - 1;
                }
                for (int j = 0; j < primes.size() && primes[j] <= spf[i] && i * primes[j] <= n; j++) {
                    int p = primes[j];
                    spf[i * p] = p;
                    if (i % p == 0) {
                        mu[i * p] = 0;
                        phi[i * p] = phi[i] * p;
                        break;
                    } else {
                        mu[i * p] = -mu[i];
                        phi[i * p] = phi[i] * (p - 1);
                    }
                }
            }
        }
        vector<int> factorize(int x) {
            vector<int> res;
            while (x > 1) { res.push_back(spf[x]); x /= spf[x]; }
            return res;
        }
    };

    /**
     * 6. Combinatoria y Catalan
     */
    struct Combinatorics {
        vector<ll> fact, invfact;
        ll MOD;
        Combinatorics(int n, ll mod = 1e9 + 7) : MOD(mod) {
            fact.assign(n + 1, 1);
            invfact.assign(n + 1, 1);
            for (int i = 1; i <= n; i++) fact[i] = fact[i - 1] * i % MOD;
            invfact[n] = modinv(fact[n], MOD);
            for (int i = n - 1; i >= 0; i--) invfact[i] = invfact[i + 1] * (i + 1) % MOD;
        }
        ll nCr(int n, int r) {
            if (r < 0 || r > n) return 0;
            return fact[n] * invfact[r] % MOD * invfact[n - r] % MOD;
        }
        ll catalan(int n) {
            if (n < 0) return 0;
            return nCr(2 * n, n) * modinv(n + 1, MOD) % MOD;
        }
    };

    /**
     * 7. Eliminación Gaussiana en R (O(N^3))
     * @param a Matriz extendida [Nx(M+1)] donde la última columna es el vector de resultados B.
     * @param ans Vector donde se GUARDARÁN los valores de (x_0, x_1, ..., x_{m-1}).
     * @return 1 (Solución única), INF (Infinitas, 'ans' devuelve una válida fijando libres a 0), 0 (Sin solución).
     */
    int gauss(vector<vector<double>> a, vector<double> &ans) {
        int n = a.size();
        int m = a[0].size() - 1; // Número de variables
        const double EPS = 1e-9;
        vector<int> where(m, -1); // Dónde está el pivote de cada columna

        for (int col = 0, row = 0; col < m && row < n; ++col) {
            int sel = row;
            for (int i = row; i < n; ++i)
                if (abs(a[i][col]) > abs(a[sel][col])) sel = i;
            if (abs(a[sel][col]) < EPS) continue;

            for (int i = col; i <= m; ++i) swap(a[sel][i], a[row][i]);
            where[col] = row;

            for (int i = 0; i < n; ++i) {
                if (i != row) {
                    double c = a[i][col] / a[row][col];
                    for (int j = col; j <= m; ++j) a[i][j] -= a[row][j] * c;
                }
            }
            ++row;
        }

        // --- EXTRACCIÓN DE SOLUCIONES ---
        ans.assign(m, 0); // Inicializa todas las variables en 0
        for (int i = 0; i < m; ++i) {
            if (where[i] != -1) // Si la variable i tiene un pivote
                ans[i] = a[where[i]][m] / a[where[i]][i];
        }

        // Verificación de consistencia
        for (int i = 0; i < n; ++i) {
            double sum = 0;
            for (int j = 0; j < m; ++j) sum += ans[j] * a[i][j];
            if (abs(sum - a[i][m]) > EPS) return 0; // Contradicción (Ej: 0 = 5)
        }

        for (int i = 0; i < m; ++i)
            if (where[i] == -1) return INF; // Variables libres = Infinitas soluciones
        return 1; // Sistema compatible determinado
    }

    /**
     * 8. Transformada Rápida de Fourier (FFT) - O(N log N)
     */
    using cd = complex<double>;
    void fft(vector<cd> & a, bool invert) {
        int n = a.size();
        for (int i = 1, j = 0; i < n; i++) {
            int bit = n >> 1;
            for (; j & bit; bit >>= 1) j ^= bit;
            j ^= bit;
            if (i < j) swap(a[i], a[j]);
        }
        for (int len = 2; len <= n; len <<= 1) {
            double angle = 2 * PI / len * (invert ? -1 : 1);
            cd wlen(cos(angle), sin(angle));
            for (int i = 0; i < n; i += len) {
                cd w(1);
                for (int j = 0; j < len / 2; j++) {
                    cd u = a[i + j], v = a[i + j + len / 2] * w;
                    a[i + j] = u + v;
                    a[i + j + len / 2] = u - v;
                    w *= wlen;
                }
            }
        }
        if (invert) for (cd & x : a) x /= n;
    }

    vector<int> multiply(vector<int> const& a, vector<int> const& b) {
        vector<cd> fa(a.begin(), a.end()), fb(b.begin(), b.end());
        int n = 1;
        while (n < a.size() + b.size()) n <<= 1;
        fa.resize(n); fb.resize(n);
        fft(fa, false); fft(fb, false);
        for (int i = 0; i < n; i++) fa[i] *= fb[i];
        fft(fa, true);
        vector<int> result(n);
        for (int i = 0; i < n; i++) result[i] = round(fa[i].real());
        while(result.size() > 1 && result.back() == 0) result.pop_back();
        return result;
    }

    void ntt(vector<ll> &a, bool invert) {
        int n = a.size();
        for (int i = 1, j = 0; i < n; i++) {
            int bit = n >> 1;
            for (; j & bit; bit >>= 1) j ^= bit;
            j ^= bit;
            if (i < j) swap(a[i], a[j]);
        }
        for (int len = 2; len <= n; len <<= 1) {
            ll wlen = modpow(root, (MOD - 1) / len);
            if (invert) wlen = modinv(wlen);
            for (int i = 0; i < n; i += len) {
                ll w = 1;
                for (int j = 0; j < len / 2; j++) {
                    ll u = a[i+j], v = a[i+j+len/2] * w % MOD;
                    a[i+j] = u + v < MOD ? u + v : u + v - MOD;
                    a[i+j+len/2] = u - v >= 0 ? u - v : u - v + MOD;
                    w = w * wlen % MOD;
                }
            }
        }
        if (invert) {
            ll n_inv = modinv(n);
            for (ll &x : a) x = x * n_inv % MOD;
        }
    }

    vector<ll> multiply_mod(vector<ll> const& a, vector<ll> const& b) {
        vector<ll> fa(a.begin(), a.end()), fb(b.begin(), b.end());
        int n = 1;
        while (n < a.size() + b.size()) n <<= 1;
        fa.resize(n); fb.resize(n);
        ntt(fa, false); ntt(fb, false);
        for (int i = 0; i < n; i++) fa[i] = fa[i] * fb[i] % MOD;
        ntt(fa, true);
        while (fa.size() > 1 && fa.back() == 0) fa.pop_back();
        return fa;
    }

    /**
     * 9. Teoría de Juegos - Nim (XOR Sum)
     */
    bool nim_game(const vector<int>& piles) {
        int xor_sum = 0;
        for (int p : piles) xor_sum ^= p;
        return xor_sum != 0;
    }

    /**
     * 10. Sumatorias Clásicas O(1) con Módulo
     */
    ll sum_n(ll n, ll m = MOD) {
        n %= m;
        return n * (n + 1) % m * modinv(2, m) % m;
    }

    ll sum_squares(ll n, ll m = MOD) {
        n %= m;
        ll res = n * (n + 1) % m;
        res = res * (2 * n + 1) % m;
        return res * modinv(6, m) % m;
    }

    ll geom_sum(ll a, ll r, ll n, ll m = MOD) {
        if (r == 1) return (n % m) * (a % m) % m;
        ll num = (modpow(r, n, m) - 1 + m) % m;
        ll den = modinv(r - 1 + m, m);
        return (a % m) * num % m * den % m;
    }

    /**
     * 11. Matrices
     */
    struct Matrix {
        vector<vector<ll>> mat;
        int r, c;
        Matrix(int r, int c) : r(r), c(c), mat(r, vector<ll>(c, 0)) {}
        
        static Matrix identity(int n) {
            Matrix res(n, n);
            for (int i = 0; i < n; i++) res.mat[i][i] = 1;
            return res;
        }
        
        Matrix operator*(const Matrix &other) const {
            Matrix res(r, other.c);
            for (int i = 0; i < r; i++)
                for (int k = 0; k < c; k++)
                    for (int j = 0; j < other.c; j++)
                        res.mat[i][j] = (res.mat[i][j] + mat[i][k] * other.mat[k][j]) % MOD;
            return res;
        }
        
        Matrix power(ll p) const {
            Matrix res = identity(r);
            Matrix base = *this;
            while (p > 0) {
                if (p & 1) res = res * base;
                base = base * base;
                p >>= 1;
            }
            return res;
        }
    };

    vector<ll> berlekamp_massey(vector<ll> s) {
        vector<ll> C = {1}, B = {1};
        int L = 0; 
        ll m = 1, b = 1;
        for (int i = 0; i < s.size(); i++) {
            ll d = 0;
            for (int j = 0; j <= L; j++) d = (d + C[j] * s[i - j]) % MOD;
            if (d == 0) {
                m++;
            } else {
                vector<ll> T = C;
                ll c = MOD - d * modpow(b, MOD - 2) % MOD;
                while (C.size() <= B.size() + m) C.push_back(0);
                for (int j = 0; j < B.size(); j++) {
                    C[j + m] = (C[j + m] + c * B[j]) % MOD;
                }
                if (2 * L <= i) {
                    L = i + 1 - L;
                    B = T;
                    b = d;
                    m = 1;
                } else {
                    m++;
                }
            }
        }
        C.resize(L + 1);
        return C; // Retorna coeficientes C[0]*S[n] + C[1]*S[n-1] ... = 0
    }

    struct PollardRho {
        using u64 = uint64_t;
        using u128 = __uint128_t;

        static u64 modpow(u64 base, u64 exp, u64 mod) {
            u64 res = 1;
            base %= mod;
            while (exp > 0) {
                if (exp % 2 == 1) res = (u128)res * base % mod;
                base = (u128)base * base % mod;
                exp /= 2;
            }
            return res;
        }

        // Test de primalidad determinista para N <= 2^64
        static bool is_prime(u64 n) {
            if (n < 2) return false;
            if (n == 2 || n == 3) return true;
            if (n % 2 == 0) return false;
            u64 d = n - 1;
            int s = 0;
            while (d % 2 == 0) { d /= 2; s++; }
            static const u64 bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
            for (u64 a : bases) {
                if (n <= a) break;
                u64 x = modpow(a, d, n);
                if (x == 1 || x == n - 1) continue;
                bool composite = true;
                for (int r = 1; r < s; r++) {
                    x = (u128)x * x % n;
                    if (x == n - 1) { composite = false; break; }
                }
                if (composite) return false;
            }
            return true;
        }

        static u64 get_factor(u64 n) {
            if (n % 2 == 0) return 2;
            if (is_prime(n)) return n;
            u64 x = 2, y = 2, d = 1, c = 1;
            auto f = [&](u64 x, u64 n, u64 c) { return (u64)(((u128)x * x % n + c) % n); };
            while (d == 1) {
                x = f(x, n, c);
                y = f(f(y, n, c), n, c);
                d = gcd(x > y ? x - y : y - x, n);
                if (d == n) { 
                    x = rand() % (n - 2) + 2; 
                    y = x; 
                    c = rand() % (n - 1) + 1; 
                    d = 1; 
                }
            }
            return d;
        }

        static void _factorize(u64 n, map<u64, int>& factors) {
            if (n == 1) return;
            if (is_prime(n)) { factors[n]++; return; }
            u64 divisor = get_factor(n);
            _factorize(divisor, factors);
            _factorize(n / divisor, factors);
        }

        // API Principal
        static map<u64, int> factorize(u64 n) {
            map<u64, int> factors;
            _factorize(n, factors);
            return factors;
        }
    };

    /**
     * 12. Teorema de Lucas (nCr mod P para N, R gigantes y P primo pequeño <= 10^6)
     * Reutiliza la estructura Combinatorics para el cálculo de los dígitos en base P.
     */
    struct Lucas {
        Combinatorics comb;
        ll P;
        
        // Inicializa la combinatoria hasta P-1 usando el módulo P
        Lucas(ll p) : P(p), comb(p - 1, p) {}
        
        // Teorema de Lucas recursivo para N, R hasta 10^18
        ll solve(ll n, ll r) {
            if (r < 0 || r > n) return 0;
            if (r == 0) return 1;
            
            ll ni = n % P;
            ll ri = r % P;
            
            // Si el dígito inferior es mayor, el combinatorio se hace 0
            if (ni < ri) return 0; 
            
            // Llamada recursiva multiplicada por el nCr del dígito actual
            return solve(n / P, r / P) * comb.nCr(ni, ri) % P;
        }
    };

    /*
    =============================================
    CHEAT SHEET FÓRMULAS RÁPIDAS (EXPLICADAS)
    =============================================

    1. STARS AND BARS (Estrellas y Barras)
       ¿Para qué sirve? Calcular cuántas formas hay de repartir N objetos idénticos entre K cajas/personas.
       - Problema Típico: "Encuentra el número de soluciones enteras no negativas a: x1 + x2 + ... + xK = N"
       - Fórmula (Cajas pueden estar vacías): nCr(N + K - 1, K - 1)
       - Fórmula (Cada caja debe tener al menos 1): nCr(N - 1, K - 1)

    2. NÚMEROS DE CATALAN (C_n)
       ¿Para qué sirve? Contar estructuras recursivas que no se cruzan o se balancean.
       - Problema Típico: "¿Cuántas secuencias válidas de N pares de paréntesis existen?
         ¿Cuántos árboles binarios distintos se pueden hacer con N nodos? ¿De cuántas formas puedo
         triangular un polígono convexo de N+2 lados?"
       - Fórmula: nCr(2n, n) / (n + 1)

    3. TEOREMA DE PICK (Geometría en cuadrícula)
       ¿Para qué sirve? Relacionar el Área de un polígono con los puntos (coordenadas enteras) que toca.
       - Problema Típico: Te dan vértices de un campo en (X, Y) enteros. "¿Cuántos árboles (puntos enteros)
         quedan estrictamente encerrados dentro del campo?"
       - Fórmula: Area = Interior + (Borde / 2) - 1  --->  Interior = Area - (Borde / 2) + 1
       *Nota:* Los puntos de borde entre dos vértices (x1, y1) y (x2, y2) son: gcd(|x1 - x2|, |y1 - y2|).

    4. FÓRMULA DE CAYLEY (Grafos)
       ¿Para qué sirve? Contar árboles recubridores (Spanning Trees).
       - Problema Típico: "Tienes N ciudades completamente conectadas entre sí (grafo completo K_n).
         ¿De cuántas formas puedes elegir los caminos para que todas queden conectadas sin ciclos?"
       - Fórmula: N^(N - 2)

    5. TEOREMA DE EULER PARA GRAFOS PLANOS
       ¿Para qué sirve? Resolver problemas de regiones creadas por líneas o aristas que no se cruzan.
       - Problema Típico: "Dibujas N nodos y los conectas con E líneas sin que se crucen. ¿En cuántas
         'islas' o regiones queda dividido el mapa (incluyendo el exterior)?"
       - Fórmula: Caras = 2 - Vértices + Aristas   (F = 2 - V + E)

    6. DERANGEMENTS (Desarreglos - D_n)
       ¿Para qué sirve? Contar permutaciones donde NINGÚN elemento queda en su posición original.
       - Problema Típico: "N personas dejan su sombrero. Se los devuelven al azar. ¿De cuántas formas
         es posible que NINGUNA persona reciba su propio sombrero?"
       - Fórmula Recursiva: D_n = (n - 1) * (D_{n-1} + D_{n-2})  (Con D_0 = 1, D_1 = 0).
    =============================================
    */
}

int main() {
    // 1. Optimización rápida
    ios_base::sync_with_stdio(0); cin.tie(0);

    // ==========================================
    // EJEMPLOS DE USO MÁGICOS
    // ==========================================

    // 1. Inverso Modular y extGCD
    // ll x, y;
    // ll g = MathAlgo::extGCD(15, 26, x, y);
    // cout << "Inverso de 15 mod 26 es: " << MathAlgo::modinv(15, 26) << "\n";

    // 2. Teorema Chino del Resto (CRT)
    // auto [ans, lcm_val] = MathAlgo::crt(2, 3, 3, 5); // x=2(mod 3), x=3(mod 5)
    // cout << "x = " << ans << " mod " << lcm_val << "\n";

    // 3. Criba Lineal (Factorización ultra rápida)
    // MathAlgo::Sieve sieve(1e6);
    // vector<int> factors = sieve.factorize(100); // Devuelve {2, 2, 5, 5}

    // 4. Combinatoria Inteligente
    // MathAlgo::Combinatorics comb(100000);
    // cout << "10 en 3: " << comb.nCr(10, 3) << "\n";

    // 5. FFT - Multiplicación de Polinomios
    // // a(x) = 1 + 2x + 3x^2
    // // b(x) = 4 + 5x + 6x^2
    // vector<int> a = {1, 2, 3};
    // vector<int> b = {4, 5, 6};
    // vector<int> c = MathAlgo::multiply(a, b);
    // // c = {4, 13, 28, 27, 18} -> c(x) = 4 + 13x + 28x^2 + 27x^3 + 18x^4

    // 6. Gauss, Sistema a resolver:
    // 2x +  y = 5
    //  x -  y = 1
    //
    // Matriz extendida:
    // [2,  1, | 5]
    // [1, -1, | 1]

    // vector<vector<double>> matrix = {
    //     {2.0,  1.0, 5.0},
    //     {1.0, -1.0, 1.0}
    // };
    // vector<double> variables;
    // int solutions = MathAlgo::gauss(matrix, variables);
    //
    // if(solutions == 1) {
    //     // variables[0] será x (2.0)
    //     // variables[1] será y (1.0)
    //     cout << "x = " << variables[0] << ", y = " << variables[1] << "\n";
    // } else if (solutions == 0) {
    //     cout << "No hay solucion\n";
    // } else {
    //     cout << "Infinitas soluciones\n";
    // }

    // 7. Pollard Rho - Factorización de números grandes
    // PollardRho::u64 N = 1000000000000000003ULL * 2ULL; // Un número gigante
    // Test rápido de primalidad
    // if (PollardRho::is_prime(N)) {
        // cout << N << " es primo.\n";
    // } else {
        // Uso directo de la API
        // map<PollardRho::u64, int> factores = PollardRho::factorize(N);
        
        // cout << "Factores de " << N << ":\n";
        // for (auto [primo, potencia] : factores) {
            // cout << primo << "^" << potencia << "\n";
        // }
    // }

    // 8. Teorema de Lucas para combinaciones gigantes
    // Ejemplo: 10^18 en 5*10^17 módulo 1000003 (que es primo)
    // MathAlgo::Lucas lucas_solver(1000003);
    // cout << "Lucas (10^18 en 5*10^17) mod 1000003: " 
    //      << lucas_solver.solve(1000000000000000000LL, 500000000000000000LL) << "\n";

    return 0;
}
