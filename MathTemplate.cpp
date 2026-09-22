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
    struct Sieve {
        vector<int> primes, spf;
        Sieve(int n) {
            spf.assign(n + 1, 0);
            for (int i = 2; i <= n; i++) {
                if (spf[i] == 0) { spf[i] = i; primes.push_back(i); }
                for (int j = 0; j < primes.size() && primes[j] <= spf[i] && i * primes[j] <= n; j++) {
                    spf[i * primes[j]] = primes[j];
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

    return 0;
}
