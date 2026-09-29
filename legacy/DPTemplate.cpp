/**
 * @file DPTemplateICPC.cpp
 * @brief Plantilla de Programación Dinámica (Optimizaciones y Clásicos).
 * @details Incluye Li Chao Tree (CHT dinámico), Divide & Conquer, Knuth, Digit DP, SOS DP, LIS y Bitset Knapsack.
 * @note Árboles (Rerooting) excluidos por estar en plantilla separada.
 */
//  ____  ____    _____ _____ __  __ ____  _        _  _____ _____ 
// |  _ \|  _ \  |_   _| ____|  \/  |  _ \| |      / \|_   _| ____|
// | | | | |_) |   | | |  _| | |\/| | |_) | |     / _ \ | | |  _|  
// | |_| |  __/    | | | |___| |  | |  __/| |___ / ___ \| | | |___ 
// |____/|_|       |_| |_____|_|  |_|_|   |_____/_/   \_\_| |_____|
//
//           DYNAMIC PROGRAMMING TEMPLATE

#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

const ll INF = 1e18; // Usar INF seguro para DP

// ==========================================
// ESTRUCTURA MAESTRA DE DP (ICPC)
// ==========================================

namespace DPAlgo {

    /**
     * 1. Li Chao Tree (Convex Hull Trick Dinámico) - O(log(MAX_X))
     * Sirve para: dp[i] = min/max(m_j * x_i + c_j)
     * Soporta queries y pendientes no monótonas en rangos grandes (ej. -1e9 a 1e9).
     * Configurado para MINIMIZAR (Cambiar '>' por '<' en la evaluación para maximizar).
     */
    struct Line {
        ll m, c;
        Line() : m(0), c(INF) {} // c = INF para Min, c = -INF para Max
        Line(ll m, ll c) : m(m), c(c) {}
        ll eval(ll x) const { return m * x + c; }
    };

    struct LiChaoTree {
        struct Node {
            Line line;
            int lc = -1, rc = -1;
        };
        vector<Node> tree;
        ll min_x, max_x;

        LiChaoTree(ll min_x = -1e9, ll max_x = 1e9) : min_x(min_x), max_x(max_x) {
            tree.emplace_back();
        }

        void add_line(Line new_line, int v, ll l, ll r) {
            ll m = l + (r - l) / 2;
            bool left_better = new_line.eval(l) < tree[v].line.eval(l);
            bool mid_better = new_line.eval(m) < tree[v].line.eval(m);

            if (mid_better) swap(tree[v].line, new_line);

            if (l == r) return;

            if (left_better != mid_better) {
                if (tree[v].lc == -1) { tree[v].lc = tree.size(); tree.emplace_back(); }
                add_line(new_line, tree[v].lc, l, m);
            } else {
                if (tree[v].rc == -1) { tree[v].rc = tree.size(); tree.emplace_back(); }
                add_line(new_line, tree[v].rc, m + 1, r);
            }
        }
        void add(ll m, ll c) { add_line(Line(m, c), 0, min_x, max_x); }

        ll query(ll x, int v, ll l, ll r) {
            if (v == -1) return INF; // Retornar -INF si se está maximizando
            ll m = l + (r - l) / 2;
            ll res = tree[v].line.eval(x);
            if (l == r) return res;
            if (x <= m) return min(res, query(x, tree[v].lc, l, m));
            else return min(res, query(x, tree[v].rc, m + 1, r));
        }
        ll query(ll x) { return query(x, 0, min_x, max_x); }
    };

    /**
     * 2. Divide and Conquer DP - O(K * N log N)
     * Sirve para: dp[i][j] = min_{k < j}(dp[i-1][k] + Cost(k+1, j))
     * REQUISITO: El punto óptimo k cumple optimal_k(i, j) <= optimal_k(i, j+1).
     * @param dp_prev Fila i-1 de la matriz DP.
     * @param dp_curr Fila i de la matriz DP (se llena en esta función).
     * @param cost Lambda/Función que retorna el costo entre [l, r] en O(1) o amortizado.
     */
    template<class CostFunc>
    void dnc(int l, int r, int optl, int optr, const vector<ll>& dp_prev, vector<ll>& dp_curr, CostFunc& cost) {
        if (l > r) return;
        int mid = l + (r - l) / 2;
        pair<ll, int> best = {INF, -1};

        for (int k = optl; k <= min(mid - 1, optr); k++) { // k es el punto de partición
            ll current_cost = (dp_prev[k] == INF) ? INF : dp_prev[k] + cost(k + 1, mid);
            if (current_cost < best.first) {
                best = {current_cost, k};
            }
        }
        
        dp_curr[mid] = best.first;
        int opt = best.second;
        if(opt == -1) opt = optl; // Manejo de casos imposibles

        dnc(l, mid - 1, optl, opt, dp_prev, dp_curr, cost);
        dnc(mid + 1, r, opt, optr, dp_prev, dp_curr, cost);
    }

    /**
     * 3. Knuth Optimization - O(N^2) en lugar de O(N^3)
     * Sirve para: dp[i][j] = min_{i < k < j}(dp[i][k] + dp[k][j]) + Cost(i, j)
     * REQUISITO: optimal_k(i, j-1) <= optimal_k(i, j) <= optimal_k(i+1, j)
     * @param N Número de elementos (0-indexed hasta N).
     * @param cost Lambda/Función que calcula el costo C(i, j).
     */
    template<class CostFunc>
    vector<vector<ll>> knuth(int n, CostFunc& cost) {
        vector<vector<ll>> dp(n + 1, vector<ll>(n + 1, 0));
        vector<vector<int>> opt(n + 1, vector<int>(n + 1, 0));

        // Inicialización de intervalos de longitud 1 y 2
        for (int i = 0; i < n; i++) {
            opt[i][i + 1] = i; // Punto óptimo base
        }

        for (int len = 2; len <= n; len++) {
            for (int i = 0; i + len <= n; i++) {
                int j = i + len;
                dp[i][j] = INF;
                for (int k = opt[i][j - 1]; k <= min(j - 1, opt[i + 1][j]); k++) {
                    ll current = dp[i][k] + dp[k][j] + cost(i, j);
                    if (current < dp[i][j]) {
                        dp[i][j] = current;
                        opt[i][j] = k;
                    }
                }
            }
        }
        return dp;
    }

    /**
     * 4. Esqueleto de Digit DP
     * Sirve para: "Cuántos números entre [L, R] cumplen X condición".
     * Se llama como: solve(R) - solve(L-1).
     * Modificar `memo` y estado interno según el problema (suma, ceros, etc).
     */
    namespace DigitDP {
        string s;
        ll memo[20][2][2]; // [índice][tight][leading_zero] ... agregar estados.

        ll dp(int idx, bool tight, bool lz) {
            if (idx == s.size()) {
                return 1; // Retorna 1 si el número formado es válido, sino 0
            }
            if (memo[idx][tight][lz] != -1) return memo[idx][tight][lz];

            int limit = tight ? (s[idx] - '0') : 9;
            ll ans = 0;

            for (int d = 0; d <= limit; d++) {
                bool next_tight = tight && (d == limit);
                bool next_lz = lz && (d == 0);
                
                // Transición: ans += dp(idx + 1, next_tight, next_lz, ...)
                ans += dp(idx + 1, next_tight, next_lz);
            }
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
     * 5. Sum Over Subsets (SOS DP / FMT) - O(N * 2^N)
     * Calcula para toda máscara M: res[M] = Suma de todos sus subconjuntos.
     * Ejemplo: dp[mask] = sum(dp[submask]) donde (submask & mask) == submask.
     */
    void sos_dp(vector<ll>& dp, int n) {
        // n es la cantidad de bits. El arreglo 'dp' debe ser de tamaño 2^n.
        for (int i = 0; i < n; i++) {
            for (int mask = 0; mask < (1 << n); mask++) {
                if (mask & (1 << i)) {
                    dp[mask] += dp[mask ^ (1 << i)];
                }
            }
        }
    }

    /**
     * 6. LIS (Longest Increasing Subsequence) + Reconstrucción - O(N log N)
     * Retorna {Longitud, Secuencia Reconstruida}
     */
    pair<int, vector<int>> lis(const vector<int>& a) {
        int n = a.size();
        vector<int> d, pos; // d = valores mínimos, pos = índices en 'a'
        vector<int> parent(n, -1); // Para reconstruir el camino

        for (int i = 0; i < n; i++) {
            auto it = lower_bound(d.begin(), d.end(), a[i]); // Usar upper_bound para No Decreciente
            int idx = it - d.begin();
            
            if (it == d.end()) {
                d.push_back(a[i]);
                pos.push_back(i);
            } else {
                *it = a[i];
                pos[idx] = i;
            }
            if (idx > 0) parent[i] = pos[idx - 1];
        }

        vector<int> seq;
        int curr = pos.empty() ? -1 : pos.back();
        while (curr != -1) {
            seq.push_back(a[curr]);
            curr = parent[curr];
        }
        reverse(seq.begin(), seq.end());
        return {d.size(), seq};
    }

    /**
     * 7. Bitset Knapsack - O((N * MAX_SUM) / 64)
     * Determina si es posible formar exactamente una 'suma' usando un subconjunto.
     */
    const int MAX_SUM = 1e5 + 5; 
    bool bitset_knapsack(const vector<int>& weights, int target) {
        bitset<MAX_SUM> dp;
        dp[0] = 1;
        for (int w : weights) {
            dp |= (dp << w); // Transición ultra rápida a nivel de bits
        }
        return dp[target];
    }
}

/*
=============================================
CHEAT SHEET PARA RECONOCER OPTIMIZACIONES DP
=============================================

1. Convex Hull Trick / Li Chao Tree
   Ecuación: dp[i] = min( dp[j] + m[j] * x[i] + c[j] )
   Pista: Multiplicación de dos variables, una depende exclusívamente de 'i' y otra de 'j'.

2. Divide & Conquer Optimization
   Ecuación: dp[i][j] = min( dp[i-1][k] + Cost(k+1, j) )
   Pista: Agrupar un arreglo 1D en K particiones continuas. La función de costo (ej. suma de pares, 
   elementos únicos) cumple con la "Desigualdad Quadrangle" (si extiendes el rango, el óptimo no retrocede).

3. Knuth Optimization
   Ecuación: dp[i][j] = min( dp[i][k] + dp[k][j] ) + Cost(i, j)
   Pista: Problemas de fusionar intervalos adyacentes hasta tener un solo intervalo grande. 
   (Ej: Pegar trozos de cuerda, árboles de búsqueda óptimos).

4. SOS DP (Sum Over Subsets)
   Pista: Te piden calcular iteraciones sobre combinaciones binarias. Ejemplo: "Dado un arreglo A, 
   calcula para cada máscara X la suma de A[Y] donde Y es subconjunto de X". (Y & X == Y).

=============================================
*/


int main() {
    ios_base::sync_with_stdio(0); cin.tie(0);

    // ==========================================
    // EJEMPLOS DE USO MÁGICOS
    // ==========================================

    // 1. Li Chao Tree (CHT)
    // Minimizar y = mx + c para un 'x' dado
    // DPAlgo::LiChaoTree cht;
    // cht.add(-2, 10); // Línea: y = -2x + 10
    // cht.add(1, 2);   // Línea: y = 1x + 2
    // cout << "Minimo en x=3: " << cht.query(3) << "\n"; // Evalúa ambas rectas

    // 2. Divide and Conquer DP
    // Ejemplo: Dividir N elementos en K grupos
    // vector<ll> dp_prev(N + 1, INF), dp_curr(N + 1, 0);
    // dp_prev[0] = 0;
    // auto cost_func = [&](int l, int r) { return cost_matrix[l][r]; };
    // for(int k = 1; k <= K; k++) {
    //     DPAlgo::dnc(1, N, 1, N, dp_prev, dp_curr, cost_func);
    //     dp_prev = dp_curr;
    // }

    // 3. Knuth Optimization
    // Ejemplo: Optimal Binary Search Tree / Slimes / Matrix Chain
    // int N = 4;
    // auto cost_f = [&](int i, int j) { return prefix_sum[j] - prefix_sum[i]; };
    // vector<vector<ll>> dp = DPAlgo::knuth(N, cost_f);
    // cout << "Costo optimo: " << dp[0][N] << "\n";

    // 4. Digit DP
    // cout << "Numeros validos entre 10 y 100: " 
    //      << DPAlgo::DigitDP::solve(100) - DPAlgo::DigitDP::solve(9) << "\n";

    // 5. SOS DP
    // int bits = 3;
    // vector<ll> dp = {1, 2, 3, 4, 5, 6, 7, 8}; // Tamaño 2^3 = 8
    // DPAlgo::sos_dp(dp, bits);
    // // Ahora dp[mask] tiene la suma de sus submasks

    // 6. LIS con reconstrucción
    // vector<int> a = {10, 22, 9, 33, 21, 50, 41, 60};
    // auto [len, seq] = DPAlgo::lis(a);
    // cout << "LIS Len: " << len << " | Secuencia: ";
    // for(int x : seq) cout << x << " "; 

    // 7. Bitset Knapsack
    // vector<int> pesos = {2, 5, 8};
    // cout << (DPAlgo::bitset_knapsack(pesos, 13) ? "Posible" : "Imposible") << "\n";

    return 0;
}