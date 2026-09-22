#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")

#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace std;
using namespace __gnu_pbds;
using ll = long long;

// ============================================================================
// 1. PBDS (Policy-Based Data Structures)
// ============================================================================
namespace PBDS {
    // Set ordenado: elementos únicos, soporta consultas por índice en O(log N)
    template <typename T>
    using ordered_set = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;

    // Multiset ordenado: soporta duplicados guardando pares {valor, ID único}
    template <typename T>
    using ordered_multiset = tree<pair<T, int>, null_type, less<pair<T, int>>, rb_tree_tag, tree_order_statistics_node_update>;
}

// ============================================================================
// 2. Meet in the Middle (MitM)
// ============================================================================
struct MeetInTheMiddle {
    // Retorna true si existe un subconjunto de 'a' que sume exactamente 'target'
    // Complejidad: O(2^(N/2) * log(2^(N/2))), ideal para N <= 40
    static bool subset_sum(const vector<ll>& a, ll target) {
        int n = a.size();
        int mid = n / 2;
        vector<ll> left_half, right_half;

        // Fuerza bruta en la primera mitad
        for (int i = 0; i < (1 << mid); i++) {
            ll sum = 0;
            for (int j = 0; j < mid; j++) {
                if ((i >> j) & 1) sum += a[j];
            }
            left_half.push_back(sum);
        }

        // Fuerza bruta en la segunda mitad
        int n2 = n - mid;
        for (int i = 0; i < (1 << n2); i++) {
            ll sum = 0;
            for (int j = 0; j < n2; j++) {
                if ((i >> j) & 1) sum += a[mid + j];
            }
            right_half.push_back(sum);
        }

        // Ordenamos una de las mitades para poder hacer búsqueda binaria
        sort(right_half.begin(), right_half.end());

        // Buscamos si el complemento exacto existe en la segunda mitad
        for (ll sum1 : left_half) {
            if (binary_search(right_half.begin(), right_half.end(), target - sum1)) {
                return true;
            }
        }
        return false;
    }
};

// ============================================================================
// 3. Zobrist Hashing (XOR Hashing)
// ============================================================================
struct ZobristHash {
    using u64 = uint64_t;
    mt19937_64 rng;
    unordered_map<int, u64> val_hash;

    // Inicializamos el RNG con una semilla
    ZobristHash(u64 seed = 1337) : rng(seed) {
        val_hash[0] = 0; // El 0 puede representar un comodín o celda vacía
    }

    // Retorna el hash aleatorio único para el valor x (lo crea si no existe)
    u64 get_hash(int x) {
        if (!val_hash.count(x)) val_hash[x] = rng();
        return val_hash[x];
    }

    // Construye un arreglo de prefijos XOR en O(N)
    vector<u64> build_pref(const vector<int>& a) {
        int n = a.size();
        vector<u64> pref(n + 1, 0);
        for (int i = 0; i < n; i++) {
            pref[i + 1] = pref[i] ^ get_hash(a[i]);
        }
        return pref;
    }

    // Consulta el hash de un subsegmento 0-indexed [l, r] en O(1)
    u64 query_range(const vector<u64>& pref, int l, int r) {
        return pref[r + 1] ^ pref[l];
    }
};

// ============================================================================
// 4. Bitset Optimization (Cierre Transitivo / Reachability)
// ============================================================================
template<size_t MAXN>
struct BitsetReachability {
    bitset<MAXN> reachable[MAXN]; // Matriz de adyacencia de bits

    // Floyd-Warshall optimizado: O(V^3 / 64) en lugar de O(V^3)
    void build(int n, const vector<pair<int, int>>& edges) {
        // Inicializar: cada nodo se alcanza a sí mismo
        for (int i = 1; i <= n; i++) {
            reachable[i].reset();
            reachable[i][i] = 1; 
        }
        // Marcar aristas directas dirigidas
        for (auto& edge : edges) {
            reachable[edge.first][edge.second] = 1;
        }
        // Algoritmo principal (Nodos puente en el bucle exterior)
        for (int k = 1; k <= n; k++) {
            for (int i = 1; i <= n; i++) {
                if (reachable[i][k]) {
                    // Si 'i' llega a 'k', entonces 'i' puede llegar a todo lo que 'k' llega.
                    // El operador |= procesa 64 bits en 1 ciclo de reloj del procesador.
                    reachable[i] |= reachable[k]; 
                }
            }
        }
    }
    
    // Consulta en O(1) si hay camino de u a v
    bool can_reach(int u, int v) {
        return reachable[u][v];
    }
};

// ============================================================================
// Ejemplos de uso (Main)
// ============================================================================
int main() {
    ios_base::sync_with_stdio(0); cin.tie(0);

    // ---------------------------------------------------------
    // Ejemplo de uso: PBDS (Ordered Set / Multiset)
    // ---------------------------------------------------------
    PBDS::ordered_set<int> os;
    os.insert(10); os.insert(20); os.insert(5);
    
    // find_by_order(k) retorna un iterador al k-ésimo elemento (0-indexed)
    cout << "El elemento mas pequeno es: " << *os.find_by_order(0) << "\n"; // 5
    
    // order_of_key(x) retorna cuántos elementos son estrictamente menores a x
    cout << "Menores a 15: " << os.order_of_key(15) << "\n"; // 1 (solo el 5)

    // Para usar multiset (soportar duplicados):
    PBDS::ordered_multiset<int> oms;
    int id = 0; 
    oms.insert({5, id++}); oms.insert({5, id++}); // Insertamos dos '5' distintos
    cout << "Menores a 6 en multiset: " << oms.order_of_key({6, -1}) << "\n"; // 2


    // ---------------------------------------------------------
    // Ejemplo de uso: Meet in the Middle
    // ---------------------------------------------------------
    vector<ll> weights = {10, 20, 15, 5, 30};
    ll target = 35;
    
    // Verifica si algún subconjunto suma exactamente 35 en O(2^(N/2) * log(2^(N/2)))
    if (MeetInTheMiddle::subset_sum(weights, target)) {
        cout << "Existe un subconjunto que suma " << target << "\n";
    }


    // ---------------------------------------------------------
    // Ejemplo de uso: Zobrist Hashing (XOR Hashing)
    // ---------------------------------------------------------
    ZobristHash zh;
    vector<int> A = {1, 2, 3, 2, 1};
    vector<int> B = {3, 1, 1, 2, 2}; // Mismos elementos, orden distinto

    // Construimos los arreglos de prefijos
    vector<ZobristHash::u64> prefA = zh.build_pref(A);
    vector<ZobristHash::u64> prefB = zh.build_pref(B);

    // Comparamos si A[0...4] es un anagrama de B[0...4]
    if (zh.query_range(prefA, 0, 4) == zh.query_range(prefB, 0, 4)) {
        cout << "El arreglo B es una permutacion (anagrama) exacta de A.\n";
    }


    // ---------------------------------------------------------
    // Ejemplo de uso: Bitset Optimization
    // ---------------------------------------------------------
    int N = 5;
    vector<pair<int, int>> aristas = {{1, 2}, {2, 3}, {4, 5}};
    
    // Instanciamos el struct. MAXN debe ser una constante conocida en tiempo de compilación.
    BitsetReachability<2500> reach;
    reach.build(N, aristas); // Calculamos en O(N^3 / 64)

    // Consultamos caminos en O(1)
    if (reach.can_reach(1, 3)) cout << "Hay un camino dirigido del 1 al 3.\n";
    if (!reach.can_reach(1, 5)) cout << "No hay camino del 1 al 5.\n";

    return 0;
}