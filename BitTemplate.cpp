/**
 * @file BitTemplateICPC.cpp
 * @brief Plantilla de algoritmos a nivel de bits y máscaras.
 * @details Incluye Base de XOR (GF(2)), Trie de Bits (Max XOR), Gosper's Hack y utilidades O(1).
 * @note En problemas de 64 bits (ll) recuerda usar '1LL << i', NO '1 << i'.
 */
//  ____ ___ _____  _    _     ____  ___ 
// | __ )_ _|_   _|/ \  | |   / ___|/ _ \
// |  _ \| |  | | / _ \ | |  | |  _| | | |
// | |_) | |  | |/ ___ \| |__| |_| | |_| |
// |____/___| |_/_/   \_\____|\____|\___/
//
//           BIT TEMPLATE

#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

// ==========================================
// ESTRUCTURA MAESTRA DE BITS (ICPC)
// ==========================================

namespace BitAlgo {

    /**
     * 1. Macros y Utilidades O(1) (Built-in GCC)
     */
    // Extraer el bit i-ésimo (0-indexed)
    inline bool check_bit(ll mask, int i) { return (mask >> i) & 1; }
    // Encender el bit i-ésimo
    inline ll set_bit(ll mask, int i) { return mask | (1LL << i); }
    // Apagar el bit i-ésimo
    inline ll clear_bit(ll mask, int i) { return mask & ~(1LL << i); }
    // Cambiar (Toggle) el bit i-ésimo
    inline ll toggle_bit(ll mask, int i) { return mask ^ (1LL << i); }
    // Extraer el bit menos significativo encendido (LSB) ej: 10100 -> 00100
    inline ll lsb(ll mask) { return mask & -mask; }
    
    // Conteo de bits encendidos
    inline int count_bits(int mask) { return __builtin_popcount(mask); }
    inline int count_bits_ll(ll mask) { return __builtin_popcountll(mask); }
    
    // Contar ceros a la izquierda (Leading Zeros) - Útil para saber el bit más alto
    inline int clz(ll mask) { return mask == 0 ? 64 : __builtin_clzll(mask); }
    inline int highest_bit(ll mask) { return mask == 0 ? -1 : 63 - __builtin_clzll(mask); }
    
    // Contar ceros a la derecha (Trailing Zeros) - Útil para saber el bit más bajo
    inline int ctz(ll mask) { return mask == 0 ? 64 : __builtin_ctzll(mask); }


    /**
     * 2. Iterar sobre TODAS las sub-máscaras de una máscara dada.
     * Complejidad: O(3^N) si se hace para todas las máscaras de 0 a 2^N-1.
     */
    void iterate_submasks(ll mask) {
        // Itera pasando por todas las sub-máscaras válidas excluyendo el 0
        for (ll sub = mask; sub > 0; sub = (sub - 1) & mask) {
            // ... procesar sub-máscara ...
        }
        // Si necesitas procesar el 0, simplemente lo haces fuera del loop
    }


    /**
     * 3. Gosper's Hack
     * Genera todas las máscaras de exactamente 'K' bits encendidos de tamaño 'N'.
     * Complejidad: O(N en K) exacto. (Mucho más rápido que iterar 2^N y hacer popcount)
     */
    vector<ll> gospers_hack(int k, int n) {
        vector<ll> masks;
        if (k == 0) { masks.push_back(0); return masks; }
        ll mask = (1LL << k) - 1;
        ll limit = (1LL << n);
        
        while (mask < limit) {
            masks.push_back(mask);
            ll c = mask & -mask;
            ll r = mask + c;
            mask = (((r ^ mask) >> 2) / c) | r;
        }
        return masks;
    }


    /**
     * 4. Base de XOR (Álgebra Lineal sobre GF(2))
     * Sirve para: Maximizar XOR de un subconjunto, encontrar si un XOR objetivo es posible, 
     * o encontrar la cantidad de subconjuntos distintos que generan el mismo XOR.
     * Complejidad: O(BITS) por inserción/consulta. (Normalmente BITS = 60).
     */
    struct XorBasis {
        static const int BITS = 60;
        ll basis[BITS];
        int sz = 0; // Tamaño de la base (rango de la matriz)

        XorBasis() { memset(basis, 0, sizeof(basis)); }

        // Inserta un número a la base. Retorna true si aumentó el rango (era linealmente independiente)
        bool insert(ll mask) {
            for (int i = BITS - 1; i >= 0; i--) {
                if (!check_bit(mask, i)) continue;
                if (!basis[i]) {
                    basis[i] = mask;
                    sz++;
                    return true;
                }
                mask ^= basis[i];
            }
            return false;
        }

        // Determina si un número 'k' puede ser formado haciendo XOR de algún subconjunto
        bool can_form(ll k) {
            for (int i = BITS - 1; i >= 0; i--) {
                if (check_bit(k, i)) {
                    if (!basis[i]) return false;
                    k ^= basis[i];
                }
            }
            return k == 0;
        }

        // Devuelve el valor MÁXIMO posible haciendo XOR con un subconjunto de la base
        // Si se pasa 'initial', devuelve el máximo de (initial ^ subset)
        ll get_max(ll initial = 0) {
            ll res = initial;
            for (int i = BITS - 1; i >= 0; i--) {
                if ((res ^ basis[i]) > res) {
                    res ^= basis[i];
                }
            }
            return res;
        }
        
        // Devuelve el MÍNIMO valor posible mayor que 0 haciendo XOR con algún subconjunto
        ll get_min() {
            for (int i = 0; i < BITS; i++) {
                if (basis[i]) return basis[i];
            }
            return 0;
        }
    };


    /**
     * 5. Trie de Bits (Bit Trie)
     * Sirve para: Encontrar el máximo/mínimo XOR de un elemento 'X' contra un ARREGLO de elementos.
     * Complejidad: O(BITS) estricto por inserción/búsqueda.
     */
    struct BitTrie {
        static const int BITS = 30; // Usar 60 si los números llegan a 10^18
        struct Node {
            int next[2];
            int cnt; // Cuenta de elementos que pasan por el nodo (útil para borrar)
            Node() { next[0] = next[1] = -1; cnt = 0; }
        };

        vector<Node> t;
        
        BitTrie() { t.emplace_back(); }

        void insert(ll x, int val = 1) { // val = 1 para insertar, val = -1 para borrar
            int u = 0;
            for (int i = BITS - 1; i >= 0; i--) {
                int bit = (x >> i) & 1;
                if (t[u].next[bit] == -1) {
                    t[u].next[bit] = t.size();
                    t.emplace_back();
                }
                u = t[u].next[bit];
                t[u].cnt += val;
            }
        }

        // Devuelve el máximo valor de (x ^ y) donde 'y' está en el Trie
        ll max_xor(ll x) {
            int u = 0;
            ll ans = 0;
            for (int i = BITS - 1; i >= 0; i--) {
                int bit = (x >> i) & 1;
                int opposite = bit ^ 1;
                
                // Si el camino opuesto existe y tiene elementos, lo tomamos para maximizar el XOR
                if (t[u].next[opposite] != -1 && t[t[u].next[opposite]].cnt > 0) {
                    ans |= (1LL << i);
                    u = t[u].next[opposite];
                } else {
                    u = t[u].next[bit];
                }
            }
            return ans;
        }
    };
}

/*
=============================================
CHEAT SHEET PARA BITS
=============================================

1. TRUCOS CLÁSICOS DE MÁSCARAS:
   - x & (x - 1): Apaga el bit encendido más a la derecha (LSB). Útil para saber si 'x' es potencia de 2 (da 0).
   - x & -x: Aísla el LSB. Deja todos los ceros igual, deja el LSB igual y apaga todos los bits a la izquierda.
   - x | (x + 1): Enciende el bit apagado más a la derecha.
   - x | (x - 1): Enciende todos los bits a la derecha del LSB.
   - ~x & (x + 1): Aísla el cero más a la derecha (lo vuelve 1 y los demás 0).

2. XOR BASIS (Base de XOR):
   - Cualquier conjunto de números puede reducirse a una "base" de tamaño máximo 60 (para long long).
   - El tamaño de la base `sz` es el "rango" (Rank) de la matriz de bits.
   - Si un conjunto tiene tamaño N y su base tiene tamaño R, hay 2^(N - R) formas de formar 
     el XOR 0 usando subconjuntos (incluyendo el vacío).
   - Por extensión, si un número 'K' se puede formar, hay exactamente 2^(N - R) formas distintas 
     de formarlo con el arreglo original.

3. TRIE DE BITS vs XOR BASIS:
   - XOR Basis te da el máximo XOR operando VARIOS elementos (un subconjunto entero).
   - Bit Trie te da el máximo XOR operando UN SOLO elemento del arreglo contra tu query 'X'.

=============================================
*/


int main() {
    ios_base::sync_with_stdio(0); cin.tie(0);

    // ==========================================
    // EJEMPLOS DE USO MÁGICOS
    // ==========================================

    // 1. Funciones GCC O(1)
    // ll mask = 18; // Binario: 10010
    // cout << "Leading zeros de 18: " << BitAlgo::clz(mask) << "\n";
    // cout << "Bit mas significativo: " << BitAlgo::highest_bit(mask) << "\n";
    // cout << "Bits encendidos: " << BitAlgo::count_bits_ll(mask) << "\n";
    // cout << "Mascara sin el LSB: " << (mask ^ BitAlgo::lsb(mask)) << "\n";

    // 2. Iterar sobre sub-máscaras
    // ll m = 5; // Binario: 101. Sus submascaras validas: 101 (5), 100 (4), 001 (1)
    // // BitAlgo::iterate_submasks(m);

    // 3. Gosper's Hack (Todas las combinaciones N en K)
    // // Generar todas las máscaras de longitud 5 con exactamente 3 bits encendidos
    // vector<ll> combs = BitAlgo::gospers_hack(3, 5);
    // // combs contendrá: 00111, 01011, 01101, 01110, 10011, 10101, 10110, 11001, 11010, 11100

    // 4. Base de XOR (GF(2))
    // BitAlgo::XorBasis basis;
    // vector<ll> arr = {1, 2, 3}; 
    // for(ll x : arr) basis.insert(x); // '3' no aumentará el tamaño porque 1^2 = 3
    // 
    // cout << "Rango de la base: " << basis.sz << "\n";
    // cout << "El 3 se puede formar? " << (basis.can_form(3) ? "Si" : "No") << "\n";
    // cout << "Maximo XOR posible con subset de arr: " << basis.get_max() << "\n";

    // 5. Trie de Bits (Max XOR subarray/pair)
    // BitAlgo::BitTrie trie;
    // trie.insert(3); // 011
    // trie.insert(10); // 1010
    // trie.insert(5);  // 0101
    // 
    // // Queremos maximizar (6 ^ Y) donde Y está en el trie.
    // // 6 en binario es 0110. El óptimo será 1010 (10) -> 6 ^ 10 = 12 (1100).
    // cout << "Maximo XOR con 6: " << trie.max_xor(6) << "\n";

    return 0;
}