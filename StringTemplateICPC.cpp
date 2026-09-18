/**
 * @file StringTemplateICPC.cpp
 * @brief Plantilla de algoritmos para cadenas y patrones.
 * @details Incluye KMP, Z-function, KMP Stack, hashing doble, Trie, Aho-Corasick, Manacher y Suffix Automaton.
 * @note Verifica el alfabeto y la base del hash segun las restricciones.
 */
//   ____ ___  ____  _____   ____  _   _
//  / ___/ _ \|  _ \| ____| / ___|| | | |
// | |  | | | | | | |  _|   \___ \| | | |
// | |__| |_| | |_| | |___   ___) | |_| |
//  \____\___/|____/|_____| |____/ \___/
//
//                    STRING TEMPLATE

#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

const int INF = 1e9;

// ==========================================
// ESTRUCTURA MAESTRA DE STRINGS (ICPC)
// ==========================================

namespace StringAlgo {

    /**
     * 1. KMP (Prefix Function / Pi Array) - O(N)
     * pi[i] = longitud del prefijo propio más largo que también es sufijo en s[0..i]
     */
    vector<int> prefix_function(const string& s) {
        int n = s.size();
        vector<int> pi(n, 0);
        for (int i = 1; i < n; i++) {
            int j = pi[i - 1];
            while (j > 0 && s[i] != s[j]) j = pi[j - 1];
            if (s[i] == s[j]) j++;
            pi[i] = j;
        }
        return pi;
    }

    /**
     * 2. Algoritmo Z (Z-Array) - O(N)
     * z[i] = longitud del prefijo más largo de s que coincide con el prefijo de s[i..N-1]
     */
    vector<int> z_function(const string& s) {
        int n = s.size(), l = 0, r = 0;
        vector<int> z(n, 0);
        for (int i = 1; i < n; i++) {
            if (i <= r) z[i] = min(r - i + 1, z[i - l]);
            while (i + z[i] < n && s[z[i]] == s[i + z[i]]) z[i]++;
            if (i + z[i] - 1 > r) { l = i; r = i + z[i] - 1; }
        }
        return z;
    }

    /**
     * 3. KMP Stack (Borrar todas las ocurrencias en O(N))
     * Elimina iterativamente todas las ocurrencias del patrón P en el texto S.
     * Retorna el string resultante limpio.
     */
    string remove_all_occurrences(const string& s, const string& p) {
        vector<int> pi_p = prefix_function(p);
        string ans = "";
        vector<int> estado(s.size() + 1, 0);

        for (char c : s) {
            ans.push_back(c);
            int j = estado[ans.size() - 1];

            while (j > 0 && p[j] != c) {
                j = pi_p[j - 1];
            }
            if (p[j] == c) j++;

            estado[ans.size()] = j;

            if (j == p.size()) {
                ans.resize(ans.size() - p.size());
            }
        }
        return ans;
    }

    /**
     * 4. Algoritmo de Manacher - O(N)
     * Encuentra la longitud del palíndromo más largo centrado en cada posición.
     * Intercala '#' para manejar palíndromos de longitud par e impar uniformemente.
     * Retorna p, donde p[i] es la longitud del palíndromo en el string original.
     */
    vector<int> manacher(const string& s) {
        string t = "^#";
        for (char c : s) {
            t += c;
            t += '#';
        }
        t += '$'; // Delimitadores para evitar checar desbordamientos de índices

        int n = t.size();
        vector<int> p(n, 0);
        int c = 0, r = 0; // Centro y límite derecho del palíndromo más a la derecha

        for (int i = 1; i < n - 1; i++) {
            // Espejo de i respecto al centro c
            int i_mirror = 2 * c - i;

            if (r > i) {
                p[i] = min(r - i, p[i_mirror]);
            } else {
                p[i] = 0;
            }

            // Expandir el palíndromo centrado en i
            while (t[i + 1 + p[i]] == t[i - 1 - p[i]]) {
                p[i]++;
            }

            // Si el palíndromo centrado en i supera r, ajustar c y r
            if (i + p[i] > r) {
                c = i;
                r = i + p[i];
            }
        }

        // Devolvemos solo la parte útil, quitando el delimitador inicial '^' y el final '$'
        return vector<int>(p.begin() + 1, p.end() - 1);
    }
}

/**
 * 5. Double Rolling Hash (O(N) build, O(1) query)
 * Ultra robusto contra "Anti-Hash tests" en Codeforces.
 */
struct StringHash {
    string s;
    int n;
    vector<ll> h1, h2, p1, p2;
    const ll MOD1 = 1e9 + 7, MOD2 = 1e9 + 9;
    const ll BASE1 = 31, BASE2 = 37;

    StringHash(const string& _s) : s(_s), n(s.size()) {
        h1.assign(n + 1, 0); h2.assign(n + 1, 0);
        p1.assign(n + 1, 1); p2.assign(n + 1, 1);
        for (int i = 0; i < n; i++) {
            h1[i + 1] = (h1[i] * BASE1 + (s[i] - 'a' + 1)) % MOD1;
            h2[i + 1] = (h2[i] * BASE2 + (s[i] - 'a' + 1)) % MOD2;
            p1[i + 1] = (p1[i] * BASE1) % MOD1;
            p2[i + 1] = (p2[i] * BASE2) % MOD2;
        }
    }

    // Devuelve el Hash del substring S[l..r] (0-indexed, inclusivo) en O(1)
    // Empaquetado en un solo uint64_t (Mitad superior: hash1, Mitad inferior: hash2)
    uint64_t get_hash(int l, int r) {
        ll hash1 = (h1[r + 1] - (h1[l] * p1[r - l + 1]) % MOD1 + MOD1) % MOD1;
        ll hash2 = (h2[r + 1] - (h2[l] * p2[r - l + 1]) % MOD2 + MOD2) % MOD2;

        // Desplazamos hash1 a los 32 bits más significativos y lo unimos con hash2
        return ((uint64_t)hash1 << 32) | (uint32_t)hash2;
    }
};

/**
 * 6. Trie (Árbol de Prefijos)
 * Útil para diccionarios y conteo de prefijos. Usa vectores dinámicos (evita punteros lentos).
 */
struct Trie {
    struct Node {
        int next[26];
        bool is_end = false;
        int count = 0; // Cuántas palabras pasan por este nodo
        Node() { fill(next, next + 26, -1); }
    };

    vector<Node> t;

    Trie() { t.emplace_back(); }

    void insert(const string& s) {
        int v = 0;
        for (char ch : s) {
            int c = ch - 'a';
            if (t[v].next[c] == -1) {
                t[v].next[c] = t.size();
                t.emplace_back();
            }
            v = t[v].next[c];
            t[v].count++;
        }
        t[v].is_end = true;
    }

    bool search(const string& s) {
        int v = 0;
        for (char ch : s) {
            int c = ch - 'a';
            if (t[v].next[c] == -1) return false;
            v = t[v].next[c];
        }
        return t[v].is_end;
    }
};

/**
 * 7. Aho-Corasick (Búsqueda Simultánea de Múltiples Patrones - 0 TLE)
 * Construcción: O(Sumatoria de longitudes de patrones)
 * Búsqueda: O(|S| + Sumatoria de longitudes) ESTRICTO.
 * Utiliza DP Lazy sobre el Suffix Tree. Soporta múltiples queries sobre textos distintos.
 */
struct AhoCorasick {

    struct Node {
        int next[26];
        int link = 0;         // Enlace de fallo (suffix link)
        vector<int> word_ids; // IDs de los patrones que terminan exactamente aquí
        Node() { fill(next, next + 26, -1); }
    };

    vector<Node> t;
    vector<int> bfs_order; // Orden topológico para la DP perezosa
    vector<int> pat_len;   // Longitud de cada patrón
    int num_patterns = 0;

    AhoCorasick() {
        t.emplace_back();
    }

    // Inserta un patrón y le asigna un ID automático (0, 1, 2...).
    // Retorna el ID asignado.
    int insert(const string& s) {
        int v = 0;
        int id = num_patterns++;
        pat_len.push_back(s.size());

        for (char ch : s) {
            int c = ch - 'a';
            if (t[v].next[c] == -1) {
                t[v].next[c] = t.size();
                t.emplace_back();
            }
            v = t[v].next[c];
        }
        t[v].word_ids.push_back(id);
        return id;
    }

    // OBLIGATORIO: Llamar una única vez después de insertar TODOS los patrones.
    void build() {
        queue<int> q;
        for (int c = 0; c < 26; c++) {
            if (t[0].next[c] != -1) {
                t[t[0].next[c]].link = 0;
                q.push(t[0].next[c]);
            } else {
                t[0].next[c] = 0;
            }
        }

        bfs_order.clear();
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            bfs_order.push_back(u);

            for (int c = 0; c < 26; c++) {
                int v = t[u].next[c];
                if (v != -1) {
                    t[v].link = t[t[u].link].next[c];
                    q.push(v);
                } else {
                    t[u].next[c] = t[t[u].link].next[c]; // Cerrar el autómata
                }
            }
        }
    }

    // 1. Devuelve un vector donde ans[id] es cuántas veces aparece el patrón 'id'
    vector<int> count_occurrences(const string& text) {
        vector<int> freq(t.size(), 0);
        int u = 0;
        for (char ch : text) {
            u = t[u].next[ch - 'a'];
            freq[u]++;
        }

        // DP Lazy: Propagar frecuencias de hojas a raíz
        for (int i = (int)bfs_order.size() - 1; i >= 0; i--) {
            int curr = bfs_order[i];
            freq[t[curr].link] += freq[curr];
        }

        vector<int> ans(num_patterns, 0);
        for (int i = 1; i < t.size(); i++) {
            if (freq[i] > 0) {
                for (int id : t[i].word_ids) {
                    ans[id] += freq[i];
                }
            }
        }
        return ans;
    }

    // 2. Devuelve un vector donde ans[id] es el índice de inicio de la PRIMERA aparición
    // Si el patrón no aparece, ans[id] será -1.
    vector<int> first_occurrences(const string& text) {
        vector<int> first_end(t.size(), INF);
        int u = 0;
        for (int i = 0; i < text.size(); i++) {
            u = t[u].next[text[i] - 'a'];
            first_end[u] = min(first_end[u], i);
        }

        // DP Lazy: Propagar mínimos de hojas a raíz
        for (int i = (int)bfs_order.size() - 1; i >= 0; i--) {
            int curr = bfs_order[i];
            first_end[t[curr].link] = min(first_end[t[curr].link], first_end[curr]);
        }

        vector<int> ans(num_patterns, -1);
        for (int i = 1; i < t.size(); i++) {
            if (first_end[i] != INF) {
                for (int id : t[i].word_ids) {
                    int start_idx = first_end[i] - pat_len[id] + 1;
                    if (ans[id] == -1) ans[id] = start_idx;
                    else ans[id] = min(ans[id], start_idx);
                }
            }
        }
        return ans;
    }
};

/**
 * 8. Suffix Automaton (Autómata de Sufijos - SAM)
 * Construcción: O(N) tiempo y O(N) memoria.
 * Resuelve: Subcadenas distintas, primera aparición, conteo de ocurrencias.
 * Nota: El autómata comprime todos los sufijos en máximo 2N-1 estados.
 */
struct SuffixAutomaton {
    struct Node {
        int len = 0;        // Longitud del sufijo más largo en este estado
        int link = -1;      // Enlace de fallo (hacia el sufijo más largo que se separó)
        int next[26];       // Transiciones hacia adelante (tren de lectura)

        int firstpos = -1;  // Índice final de la primera vez que aparece
        ll cnt = 0;         // Tamaño del conjunto EndPos (frecuencia de aparición)

        Node() { fill(next, next + 26, -1); }
    };

    vector<Node> t;
    int last;
    bool occ_computed = false; // Bandera para la DP de frecuencias

    // Constructor: Arma el autómata de forma Online (letra por letra)
    SuffixAutomaton(const string& s = "") {
        t.reserve(s.size() * 2 + 1);
        t.emplace_back();
        t[0].len = 0;
        t[0].link = -1;
        last = 0;

        for (char c : s) {
            extend(c - 'a');
            // Nota: Si hay mayúsculas/símbolos, usa 'extend(c)' y cambia next[26] por map<char, int>
        }
    }

    void extend(int c) {
        int cur = t.size();
        t.emplace_back();

        t[cur].len = t[last].len + 1;
        t[cur].firstpos = t[cur].len - 1;
        t[cur].cnt = 1; // Todo prefijo real del texto aporta 1 aparición inicial

        int p = last;
        // Subimos por los links buscando si alguien ya conocía esta letra
        while (p != -1 && t[p].next[c] == -1) {
            t[p].next[c] = cur;
            p = t[p].link;
        }

        if (p == -1) {
            // Caso base: La letra es totalmente nueva
            t[cur].link = 0;
        } else {
            int q = t[p].next[c];
            if (t[p].len + 1 == t[q].len) {
                // Caso Ideal: Las longitudes encajan, q es el link perfecto
                t[cur].link = q;
            } else {
                // Caso Difícil (Ruptura de EndPos): Creamos el CLON
                int clone = t.size();
                t.push_back(t[q]); // Copiamos transiciones y link de q

                t[clone].len = t[p].len + 1; // Forzamos la longitud exclusiva
                t[clone].cnt = 0;            // Los clones NO aportan aparición base (no son prefijos reales)
                // El firstpos del clon se hereda de q automáticamente al copiar

                // Redirigimos el tráfico que apuntaba a q para que apunte al clon
                while (p != -1 && t[p].next[c] == q) {
                    t[p].next[c] = clone;
                    p = t[p].link;
                }
                // Ajustamos las redes de seguridad (links)
                t[q].link = t[cur].link = clone;
            }
        }
        last = cur;
    }

    // DP Topológica - O(N): Propaga las frecuencias de hojas a raíz por los links
    // OBLIGATORIO llamarlo antes de count_occurrences()
    void compute_occurrences() {
        if (occ_computed) return;
        occ_computed = true;

        int n = t.size();
        // Ordenamiento lineal (Counting Sort) basado en la longitud (len)
        vector<int> c(n + 1, 0), order(n);
        for (int i = 0; i < n; i++) c[t[i].len]++;
        for (int i = 1; i <= n; i++) c[i] += c[i - 1];
        for (int i = 0; i < n; i++) order[--c[t[i].len]] = i;

        // Propagamos 'cnt' de los estados más largos a los más cortos
        for (int i = n - 1; i > 0; i--) {
            int u = order[i];
            if (t[u].link != -1) {
                t[t[u].link].cnt += t[u].cnt;
            }
        }
    }

    // 1. Número de subcadenas distintas (Number of Different Substrings)
    ll distinct_substrings() {
        ll ans = 0;
        for (int i = 1; i < t.size(); i++) {
            ans += t[i].len - t[t[i].link].len;
        }
        return ans;
    }

    // 2. Retorna el índice de INICIO de la primera aparición del patrón 'p'
    // Retorna -1 si no existe. Tiempo: O(|p|)
    int first_occurrence(const string& p) {
        int u = 0;
        for (char ch : p) {
            int c = ch - 'a';
            if (t[u].next[c] == -1) return -1;
            u = t[u].next[c];
        }
        return t[u].firstpos - p.size() + 1;
    }

    // 3. Retorna cuántas veces aparece exactamente el patrón 'p' en el texto
    // Tiempo: O(|p|)
    ll count_occurrences(const string& p) {
        compute_occurrences();
        int u = 0;
        for (char ch : p) {
            int c = ch - 'a';
            if (t[u].next[c] == -1) return 0;
            u = t[u].next[c];
        }
        return t[u].cnt;
    }
};

/**
 * 9. Suffix Array (Arreglo de Sufijos) + LCP Array
 * Construcción: SA en O(N log N) con Radix Sort, LCP en O(N) con Kasai.
 * Agrega automáticamente el centinela '$' (ASCII menor que 'a'-'z').
 * p[i]   = Índice donde comienza el i-ésimo sufijo lexicográficamente menor.
 * lcp[i] = Longitud del prefijo común más largo entre p[i] y p[i-1].
 */
struct SuffixArray {
    string s;
    int n;
    vector<int> p, c, lcp;

    SuffixArray(string _s) : s(_s + "$"), n(s.size()) {
        build_sa();
        build_lcp();
    }

    void build_sa() {
        const int alphabet = 256;
        p.assign(n, 0);
        c.assign(n, 0);
        vector<int> cnt(max(alphabet, n), 0), p_new(n), c_new(n);

        // K = 0 (Ordenar por el primer caracter)
        for (int i = 0; i < n; i++) cnt[s[i]]++;
        for (int i = 1; i < alphabet; i++) cnt[i] += cnt[i - 1];
        for (int i = 0; i < n; i++) p[--cnt[s[i]]] = i;
        
        c[p[0]] = 0;
        int classes = 1;
        for (int i = 1; i < n; i++) {
            if (s[p[i]] != s[p[i - 1]]) classes++;
            c[p[i]] = classes - 1;
        }

        // Transiciones K -> K + 1 (Longitudes 2^k)
        for (int k = 0; (1 << k) < n; ++k) {
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

    // Algoritmo de Kasai para construir el LCP en O(N)
    void build_lcp() {
        lcp.assign(n, 0);
        vector<int> rank(n, 0);
        for (int i = 0; i < n; i++) rank[p[i]] = i;
        
        int k = 0;
        for (int i = 0; i < n - 1; i++) { // n-1 para ignorar el sufijo "$"
            int pi = rank[i];
            int j = p[pi - 1];
            // Expandimos el match actual
            while (s[i + k] == s[j + k]) k++;
            lcp[pi] = k;
            if (k > 0) k--; // Descontamos el primer caracter para el siguiente sufijo
        }
    }
};


int main() {
    // 1. Optimización rápida
    ios_base::sync_with_stdio(0); cin.tie(0);

    // ==========================================
    // EJEMPLOS DE USO MÁGICOS
    // ==========================================

    // 1. Z-Algorithm
    // string s = "abracadabra";
    // vector<int> z = StringAlgo::z_function(s);

    // 2. KMP (búsqueda de patrones)
    // string s = "ababcababc";
    // vector<int> pi = StringAlgo::prefix_function(s);

    // 3. KMP Stack (eliminar todas las ocurrencias)
    // string s = "ababcababc";
    // string p = "abc";
    // string result = StringAlgo::remove_all_occurrences(s, p);
    
    // 4. Manacher (palíndromo más largo)
    // string s = "babad";
    // vector<int> p = StringAlgo::manacher(s);
    // int max_len = 0, center_idx = 0;
    // for(int i = 0; i < p.size(); i++) {
    //     if(p[i] > max_len) {
    //         max_len = p[i];
    //         center_idx = i;
    //     }
    // }
    // int start = (center_idx - max_len) / 2;
    // // cout << "Max palindrome: " << s.substr(start, max_len) << "\n";

    // 5. Hash O(1)
    // StringHash sh(s);
    // if(sh.get_hash(0, 3) == sh.get_hash(7, 10)) {
    //     // ¡Detectó que "abra" es igual a "abra" en O(1)!
    // }

    // 6. Trie
    // Trie trie;
    // trie.insert("codeforces");
    // if(trie.search("codeforces")) cout << "YES\n";
    
    // 7. Aho-Corasick (múltiples patrones sobre un texto)
    // AhoCorasick ac;
    // ac.insert("aba");
    // ac.insert("ba");
    // ac.build();
    // vector<int> cnt = ac.count_occurrences("ababa");
    // vector<int> first = ac.first_occurrences("ababa");

    // 8. Suffix Automaton
    // SuffixAutomaton sam(s);
    // ll distinct = sam.distinct_substrings();
    // ll veces = sam.count_occurrences("aba");
    // int pos = sam.first_occurrence("aba");

    // 9. Suffix Array + LCP
    // string s = "banana";
    // SuffixArray sa(s);
    //
    // // a) Imprimir los sufijos ordenados (Ignorando i=0 que es "$")
    // // for(int i = 1; i < sa.n; i++) {
    // //     cout << sa.p[i] << ": " << sa.s.substr(sa.p[i]) << "\n";
    // // }
    //
    // // b) Substring repetido más largo (Max del LCP array)
    // // int max_lcp = 0, start_idx = 0;
    // // for(int i = 1; i < sa.n; i++) {
    // //     if(sa.lcp[i] > max_lcp) {
    // //         max_lcp = sa.lcp[i];
    // //         start_idx = sa.p[i];
    // //     }
    // // }
    // // cout << "Repetido mas largo: " << sa.s.substr(start_idx, max_lcp) << "\n";

    return 0;
}