/**
 * @file StringTemplateICPC.cpp
 * @brief Plantilla de algoritmos para cadenas y patrones.
 * @details Incluye utilidades básicas (split, trim, replace), KMP, Z-function, hashing doble, Trie, Aho-Corasick, Manacher y Suffix Automaton.
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

    // ==========================================
    // 0. UTILIDADES BÁSICAS (O(N))
    // ==========================================
    
    inline void to_lower(string &s) { transform(s.begin(), s.end(), s.begin(), ::tolower); }
    inline void to_upper(string &s) { transform(s.begin(), s.end(), s.begin(), ::toupper); }
    
    inline void trim(string &s) {
        s.erase(s.begin(), find_if(s.begin(), s.end(), [](unsigned char ch) { return !isspace(ch); }));
        s.erase(find_if(s.rbegin(), s.rend(), [](unsigned char ch) { return !isspace(ch); }).base(), s.end());
    }

    inline void replace_all(string &s, const string &old_val, const string &new_val) {
        size_t pos = 0;
        while ((pos = s.find(old_val, pos)) != string::npos) {
            s.replace(pos, old_val.length(), new_val);
            pos += new_val.length();
        }
    }

    inline vector<string> split(const string &s, char delim) {
        vector<string> elems; string item;
        istringstream iss(s);
        while (getline(iss, item, delim)) elems.push_back(item);
        return elems;
    }

    inline string join(const vector<string> &v, const string &delim) {
        string res;
        for (size_t i = 0; i < v.size(); ++i) {
            res += v[i];
            if (i < v.size() - 1) res += delim;
        }
        return res;
    }

    // ==========================================
    // ALGORITMOS DE BÚSQUEDA Y PATRONES
    // ==========================================

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
     */
    string remove_all_occurrences(const string& s, const string& p) {
        vector<int> pi_p = prefix_function(p);
        string ans = "";
        vector<int> estado(s.size() + 1, 0);

        for (char c : s) {
            ans.push_back(c);
            int j = estado[ans.size() - 1];
            while (j > 0 && p[j] != c) j = pi_p[j - 1];
            if (p[j] == c) j++;
            estado[ans.size()] = j;
            if (j == p.size()) ans.resize(ans.size() - p.size());
        }
        return ans;
    }

    /**
     * 4. Rabin-Karp (Búsqueda de patrón único) - O(N + M)
     * Retorna índices de inicio donde 'pattern' aparece en 'text'.
     */
    vector<int> rabin_karp(const string& text, const string& pattern) {
        const int p = 31;
        const ll m = 1e9 + 9;
        int S = pattern.size(), T = text.size();
        if (S > T || S == 0) return {};

        vector<ll> p_pow(max(S, T));
        p_pow[0] = 1;
        for (int i = 1; i < (int)p_pow.size(); i++)
            p_pow[i] = (p_pow[i - 1] * p) % m;

        vector<ll> h(T + 1, 0);
        for (int i = 0; i < T; i++)
            h[i + 1] = (h[i] + (text[i] - 'a' + 1) * p_pow[i]) % m;
            
        ll h_s = 0;
        for (int i = 0; i < S; i++)
            h_s = (h_s + (pattern[i] - 'a' + 1) * p_pow[i]) % m;

        vector<int> occ;
        for (int i = 0; i + S - 1 < T; i++) {
            ll cur_h = (h[i + S] + m - h[i]) % m;
            if (cur_h == h_s * p_pow[i] % m) occ.push_back(i);
        }
        return occ;
    }

    /**
     * 5. Algoritmo de Manacher - O(N)
     * Encuentra la longitud del palíndromo más largo centrado en cada posición.
     */
    vector<int> manacher(const string& s) {
        string t = "^#";
        for (char c : s) { t += c; t += '#'; }
        t += '$'; 

        int n = t.size(), c = 0, r = 0; 
        vector<int> p(n, 0);

        for (int i = 1; i < n - 1; i++) {
            int i_mirror = 2 * c - i;
            if (r > i) p[i] = min(r - i, p[i_mirror]);
            else p[i] = 0;

            while (t[i + 1 + p[i]] == t[i - 1 - p[i]]) p[i]++;

            if (i + p[i] > r) {
                c = i; r = i + p[i];
            }
        }
        return vector<int>(p.begin() + 1, p.end() - 1);
    }
}

/**
 * 6. Double Rolling Hash (O(N) build, O(1) query)
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

    uint64_t get_hash(int l, int r) {
        ll hash1 = (h1[r + 1] - (h1[l] * p1[r - l + 1]) % MOD1 + MOD1) % MOD1;
        ll hash2 = (h2[r + 1] - (h2[l] * p2[r - l + 1]) % MOD2 + MOD2) % MOD2;
        return ((uint64_t)hash1 << 32) | (uint32_t)hash2;
    }
};

/**
 * 7. Trie (Árbol de Prefijos)
 */
struct Trie {
    struct Node {
        int next[26];
        bool is_end = false;
        int count = 0; 
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
 * 8. Aho-Corasick (Búsqueda Simultánea de Múltiples Patrones)
 */
struct AhoCorasick {
    struct Node {
        int next[26];
        int link = 0;       
        vector<int> word_ids; 
        Node() { fill(next, next + 26, -1); }
    };

    vector<Node> t;
    vector<int> bfs_order; 
    vector<int> pat_len;   
    int num_patterns = 0;

    AhoCorasick() { t.emplace_back(); }

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
            int u = q.front(); q.pop();
            bfs_order.push_back(u);

            for (int c = 0; c < 26; c++) {
                int v = t[u].next[c];
                if (v != -1) {
                    t[v].link = t[t[u].link].next[c];
                    q.push(v);
                } else {
                    t[u].next[c] = t[t[u].link].next[c]; 
                }
            }
        }
    }

    vector<int> count_occurrences(const string& text) {
        vector<int> freq(t.size(), 0);
        int u = 0;
        for (char ch : text) {
            u = t[u].next[ch - 'a'];
            freq[u]++;
        }
        for (int i = (int)bfs_order.size() - 1; i >= 0; i--) {
            int curr = bfs_order[i];
            freq[t[curr].link] += freq[curr];
        }
        vector<int> ans(num_patterns, 0);
        for (int i = 1; i < t.size(); i++) {
            if (freq[i] > 0) {
                for (int id : t[i].word_ids) ans[id] += freq[i];
            }
        }
        return ans;
    }

    vector<int> first_occurrences(const string& text) {
        vector<int> first_end(t.size(), INF);
        int u = 0;
        for (int i = 0; i < text.size(); i++) {
            u = t[u].next[text[i] - 'a'];
            first_end[u] = min(first_end[u], i);
        }
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
 * 9. Suffix Automaton (Autómata de Sufijos - SAM)
 */
struct SuffixAutomaton {
    struct Node {
        int len = 0;        
        int link = -1;      
        int next[26];       
        int firstpos = -1;  
        ll cnt = 0;         
        Node() { fill(next, next + 26, -1); }
    };

    vector<Node> t;
    int last;
    bool occ_computed = false; 

    SuffixAutomaton(const string& s = "") {
        t.reserve(s.size() * 2 + 1);
        t.emplace_back();
        t[0].len = 0; t[0].link = -1; last = 0;
        for (char c : s) extend(c - 'a');
    }

    void extend(int c) {
        int cur = t.size(); t.emplace_back();
        t[cur].len = t[last].len + 1;
        t[cur].firstpos = t[cur].len - 1;
        t[cur].cnt = 1; 

        int p = last;
        while (p != -1 && t[p].next[c] == -1) {
            t[p].next[c] = cur; p = t[p].link;
        }

        if (p == -1) {
            t[cur].link = 0;
        } else {
            int q = t[p].next[c];
            if (t[p].len + 1 == t[q].len) {
                t[cur].link = q;
            } else {
                int clone = t.size();
                t.push_back(t[q]); 
                t[clone].len = t[p].len + 1; 
                t[clone].cnt = 0;            
                while (p != -1 && t[p].next[c] == q) {
                    t[p].next[c] = clone; p = t[p].link;
                }
                t[q].link = t[cur].link = clone;
            }
        }
        last = cur;
    }

    void compute_occurrences() {
        if (occ_computed) return;
        occ_computed = true;
        int n = t.size();
        vector<int> c(n + 1, 0), order(n);
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
            if (t[u].next[c] == -1) return -1;
            u = t[u].next[c];
        }
        return t[u].firstpos - p.size() + 1;
    }

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
* 10. Suffix Array + LCP Array + RMQ
*/
struct SuffixArray {
    string s;
    int n, log_n;
    vector<int> p, c, lcp, rank;
    vector<vector<int>> st; 

    SuffixArray(string _s) : s(_s + "$"), n(s.size()) {
        build_sa(); build_lcp(); build_rmq();
    }

    void build_sa() {
        const int alphabet = 256;
        p.assign(n, 0); c.assign(n, 0);
        vector<int> cnt(max(alphabet, n), 0), p_new(n), c_new(n);

        for (int i = 0; i < n; i++) cnt[s[i]]++;
        for (int i = 1; i < alphabet; i++) cnt[i] += cnt[i - 1];
        for (int i = 0; i < n; i++) p[--cnt[s[i]]] = i;
        
        c[p[0]] = 0;
        int classes = 1;
        for (int i = 1; i < n; i++) {
            if (s[p[i]] != s[p[i - 1]]) classes++;
            c[p[i]] = classes - 1;
        }

        for (int k = 0; (1 << k) < n; ++k) {
            for (int i = 0; i < n; i++) {
                p_new[i] = p[i] - (1 << k);
                if (p_new[i] < 0) p_new[i] += n;
            }
            fill(cnt.begin(), cnt.begin() + classes, 0);
            for (int i = 0; i < n; i++) cnt[c[p_new[i]]]++;
            for (int i = 1; i < classes; i++) cnt[i] += cnt[i - 1];
            for (int i = n - 1; i >= 0; i--) p[--cnt[c[p_new[i]]]] = p_new[i];

            c_new[p[0]] = 0; classes = 1;
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
        lcp.assign(n, 0); rank.assign(n, 0);
        for (int i = 0; i < n; i++) rank[p[i]] = i; 
        
        int k = 0;
        for (int i = 0; i < n - 1; i++) { 
            int pi = rank[i];
            int j = p[pi - 1];
            while (s[i + k] == s[j + k]) k++;
            lcp[pi] = k;
            if (k > 0) k--; 
        }
    }

    void build_rmq() {
        log_n = log2(n) + 1;
        st.assign(n, vector<int>(log_n));
        for (int i = 0; i < n; i++) st[i][0] = lcp[i];
        for (int j = 1; j < log_n; j++) {
            for (int i = 0; i + (1 << j) <= n; i++) {
                st[i][j] = min(st[i][j - 1], st[i + (1 << (j - 1))][j - 1]);
            }
        }
    }

    int get_lcp(int i, int j) {
        if (i == j) return n - 1 - i; 
        int u = rank[i], v = rank[j];
        if (u > v) swap(u, v);
        u++; 
        int k = log2(v - u + 1);
        return min(st[u][k], st[v - (1 << k) + 1][k]);
    }
};


int main() {
    ios_base::sync_with_stdio(0); cin.tie(0);

    // ==========================================
    // EJEMPLOS DE USO MÁGICOS
    // ==========================================

    /*
    // --- 0. Utilidades Básicas ---
    string input = "  manzana, pera, uva  ";
    StringAlgo::trim(input);
    StringAlgo::to_upper(input);
    StringAlgo::replace_all(input, " ", "");
    vector<string> frutas = StringAlgo::split(input, ',');
    cout << "Unidas de nuevo: " << StringAlgo::join(frutas, " - ") << "\n";
    
    // --- 4. Rabin-Karp ---
    vector<int> ocurrencias = StringAlgo::rabin_karp("abracadabra", "abra");
    // ocurrencias = {0, 7}
    */

    return 0;
}