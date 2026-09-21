
/**
 * @file RQTemplateICPC.cpp
 * @brief Plantilla de estructuras para consultas y actualizaciones en rangos.
 * @details Incluye Sparse Table, Fenwick Tree y varias variantes de Segment Tree.
 * @note Cambia la operacion de combinacion y el valor neutro segun el problema.
 */
//   ____ ___  ____  _____   ____  _   _
//  / ___/ _ \|  _ \| ____| / ___|| | | |
// | |  | | | | | | |  _|   \___ \| | | |
// | |__| |_| | |_| | |___   ___) | |_| |
//  \____\___/|____/|_____| |____/ \___/
//
//                RANGE QUERIES TEMPLATE

#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
const ll LINF = 1e18;

// ====================================================================
// SPARSE TABLE (Consultas en Rango Estáticas)
// ====================================================================
class SparseTable {
private:
    int n, log_n;
    vector<vector<ll>> st;

public:
    /**
     * Inicializa y construye la Sparse Table.
     * @param a El arreglo de entrada (0-indexed).
     */
    SparseTable(const vector<ll>& a) {
        n = a.size();
        log_n = log2(n) + 1;
        st.assign(n, vector<ll>(log_n));

        for (int i = 0; i < n; i++) {
            st[i][0] = a[i];
        }

        for (int j = 1; j < log_n; j++) {
            for (int i = 0; i + (1 << j) <= n; i++) {
                // Cambiar min a max o __gcd dependiendo del problema
                st[i][j] = min(st[i][j - 1], st[i + (1 << (j - 1))][j - 1]);
            }
        }
    }

    /**
     * Consulta el valor mínimo en el rango [L, R].
     * @param L Índice izquierdo del rango (0-indexed, inclusivo).
     * @param R Índice derecho del rango (0-indexed, inclusivo).
     * @return El valor mínimo en el rango.
     */
    ll query(int L, int R) {
        int j = log2(R - L + 1);
        // Cambiar min a max o __gcd dependiendo del problema
        return min(st[L][j], st[R - (1 << j) + 1][j]);
    }
};

// ====================================================================
// FENWICK TREE (Binary Indexed Tree)
// ====================================================================
class FenwickTree {
private:
    int n;
    vector<ll> bit;

public:
    /**
     * Inicializa un Fenwick Tree de tamaño n.
     * @param n Tamaño del árbol (internamente el arreglo es 1-indexed).
     */
    FenwickTree(int n) : n(n), bit(n + 1, 0) {}

    /**
     * Suma 'val' al elemento en el índice 'idx'.
     * @param idx Índice a actualizar (1-indexed).
     * @param val Valor a sumar.
     */
    void add(int idx, ll val) {
        for (; idx <= n; idx += idx & -idx) {
            bit[idx] += val;
        }
    }

    /**
     * Consulta la suma del prefijo desde 1 hasta 'idx'.
     * @param idx Índice final del prefijo (1-indexed).
     * @return Suma en el rango [1, idx].
     */
    ll sum(int idx) {
        ll res = 0;
        for (; idx > 0; idx -= idx & -idx) {
            res += bit[idx];
        }
        return res;
    }

    /**
     * Consulta la suma en el rango [L, R].
     * @param L Índice izquierdo (1-indexed, inclusivo).
     * @param R Índice derecho (1-indexed, inclusivo).
     * @return Suma en el rango [L, R].
     */
    ll query(int L, int R) {
        return sum(R) - sum(L - 1);
    }
};

// ====================================================================
// ITERATIVE SEGMENT TREE (Point Update, Range Query)
// ====================================================================
class IterativeSegTree {
private:
    int n;
    vector<ll> st;
    const ll NEUTRAL = 0; // Cambiar a LINF para min, -LINF para max

    ll combine(ll a, ll b) {
        return a + b; // Cambiar a min(a,b) o max(a,b)
    }

public:
    /**
     * Inicializa y construye el Segment Tree Iterativo.
     * @param a El arreglo de entrada (0-indexed).
     */
    IterativeSegTree(const vector<ll>& a) {
        n = a.size();
        st.assign(2 * n, NEUTRAL);
        for (int i = 0; i < n; i++) {
            st[n + i] = a[i];
        }
        for (int i = n - 1; i > 0; i--) {
            st[i] = combine(st[i << 1], st[i << 1 | 1]);
        }
    }

    /**
     * Actualiza el elemento en el índice 'pos' al valor 'val'.
     * @param pos Índice a actualizar (0-indexed).
     * @param val El nuevo valor.
     */
    void update(int pos, ll val) {
        for (st[pos += n] = val; pos > 1; pos >>= 1) {
            st[pos >> 1] = combine(st[pos], st[pos ^ 1]);
        }
    }

    /**
     * Consulta el valor combinado en el rango [L, R].
     * @param L Índice izquierdo (0-indexed, inclusivo).
     * @param R Índice derecho (0-indexed, inclusivo).
     * @return Resultado combinado en el rango.
     */
    ll query(int L, int R) {
        ll resL = NEUTRAL, resR = NEUTRAL;
        for (L += n, R += n + 1; L < R; L >>= 1, R >>= 1) {
            if (L & 1) resL = combine(resL, st[L++]);
            if (R & 1) resR = combine(st[--R], resR);
        }
        return combine(resL, resR);
    }
};

// ====================================================================
// RECURSIVE SEGMENT TREE (Point Update, Range Query)
// Ideal para búsqueda binaria sobre el ST (find_first)
// ====================================================================
class RecursiveSegTree {
private:
    int n;
    vector<ll> st;
    const ll NEUTRAL = 0; // Cambiar según operación (0 para max, LINF para min)

    ll combine(ll a, ll b) {
        return max(a, b); // Cambiar según problema. Para Hotel Queries es max
    }

    void build(int p, int L, int R, const vector<ll>& a) {
        if (L == R) {
            st[p] = a[L];
            return;
        }
        int mid = L + (R - L) / 2;
        build(p << 1, L, mid, a);
        build((p << 1) | 1, mid + 1, R, a);
        st[p] = combine(st[p << 1], st[(p << 1) | 1]);
    }

    void update(int p, int L, int R, int pos, ll val) {
        if (L == R) {
            st[p] += val; // Si quieres reasignar valor, usa st[p] = val
            return;
        }
        int mid = L + (R - L) / 2;
        if (pos <= mid) update(p << 1, L, mid, pos, val);
        else update((p << 1) | 1, mid + 1, R, pos, val);
        st[p] = combine(st[p << 1], st[(p << 1) | 1]);
    }

    ll query(int p, int L, int R, int qL, int qR) {
        if (qL > R || qR < L) return NEUTRAL;
        if (qL <= L && R <= qR) return st[p];
        int mid = L + (R - L) / 2;
        return combine(query(p << 1, L, mid, qL, qR),
                       query((p << 1) | 1, mid + 1, R, qL, qR));
    }

    // Retorna el índice de la primera hoja cuyo valor cumpla con la condición >= val
    // Retorna -1 si no existe. ASUME QUE EL ÁRBOL ES DE MÁXIMOS.
    int find_first(int p, int L, int R, ll val) {
        if (st[p] < val) return -1; // Ningún elemento en este subárbol cumple
        if (L == R) return L;       // ¡Encontramos la hoja!
        
        int mid = L + (R - L) / 2;
        // Revisamos primero el hijo izquierdo para garantizar que sea el PRIMER índice
        if (st[p << 1] >= val) {
            return find_first(p << 1, L, mid, val);
        } else {
            return find_first((p << 1) | 1, mid + 1, R, val);
        }
    }

public:
    RecursiveSegTree(const vector<ll>& a) {
        n = a.size();
        st.assign(4 * n, NEUTRAL);
        build(1, 0, n - 1, a);
    }

    void update(int pos, ll val) {
        update(1, 0, n - 1, pos, val);
    }

    ll query(int L, int R) {
        return query(1, 0, n - 1, L, R);
    }

    int find_first(ll val) {
        return find_first(1, 0, n - 1, val);
    }
};

// ====================================================================
// LAZY SEGMENT TREE (Range Update, Range Query)
// ====================================================================
class LazySegTree {
private:
    int n;
    vector<ll> st, lazy_add, lazy_set;
    vector<bool> marked;
    const ll NEUTRAL = 0; // Cambiar a LINF para consultas de min

    // --- 1. OPERACIÓN DE CONSULTA ---
    ll combine(ll a, ll b) {
        return a + b; // Cambiar a min(a,b) o max(a,b) según el problema
    }

    // --- 2. OPERACIÓN DE ASIGNACIÓN (SET) ---
    void apply_set(int p, ll val, int L, int R) {
        // Para RMQ (Min/Max): st[p] = val;
        // Para Suma:
        st[p] = val * (R - L + 1); 
        
        lazy_set[p] = val;
        lazy_add[p] = 0; // El SET borra cualquier ADD pendiente
        marked[p] = true;
    }

    // --- 3. OPERACIÓN DE SUMA (ADD) ---
    void apply_add(int p, ll val, int L, int R) {
        // Para RMQ (Min/Max): st[p] += val;
        // Para Suma:
        st[p] += val * (R - L + 1);
        
        if (marked[p]) {
            lazy_set[p] += val; // Si ya había un SET, se vuelve un SET mayor
        } else {
            lazy_add[p] += val; // Si no, es un ADD normal
        }
    }

    // --- LOGICA INTERNA (No tocar en competencia) ---
    void push(int p, int L, int R) {
        int mid = L + (R - L) / 2;
        int left = p << 1;
        int right = (p << 1) | 1;

        if (marked[p]) {
            apply_set(left, lazy_set[p], L, mid);
            apply_set(right, lazy_set[p], mid + 1, R);
            marked[p] = false;
        }
        if (lazy_add[p] != 0) {
            apply_add(left, lazy_add[p], L, mid);
            apply_add(right, lazy_add[p], mid + 1, R);
            lazy_add[p] = 0;
        }
    }

    void build(int p, int L, int R, const vector<ll>& a) {
        if (L == R) {
            st[p] = a[L];
            return;
        }
        int mid = L + (R - L) / 2;
        build(p << 1, L, mid, a);
        build((p << 1) | 1, mid + 1, R, a);
        st[p] = combine(st[p << 1], st[(p << 1) | 1]);
    }

    void update(int p, int L, int R, int qL, int qR, ll val, int type) {
        if (qL > R || qR < L) return;
        if (qL <= L && R <= qR) {
            if (type == 1) apply_add(p, val, L, R); // Tipo 1: Sumar
            if (type == 2) apply_set(p, val, L, R); // Tipo 2: Asignar
            return;
        }
        push(p, L, R);
        int mid = L + (R - L) / 2;
        update(p << 1, L, mid, qL, qR, val, type);
        update((p << 1) | 1, mid + 1, R, qL, qR, val, type);
        st[p] = combine(st[p << 1], st[(p << 1) | 1]);
    }

    ll query(int p, int L, int R, int qL, int qR) {
        if (qL > R || qR < L) return NEUTRAL;
        if (qL <= L && R <= qR) return st[p];
        push(p, L, R);
        int mid = L + (R - L) / 2;
        return combine(query(p << 1, L, mid, qL, qR),
                       query((p << 1) | 1, mid + 1, R, qL, qR));
    }

    ll get_point(int p, int L, int R, int pos) {
        if (L == R) return st[p];
        push(p, L, R);
        int mid = L + (R - L) / 2;
        if (pos <= mid) return get_point(p << 1, L, mid, pos);
        else return get_point((p << 1) | 1, mid + 1, R, pos);
    }

public:
    LazySegTree(const vector<ll>& a) {
        n = a.size();
        st.assign(4 * n, NEUTRAL);
        lazy_add.assign(4 * n, 0);
        lazy_set.assign(4 * n, 0);
        marked.assign(4 * n, false);
        build(1, 0, n - 1, a);
    }

    void add_range(int L, int R, ll val) { update(1, 0, n - 1, L, R, val, 1); }
    void set_range(int L, int R, ll val) { update(1, 0, n - 1, L, R, val, 2); }
    ll query(int L, int R) { return query(1, 0, n - 1, L, R); }
    
    // Novedad: Consulta el valor final de un solo nodo en O(log N) directo
    ll get(int pos) { return get_point(1, 0, n - 1, pos); }
};

// ====================================================================
// DYNAMIC SEGMENT TREE (Point Update, Range Query on Large Ranges)
// ====================================================================
class DynamicSegTree {
private:
    struct Node {
        ll sum;
        int left, right;
        Node() : sum(0), left(0), right(0) {}
    };

    vector<Node> st;
    ll MAX_RANGE;
    int root;

    int new_node() {
        st.emplace_back();
        return st.size() - 1;
    }

    int update(int node, ll L, ll R, ll pos, ll val) {
        if (!node) node = new_node();
        st[node].sum += val; // O cambiar a combinación para min/max
        if (L == R) return node;
        
        ll mid = L + (R - L) / 2;
        if (pos <= mid) {
            int left_child = update(st[node].left, L, mid, pos, val);
            st[node].left = left_child; // Evita bugs de memoria
        } else {
            int right_child = update(st[node].right, mid + 1, R, pos, val);
            st[node].right = right_child;
        }
        return node;
    }

    ll query(int node, ll L, ll R, ll qL, ll qR) {
        if (!node || qL > R || qR < L) return 0; // Cambiar a NEUTRAL según la operación
        if (qL <= L && R <= qR) return st[node].sum;
        
        ll mid = L + (R - L) / 2;
        return query(st[node].left, L, mid, qL, qR) + 
               query(st[node].right, mid + 1, R, qL, qR);
    }

public:
    /**
     * Inicializa el Segment Tree Dinámico.
     * @param max_range El valor máximo que puede tomar el índice (ej. 1e9).
     */
    DynamicSegTree(ll max_range = 1e9) : MAX_RANGE(max_range) {
        st.emplace_back(); // Nodo 0 es el nulo/dummy
        root = 0;
    }

    /**
     * Suma 'val' a la posición 'pos' (puede ser valor negativo para restar).
     * @param pos Índice a actualizar (hasta max_range).
     * @param val Valor a sumar.
     */
    void update(ll pos, ll val) {
        root = update(root, 1, MAX_RANGE, pos, val);
    }

    /**
     * Consulta la suma en el rango [L, R].
     * @param L Índice izquierdo (inclusivo).
     * @param R Índice derecho (inclusivo).
     * @return Suma en el rango [L, R].
     */
    ll query(ll L, ll R) {
        return query(root, 1, MAX_RANGE, L, R);
    }
};

// ====================================================================
// PERSISTENT SEGMENT TREE (Point Update, Range Query histórico)
// ====================================================================
class PersistentSegTree {
private:
    struct Node {
        ll sum;
        int left, right;
        Node(ll _sum = 0, int _l = 0, int _r = 0) : sum(_sum), left(_l), right(_r) {}
    };

    vector<Node> st;
    vector<int> roots; // Guarda las raíces de cada versión histórica
    int n;

    int new_node(ll sum, int left, int right) {
        st.emplace_back(sum, left, right);
        return st.size() - 1;
    }

    int build(int L, int R, const vector<ll>& a) {
        if (L == R) return new_node(a[L], 0, 0);
        int mid = L + (R - L) / 2;
        int left_child = build(L, mid, a);
        int right_child = build(mid + 1, R, a);
        return new_node(st[left_child].sum + st[right_child].sum, left_child, right_child);
    }

    int update(int old_node, int L, int R, int pos, ll val) {
        if (L == R) return new_node(st[old_node].sum + val, 0, 0); // Copia modificada
        int mid = L + (R - L) / 2;
        int left_child = st[old_node].left;
        int right_child = st[old_node].right;
        
        if (pos <= mid) left_child = update(st[old_node].left, L, mid, pos, val);
        else right_child = update(st[old_node].right, mid + 1, R, pos, val);
        
        return new_node(st[left_child].sum + st[right_child].sum, left_child, right_child);
    }

    ll query(int node, int L, int R, int qL, int qR) {
        if (!node || qL > R || qR < L) return 0;
        if (qL <= L && R <= qR) return st[node].sum;
        int mid = L + (R - L) / 2;
        return query(st[node].left, L, mid, qL, qR) + 
               query(st[node].right, mid + 1, R, qL, qR);
    }

public:
    PersistentSegTree(const vector<ll>& a) {
        n = a.size();
        st.emplace_back(); // Nodo 0 nulo
        roots.push_back(build(0, n - 1, a)); // Versión 0
    }

    // Crea una nueva versión actualizando 'pos' con 'val' basándose en 'prev_version'
    int update(int prev_version, int pos, ll val) {
        int new_root = update(roots[prev_version], 0, n - 1, pos, val);
        roots.push_back(new_root);
        return roots.size() - 1; // Retorna el ID de la nueva versión
    }

    // Consulta el rango [L, R] en una versión específica del tiempo
    ll query(int version, int L, int R) {
        return query(roots[version], 0, n - 1, L, R);
    }
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    return 0;
}
