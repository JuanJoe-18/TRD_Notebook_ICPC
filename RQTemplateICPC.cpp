
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

// ====================================================================
// FENWICK TREE 2D (Point Update, Submatrix Sum Query)
// ====================================================================
class FenwickTree2D {
private:
    int n, m;
    vector<vector<ll>> bit;

public:
    /**
     * Inicializa un Fenwick Tree 2D (Matriz).
     * @param n Número de filas (internamente es 1-indexed).
     * @param m Número de columnas (internamente es 1-indexed).
     */
    FenwickTree2D(int n, int m) : n(n), m(m) {
        bit.assign(n + 1, vector<ll>(m + 1, 0));
    }

    /**
     * Suma 'val' a la celda en la fila 'r' y columna 'c'.
     * @param r Fila (1-indexed).
     * @param c Columna (1-indexed).
     * @param val Valor a sumar.
     */
    void add(int r, int c, ll val) {
        for (int i = r; i <= n; i += i & -i) {
            for (int j = c; j <= m; j += j & -j) {
                bit[i][j] += val;
            }
        }
    }

    /**
     * Consulta la suma del sub-rectángulo desde (1, 1) hasta (r, c).
     * @param r Fila final (1-indexed).
     * @param c Columna final (1-indexed).
     * @return Suma en el prefijo 2D.
     */
    ll sum(int r, int c) {
        ll res = 0;
        for (int i = r; i > 0; i -= i & -i) {
            for (int j = c; j > 0; j -= j & -j) {
                res += bit[i][j];
            }
        }
        return res;
    }

    /**
     * Consulta la suma del sub-rectángulo entre (r1, c1) y (r2, c2).
     * @param r1 Fila superior izquierda (1-indexed, inclusivo).
     * @param c1 Columna superior izquierda (1-indexed, inclusivo).
     * @param r2 Fila inferior derecha (1-indexed, inclusivo).
     * @param c2 Columna inferior derecha (1-indexed, inclusivo).
     * @return Suma en el rectángulo.
     */
    ll query(int r1, int c1, int r2, int c2) {
        return sum(r2, c2) - sum(r1 - 1, c2) - sum(r2, c1 - 1) + sum(r1 - 1, c1 - 1);
    }
};

// ====================================================================
// MERGE SORT TREE (Static Range Queries for Order Statistics)
// ====================================================================
class MergeSortTree {
private:
    int n;
    vector<vector<ll>> st;

    void build(int p, int L, int R, const vector<ll>& a) {
        if (L == R) {
            st[p].push_back(a[L]);
            return;
        }
        int mid = L + (R - L) / 2;
        build(p << 1, L, mid, a);
        build((p << 1) | 1, mid + 1, R, a);

        // Combina los dos vectores ordenados de los hijos en O(R - L)
        st[p].resize(st[p << 1].size() + st[(p << 1) | 1].size());
        merge(st[p << 1].begin(), st[p << 1].end(),
              st[(p << 1) | 1].begin(), st[(p << 1) | 1].end(),
              st[p].begin());
    }

    int query_less_than(int p, int L, int R, int qL, int qR, ll val) {
        if (qL > R || qR < L) return 0;
        if (qL <= L && R <= qR) {
            // Encuentra la cantidad de elementos estrictamente menores que 'val'
            return lower_bound(st[p].begin(), st[p].end(), val) - st[p].begin();
        }
        int mid = L + (R - L) / 2;
        return query_less_than(p << 1, L, mid, qL, qR, val) +
               query_less_than((p << 1) | 1, mid + 1, R, qL, qR, val);
    }

    int query_less_equal(int p, int L, int R, int qL, int qR, ll val) {
        if (qL > R || qR < L) return 0;
        if (qL <= L && R <= qR) {
            // Encuentra la cantidad de elementos menores o iguales que 'val'
            return upper_bound(st[p].begin(), st[p].end(), val) - st[p].begin();
        }
        int mid = L + (R - L) / 2;
        return query_less_equal(p << 1, L, mid, qL, qR, val) +
               query_less_equal((p << 1) | 1, mid + 1, R, qL, qR, val);
    }

public:
    /**
     * Inicializa y construye el Merge Sort Tree.
     * @param a El arreglo de entrada (0-indexed).
     */
    MergeSortTree(const vector<ll>& a) {
        n = a.size();
        st.assign(4 * n, vector<ll>());
        build(1, 0, n - 1, a);
    }

    /**
     * Cuenta cuántos elementos en el rango [L, R] son estrictamente menores a 'val'.
     * @param L Índice izquierdo (0-indexed, inclusivo).
     * @param R Índice derecho (0-indexed, inclusivo).
     * @param val Valor límite.
     * @return Cantidad de elementos < val.
     */
    int count_less_than(int L, int R, ll val) {
        return query_less_than(1, 0, n - 1, L, R, val);
    }

    /**
     * Cuenta cuántos elementos en el rango [L, R] son menores o iguales a 'val'.
     * @param L Índice izquierdo (0-indexed, inclusivo).
     * @param R Índice derecho (0-indexed, inclusivo).
     * @param val Valor límite.
     * @return Cantidad de elementos <= val.
     */
    int count_less_equal(int L, int R, ll val) {
        return query_less_equal(1, 0, n - 1, L, R, val);
    }
};

// ====================================================================
// MO'S ALGORITHM (Offline Range Queries en O((N+Q) * sqrt(N)))
// ====================================================================
// Estructura para almacenar y ordenar las consultas offline.
struct MoQuery {
    int l, r, id, block_size;

    // Optimización Par/Impar para el ordenamiento de Mo
    bool operator<(const MoQuery& other) const {
        int b1 = l / block_size;
        int b2 = other.l / block_size;
        if (b1 != b2) return b1 < b2;
        return (b1 & 1) ? (r < other.r) : (r > other.r);
    }
};

class MosAlgorithm {
private:
    int n, q;
    vector<ll> a, ans;
    vector<MoQuery> queries;
    ll current_ans; // Ajustar tipo según el problema

    // --- LÓGICA ESPECÍFICA DEL PROBLEMA ---
    void add(int idx) {
        // TODO: Lógica para AÑADIR a[idx] al rango
        // ej: if (++freq[a[idx]] == 1) current_ans++;
    }

    void remove(int idx) {
        // TODO: Lógica para REMOVER a[idx] del rango
        // ej: if (--freq[a[idx]] == 0) current_ans--;
    }

public:
    /**
     * Inicializa el algoritmo de Mo.
     * @param arr El arreglo original (0-indexed).
     * @param num_queries La cantidad total de consultas que se procesarán.
     */
    MosAlgorithm(const vector<ll>& arr, int num_queries) {
        n = arr.size();
        q = num_queries;
        a = arr;
        ans.assign(q, 0);
        queries.reserve(q);
        current_ans = 0; // Inicializar neutro
    }

    /**
     * Agrega una consulta al sistema offline.
     * @param L Índice izquierdo (0-indexed).
     * @param R Índice derecho (0-indexed).
     * @param id Identificador original de la consulta (para recuperar el orden).
     */
    void add_query(int L, int R, int id) {
        // El block_size ideal clásico es N / sqrt(Q)
        int block_size = max(1, (int)(n / sqrt(max(1, q))));
        queries.push_back({L, R, id, block_size});
    }

    /**
     * Procesa todas las consultas y retorna el arreglo de respuestas ordenado.
     * @return Vector con las respuestas en el orden original de las consultas.
     */
    vector<ll> solve() {
        sort(queries.begin(), queries.end());

        int curr_l = 0, curr_r = -1; // Rango inicial vacío

        for (const auto& mq : queries) {
            // Expandir rango
            while (curr_l > mq.l) add(--curr_l);
            while (curr_r < mq.r) add(++curr_r);
            // Contraer rango
            while (curr_l < mq.l) remove(curr_l++);
            while (curr_r > mq.r) remove(curr_r--);

            ans[mq.id] = current_ans;
        }

        return ans;
    }
};

// ====================================================================
// SQRT DECOMPOSITION (Point Update, Range Query genérico)
// ====================================================================
class SqrtDecomposition {
private:
    int n, block_size;
    vector<ll> a, blocks;
    const ll NEUTRAL = 0; // Cambiar según operación (0 para suma, LINF para min)

    ll combine(ll x, ll y) {
        return x + y; // Cambiar a min(x, y), max(x, y), etc.
    }

public:
    /**
     * Inicializa y construye la estructura en bloques.
     * @param arr El arreglo original (0-indexed).
     */
    SqrtDecomposition(const vector<ll>& arr) {
        n = arr.size();
        block_size = max(1, (int)sqrt(n)); // Tamaño de bloque estándar
        int num_blocks = (n + block_size - 1) / block_size;

        a = arr;
        blocks.assign(num_blocks, NEUTRAL);

        // Precomputar valores por bloque
        for (int i = 0; i < n; i++) {
            // Para RMQ usar combine, para suma usar +=, etc.
            if (i % block_size == 0) blocks[i / block_size] = a[i];
            else blocks[i / block_size] = combine(blocks[i / block_size], a[i]);
        }
    }

    /**
     * Actualiza el valor en el índice 'pos' a 'val'.
     * @param pos Índice a actualizar (0-indexed).
     * @param val Nuevo valor.
     */
    void update(int pos, ll val) {
        int b = pos / block_size;
        a[pos] = val; // O a[pos] += val según el problema

        // Recomputar el bloque completo en O(sqrt(N))
        // (Para la suma simple bastaría con blocks[b] += diferencia,
        // pero recomputar todo asegura que funcione para RMQ u otras opers).
        blocks[b] = a[b * block_size];
        int end_idx = min(n, (b + 1) * block_size);
        for (int i = b * block_size + 1; i < end_idx; i++) {
            blocks[b] = combine(blocks[b], a[i]);
        }
    }

    /**
     * Consulta el rango [L, R] en O(sqrt(N)).
     * @param L Índice izquierdo (0-indexed, inclusivo).
     * @param R Índice derecho (0-indexed, inclusivo).
     * @return Resultado combinado en el rango.
     */
    ll query(int L, int R) {
        ll res = NEUTRAL;
        int b_l = L / block_size;
        int b_r = R / block_size;

        if (b_l == b_r) {
            // Están en el mismo bloque: Fuerza bruta directa
            res = a[L];
            for (int i = L + 1; i <= R; i++) res = combine(res, a[i]);
        } else {
            // 1. Cola izquierda (parcial)
            res = a[L];
            for (int i = L + 1; i < (b_l + 1) * block_size; i++) {
                res = combine(res, a[i]);
            }
            // 2. Bloques completos en el medio
            for (int b = b_l + 1; b < b_r; b++) {
                res = combine(res, blocks[b]);
            }
            // 3. Cola derecha (parcial)
            for (int i = b_r * block_size; i <= R; i++) {
                res = combine(res, a[i]);
            }
        }
        return res;
    }
};

// ====================================================================
// NORMAL TREAP (Cartesian Tree / Balanced BST)
// ====================================================================
// Generador aleatorio para las prioridades (evita anti-hash tests)
mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());

struct TreapNode {
    ll key, prior;
    int sz;
    TreapNode *l, *r;
    TreapNode(ll _key) : key(_key), prior(rng()), sz(1), l(nullptr), r(nullptr) {}
};
typedef TreapNode* pNode;

class Treap {
private:
    pNode root;

    int sz(pNode t) { return t ? t->sz : 0; }

    void upd(pNode t) {
        if (t) t->sz = 1 + sz(t->l) + sz(t->r);
    }

    // Divide el treap 't' en 'l' (keys <= key) y 'r' (keys > key)
    void split(pNode t, ll key, pNode &l, pNode &r) {
        if (!t) return void(l = r = nullptr);
        if (key >= t->key) {
            split(t->r, key, t->r, r);
            l = t;
        } else {
            split(t->l, key, l, t->l);
            r = t;
        }
        upd(t);
    }

    // Une los treaps 'l' y 'r' asumiendo que max_key(l) <= min_key(r)
    void merge(pNode &t, pNode l, pNode r) {
        if (!l || !r) t = l ? l : r;
        else if (l->prior > r->prior) {
            merge(l->r, l->r, r);
            t = l;
        } else {
            merge(r->l, l, r->l);
            t = r;
        }
        upd(t);
    }

public:
    Treap() : root(nullptr) {}

    void insert(ll key) {
        pNode l, r;
        split(root, key, l, r);
        pNode node = new TreapNode(key);
        merge(l, l, node);
        merge(root, l, r);
    }

    void erase(ll key) {
        pNode l, m, r;
        split(root, key - 1, l, r);
        split(r, key, m, r);
        if (m) {
            pNode temp = m;
            merge(m, m->l, m->r); // Borra solo 1 instancia si hay duplicados
            delete temp;
        }
        merge(root, l, m);
        merge(root, root, r);
    }

    // K-ésimo elemento más pequeño (0-indexed)
    ll kth_element(int k) {
        pNode curr = root;
        while (curr) {
            int left_size = sz(curr->l);
            if (left_size == k) return curr->key;
            if (left_size > k) curr = curr->l;
            else {
                k -= left_size + 1;
                curr = curr->r;
            }
        }
        return LINF; // Fuera de rango
    }

    int count_less_than(ll key) {
        pNode l, r;
        split(root, key - 1, l, r);
        int res = sz(l);
        merge(root, l, r);
        return res;
    }
};

// ====================================================================
// IMPLICIT TREAP (Dynamic Array with Range Queries & Reversals)
// ====================================================================
// Usa el mismo 'rng' global definido en el Treap normal.

struct ImplicitNode {
    ll val, sum, prior;
    int sz;
    bool lazy_rev;
    ImplicitNode *l, *r;

    ImplicitNode(ll _val) : val(_val), sum(_val), prior(rng()), sz(1),
                            lazy_rev(false), l(nullptr), r(nullptr) {}
};
typedef ImplicitNode* pINode;

class ImplicitTreap {
private:
    pINode root;

    int sz(pINode t) { return t ? t->sz : 0; }
    ll sum(pINode t) { return t ? t->sum : 0; }

    void push(pINode t) {
        if (t && t->lazy_rev) {
            t->lazy_rev = false;
            swap(t->l, t->r);
            if (t->l) t->l->lazy_rev ^= true;
            if (t->r) t->r->lazy_rev ^= true;
        }
    }

    void upd(pINode t) {
        if (t) {
            t->sz = 1 + sz(t->l) + sz(t->r);
            t->sum = t->val + sum(t->l) + sum(t->r); // Cambiar para Min/Max
        }
    }

    // Divide el treap dejando 'k' elementos en 'l' y el resto en 'r'
    void split(pINode t, int k, pINode &l, pINode &r) {
        if (!t) return void(l = r = nullptr);
        push(t);
        int implicit_key = sz(t->l) + 1;

        if (k >= implicit_key) {
            split(t->r, k - implicit_key, t->r, r);
            l = t;
        } else {
            split(t->l, k, l, t->l);
            r = t;
        }
        upd(t);
    }

    void merge(pINode &t, pINode l, pINode r) {
        push(l); push(r);
        if (!l || !r) t = l ? l : r;
        else if (l->prior > r->prior) {
            merge(l->r, l->r, r);
            t = l;
        } else {
            merge(r->l, l, r->l);
            t = r;
        }
        upd(t);
    }

public:
    ImplicitTreap() : root(nullptr) {}

    // Inserta 'val' en la posición 'pos' (0-indexed)
    void insert(int pos, ll val) {
        pINode l, r; // Cambio crítico: usar pINode en lugar de pNode
        split(root, pos, l, r);
        pINode node = new ImplicitNode(val);
        merge(l, l, node);
        merge(root, l, r);
    }

    // Elimina el elemento en la posición 'pos' (0-indexed)
    void erase(int pos) {
        pINode l, mid, r;
        split(root, pos, l, r);
        split(r, 1, mid, r);
        delete mid;
        merge(root, l, r);
    }

    // Revierte (Reverse) el segmento [L, R] (0-indexed)
    void reverse_range(int L, int R) {
        pINode l, mid, r;
        split(root, L, l, mid);
        split(mid, R - L + 1, mid, r);

        mid->lazy_rev ^= true;

        merge(root, l, mid);
        merge(root, root, r);
    }

    // Consulta la suma en el segmento [L, R] (0-indexed)
    ll query_range(int L, int R) {
        pINode l, mid, r;
        split(root, L, l, mid);
        split(mid, R - L + 1, mid, r);

        ll ans = sum(mid);

        merge(root, l, mid);
        merge(root, root, r);
        return ans;
    }
};


int main() {

    // Arreglo de prueba (0-indexed)
    // vector<ll> arr = {5, 2, 8, 1, 9, 3, 7, 4};
    // int n = arr.size();

    // cout << "--- 1. SPARSE TABLE (Minimo) ---\n";
    // SparseTable st(arr);
    // Mínimo en el rango [1, 4] -> {2, 8, 1, 9} => 1
    // cout << "Minimo [1, 4]: " << st.query(1, 4) << "\n\n";

    // cout << "--- 2. FENWICK TREE 1D (1-indexed) ---\n";
    // FenwickTree bit(n);
    // for (int i = 0; i < n; i++) {
        // bit.add(i + 1, arr[i]); // Pasar a 1-indexed
    // }
    // Suma en [2, 5] -> arr[1] + arr[2] + arr[3] + arr[4] = 2 + 8 + 1 + 9 = 20
    // cout << "Suma [2, 5]: " << bit.query(2, 5) << "\n\n";

    // cout << "--- 3. LAZY SEGMENT TREE (0-indexed) ---\n";
    // LazySegTree lazy_st(arr);
    // Rango original [1, 3]: {2, 8, 1} -> Suma = 11
    // cout << "Suma original [1, 3]: " << lazy_st.query(1, 3) << "\n";
    // lazy_st.add_range(1, 3, 10); // Ahora: {12, 18, 11} -> Suma = 41
    // cout << "Suma despues de sumar 10 al rango [1, 3]: " << lazy_st.query(1, 3) << "\n";
    // lazy_st.set_range(1, 2, 5);  // Ahora: {5, 5, 11} -> Suma = 21
    // cout << "Suma despues de asignar 5 al rango [1, 2]: " << lazy_st.query(1, 3) << "\n\n";

    // cout << "--- 4. MERGE SORT TREE ---\n";
    // MergeSortTree mst(arr);
    // ¿Cuántos elementos en el rango [0, n-1] son estrictamente menores a 5?
    // Elementos < 5: {2, 1, 3, 4} => 4
    // cout << "Elementos menores a 5 en todo el arreglo: " << mst.count_less_than(0, n - 1, 5) << "\n\n";

    // cout << "--- 5. FENWICK TREE 2D (1-indexed) ---\n";
    // FenwickTree2D bit2d(3, 3); // Matriz 3x3
    // bit2d.add(2, 2, 10);
    // bit2d.add(3, 1, 5);
    // Consulta submatriz desde (1,1) hasta (3,3)
    // cout << "Suma submatriz [1..3][1..3]: " << bit2d.query(1, 1, 3, 3) << "\n\n";

    // cout << "--- 6. IMPLICIT TREAP (Reverse & Range Query) ---\n";
    // ImplicitTreap treap;
    // for (int i = 0; i < n; i++) {
        // treap.insert(i, arr[i]); // Inserta en orden construyendo el arreglo
    // }
    // Suma en [0, 2] -> {5, 2, 8} => 15
    // cout << "Suma original [0, 2]: " << treap.query_range(0, 2) << "\n";
    // Revertir todo el arreglo
    // treap.reverse_range(0, n - 1);
    // Ahora el rango [0, 2] tiene los ultimos 3 elementos originales invertidos: {4, 7, 3} => 14
    // cout << "Suma en [0, 2] despues de revertir todo el arreglo: " << treap.query_range(0, 2) << "\n";

    return 0;
}
