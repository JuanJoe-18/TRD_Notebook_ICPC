# Guía Maestra TRD: Uso de Plantillas, Modificabilidad y Reconocimiento de Problemas

Esta guía está diseñada como documento de referencia y consulta rápida para el equipo durante competencias ICPC. Aborda tres pilares esenciales:
1. **Consistencia de Índices (0-based vs 1-based)** en la biblioteca del TRD.
2. **Modificabilidad Avanzada**: Cómo transformar plantillas existentes paso a paso (ej. caminos, estados extendidos, lazy, aristas en árboles).
3. **Reconocimiento de Problemas Complejos e Implícitos**: Cómo detectar la técnica cuando el enunciado oculta el grafo, el flujo o la estructura.
4. **Técnicas sobre Grillas**: Aplanamiento, bipartición y restricciones de celdas.

---

## 1. Regla de Oro de Índices: 0-indexed vs 1-indexed

Uno de los errores más costosos y comunes en ICPC es mezclar estructuras con indexación incompatible (desbordamientos `out-of-bounds` o bucles infinitos). 

### Tabla de Índices de la Biblioteca TRD

| Estructura / Algoritmo | Archivo | Base | Motivo / Cuidado Crítico |
| :--- | :--- | :---: | :--- |
| **`Graph<T>`** | `graph.hpp` | **1-based** | `adj` tiene tamaño `n + 1`. Admite nodos `1 .. n`. |
| **`Dijkstra<T>`** | `graph.hpp` | **1-based** | Devuelve vector de tamaño `g.n + 1`. Nodo `src` debe ser $\ge 1$. |
| **`ZeroOneBfs<T>`** | `graph.hpp` | **1-based** | Compatible con `Graph`. |
| **`DSU`** | `graph.hpp` | **0-based** | `iota(all(p), 0)`. Nodos `0 .. n - 1`. Si el problema es 1-based: inicializar con `n + 1` o restar `1`. |
| **`FenwickTree` (BIT)** | `rq.hpp` | **1-based** | **CRÍTICO:** `idx += idx & -idx`. Si consultas o insertas en `0`, se produce un **bucle infinito**. |
| **`SparseTable`** | `rq.hpp` | **0-based** | Rangos `[L, R]` con $0 \le L \le R < n$. |
| **`IterativeSegTree`** | `rq.hpp` | **0-based** | Rangos semiabiertos o cerrados sobre índices `0 .. n - 1`. |
| **`HLD`** | `tree.hpp` | **1-based** | Nodos del árbol `1 .. n`. Raíz por defecto en `1`. |
| **`LCA` / Binlifting** | `tree.hpp` | **1-based** | Nodos `1 .. n`. `up[node][0] = parent`. |

### Patrón de Normalización Inmediata al Leer
Para no dudar durante el contest, define una convención mental al leer la entrada:
```cpp
// Si el problema da aristas 1-based y tu estructura es 0-based:
int u, v; cin >> u >> v;
u--; v--; // Normalizado a 0-based

// Si usas Graph/HLD/Fenwick de la plantilla (que son 1-based):
int u, v; cin >> u >> v;
g.add_undirected_edge(u, v, w); // Directo, sin restar 1
```

---

## 2. Modificabilidad Avanzada de Funciones

### A. De Dijkstra Estándar a Dijkstra con Camino (`DijkstraPath`)
En `graph.hpp`, la plantilla ya cuenta con `run_with_parents` y `restore_path`:

```cpp
auto [dist, parent] = Dijkstra<ll>().run_with_parents(g, start_node);

if (dist[end_node] == LINF) {
    cout << -1 << "\n"; // Inalcanzable
} else {
    vi path = Dijkstra<ll>::restore_path(end_node, parent);
    // path contiene: {start_node, ..., end_node}
    for (int node : path) cout << node << " ";
    cout << "\n";
}
```

### B. Dijkstra en Grafo de Estados / Multi-Estado (Estado Extendido)
**Problema típico:** "Halla el camino mínimo, pero puedes usar hasta $K$ descuentos a mitad de precio" o "cada paso gasta combustible".

**¿Cómo modificar Dijkstra sin reescribir todo?:**
En lugar de un nodo simple $u$, el estado es $(u, k)$, donde $u$ es el vértice actual y $k$ es el recurso restante ($0 \le k \le K$).

```cpp
// 1. Matriz de distancias multidimensional
vector<vll> dist(n + 1, vll(K + 1, LINF));

// 2. Priority queue guarda: {costo_acumulado, {u, estado_k}}
using State = pair<ll, pair<int, int>>;
priority_queue<State, vector<State>, greater<State>> pq;

dist[src][K] = 0;
pq.push({0, {src, K}});

while (!pq.empty()) {
    auto [d, state] = pq.top(); pq.pop();
    auto [u, rem] = state;
    if (d > dist[u][rem]) continue;

    for (int id : g.adj[u]) {
        auto& e = g.edges[id];
        
        // Transición 1: Avanzar normalmente
        if (dist[u][rem] + e.weight < dist[e.to][rem]) {
            dist[e.to][rem] = dist[u][rem] + e.weight;
            pq.push({dist[e.to][rem], {e.to, rem}});
        }
        
        // Transición 2: Usar descuento/habilidad (si queda recurso)
        if (rem > 0) {
            ll cost_with_skill = dist[u][rem] + e.weight / 2;
            if (cost_with_skill < dist[e.to][rem - 1]) {
                dist[e.to][rem - 1] = cost_with_skill;
                pq.push({dist[e.to][rem - 1], {e.to, rem - 1}});
            }
        }
    }
}
```

### C. Segment Tree: Desacoplamiento de `Node` y Operación `merge`
Para evitar rehacer el árbol cuando la consulta no es una suma simple (ej. consultar el subsegmento de suma máxima, contar ceros, o mantener pares ordenados):

```cpp
struct Node {
    ll sum, pref, suff, ans;
    Node() : sum(0), pref(0), suff(0), ans(0) {} // Valor neutro
    Node(ll val) {
        sum = val;
        pref = suff = ans = max(0LL, val);
    }
};

Node merge(const Node& a, const Node& b) {
    Node res;
    res.sum = a.sum + b.sum;
    res.pref = max(a.pref, a.sum + b.pref);
    res.suff = max(b.suff, b.sum + a.suff);
    res.ans = max({a.ans, b.ans, a.suff + b.pref});
    return res;
}
```
*Regla:* La lógica del Segment Tree (`update`, `query`) solo invoca `merge(left_child, right_child)`. Nunca toques los índices del árbol; solo redefine `struct Node` y `merge`.

### D. HLD: Modificación para Valores en las Aristas (Edge Weights)
HLD por defecto coloca los valores en los **vértices**. Si el problema tiene pesos en las **aristas**:
1. Empuja el peso de cada arista $(u, v)$ al vértice más profundo:
   ```cpp
   int child = (depth[u] > depth[v] ? u : v);
   // Asociar peso de la arista al nodo 'child'
   ```
2. Al consultar el camino entre $u$ y $v$, procesa las cadenas como siempre, pero **excluye el LCA**:
   ```cpp
   // Cuando u y v quedan en la misma cadena pesada (asumiendo pos[u] <= pos[v]):
   // En lugar de query(pos[u], pos[v]):
   if (u != v) {
       res = merge(res, segtree.query(pos[u] + 1, pos[v])); 
   }
   ```
   *(El `+1` evita incluir el LCA, ya que el peso asociado al LCA pertenece a la arista hacia su padre, no al camino entre $u$ y $v$).*

---

## 3. Problemas sobre Grillas ($R \times C$) y Extrapolaciones

### A. Aplanamiento de 2D a 1D
Muchos algoritmos (Dinic, DSU, Dijkstra) trabajan mejor con nodos enumerados linealmente de $0$ a $R \times C - 1$:
```cpp
auto get_id = [&](int r, int c) { return r * C + c; };
auto get_coords = [&](int id) { return make_pair(id / C, id % C); };
```

### B. Movimiento y Validación
Usa siempre el bloque provisto en `graph.hpp`:
```cpp
for (int d = 0; d < 4; d++) {
    int nr = r + dx4[d];
    int nc = c + dy4[d];
    if (isValid(nr, nc, R, C) && grid[nr][nc] != '#') {
        // Conectar o encolar
    }
}
```

### C. La Grilla como Grafo Bipartito (Coloración de Tablero)
**Propiedad matemática:** Toda grilla con movimientos 4-direccionales es **bipartita**.
- Celda $(r, c)$ es Blanca si $(r + c) \% 2 == 0$.
- Celda $(r, c)$ es Negra si $(r + c) \% 2 == 1$.
- Ninguna celda blanca es vecina de otra blanca.

**Cuándo usarlo:**
- Problemas de colocar fichas de dominó $1 \times 2$ sin solaparse $\to$ **Maximum Bipartite Matching**.
- Destruir el mínimo número de celdas para que no queden parejas adyacentes $\to$ **Minimum Vertex Cover** (igual al Max Matching por teorema de König).
- Red de flujo: Fuente $\to$ Celdas Blancas $\to$ Celdas Negras adyacentes $\to$ Sumidero.

### D. Restricción de Vértices en Grillas (Flujo con Capacidad en Celdas)
Si el problema dice: "Cada celda solo puede ser atravesada por un camino una sola vez":
- Divide la celda $u$ en dos nodos: $u_{in}$ y $u_{out}$.
- Agrega una arista interna: $u_{in} \to u_{out}$ con capacidad $1$ (o el límite de visitas de la celda).
- Para celdas adyacentes $u$ y $v$, conecta $u_{out} \to v_{in}$ con capacidad $\infty$.

---

## 4. Reconocimiento de Problemas Complejos e Implícitos

### A. Min-Cut Disfrazado: Selección de Proyectos (Project Selection)
**Patrón:** 
- Tienes opciones con ganancias positivas y costos asociados.
- Condición de dependencia: *"Si realizas el proyecto $A$, estás obligado a comprar la máquina $B$."*
- O bien: *"Si el elemento $X$ pertenece al conjunto $1$ y el elemento $Y$ al conjunto $2$, pagas una penalización."*

**Construcción de la Red:**
1. Fuente $S$ a cada proyecto con ganancia: arista con capacidad = $\text{Ganancia}$.
2. Cada máquina/costo al Sumidero $T$: arista con capacidad = $\text{Costo}$.
3. Aristas de dependencia $A \to B$ con capacidad = $\infty$.
4. $\text{Ganancia Máxima} = \sum (\text{Todas las Ganancias}) - \text{MinCut}(S, T)$.

### B. 2-SAT Implícito (Decisiones Binarias con Conflictos)
**Patrón:**
- $N$ variables que deben tomar exactamente uno de dos estados (ej. Interruptor ON/OFF, persona va a la Sala 1 o Sala 2, punto en posición $X$ o $Y$).
- Restricciones binarias de incompatibilidad: *"No pueden estar ambos apagados"* o *"Si seleccionas $A$, no puedes tomar $B$"*.

**Construcción de Implicaciones:**
- La restricción *"al menos uno debe ser verdadero"* ($A \lor B$) equivale a:
  - $(\neg A \implies B)$ y $(\neg B \implies A)$.
- *"No pueden ocurrir ambos"* ($\neg A \lor \neg B$):
  - $(A \implies \neg B)$ y $(B \implies \neg A)$.
- Resuelve con `TwoSat` de la plantilla. Si $A$ y $\neg A$ quedan en la misma SCC $\to$ Imposible.

### C. Grafos Implícitos y BFS 0-1
**Patrón:**
- No te dan un grafo explícito de aristas. Te dan un estado (ej. un cubo Rubik $2 \times 2$, un tablero de números, o una cadena que permuta).
- Los costos de transición solo pueden ser $0$ o $1$:
  - Ej. *"Girar a la derecha es gratis (costo 0), pero avanzar cuesta 1 segundo"*.
- **Acción:** No uses Dijkstra ($O(E \log V)$). Usa `ZeroOneBfs` con `std::deque` ($O(V + E)$). Costo 0 va a `push_front`, costo 1 a `push_back`.

### D. Meet-in-the-Middle
**Patrón:**
- $N \approx 36 \dots 44$.
- Demasiado grande para $O(2^N)$ ($\approx 10^{12}$ operaciones, TLE), pero demasiado disperso o sin subestructura óptima para DP clásico.
- **Acción:**
  1. Divide el conjunto en dos mitades de tamaño $N/2 \approx 20$.
  2. Genera todos los subconjuntos de la mitad 1 en un vector ($2^{20} \approx 10^6$ estados).
  3. Ordena el vector resultante.
  4. Para cada uno de los $2^{20}$ estados de la mitad 2, usa búsqueda binaria (`std::lower_bound`) o Two Pointers en el primer vector.

### E. Búsqueda Binaria sobre la Respuesta (BS on Answer)
**Patrón:**
- Enunciados con frases como:
  - *"Halla el valor máximo posible tal que el mínimo sea..."*
  - *"Minimiza la máxima distancia..."*
  - *"¿Cuál es el tiempo mínimo para completar la tarea?"*
- La función de verificación `check(mid)` es monótona: si es posible con $mid$, también es posible para valores mayores (o menores).
- **Acción:** La verificación suele reducirse a un Greedy o a una comprobación de caminos de costo $\le mid$.

---

## 5. Checklist de Verificación Rápida Pre-Submit

- [ ] **Límites de Overflow:** ¿Las sumas o distancias superan $2 \times 10^9$? Si sí, ¿está todo en `long long` (`ll`)? ¿Se escribió `1LL * a * b` en multiplicaciones intermedias?
- [ ] **Infinito Suficiente:** En grafos con pesos grandes, ¿es `INF` ($10^9$) suficiente o necesitas `LINF` ($10^{18}$)?
- [ ] **Casos Extremos (Corner Cases):**
  - $N = 1$ o $N = 2$.
  - Grafo no conexo o nodo inalcanzable (revisar si la respuesta debe ser `-1`).
  - Cadenas vacías o con todos los caracteres iguales.
  - Árboles que son una línea recta (degenerados).
- [ ] **Casos de Múltiples Test Cases (`T > 1`):**
  - ¿Limpiaste todas las estructuras globales? (`g.adj.clear()`, `dsu.p.clear()`, `vector<vector<int>>().swap(...)`).
  - ¿Evitaste el uso de `memset` sobre tamaños excesivos fuera de $N$?
