/**
 * @file GeometryTemplateICPC.cpp
 * @brief Plantilla Integral de Geometría Computacional (ICPC / CSES)
 * @details Soporta operaciones con enteros (evita problemas de precisión),
 * Intersección de Segmentos, Área de Polígonos, Point in Polygon,
 * Puntos del Retículo (Pick's Theorem), Convex Hull, Distancia Mínima (Closest Pair),
 * y estructuras para Sweep Line 1D y 2D.
 */
//   ____ ___  ____  _____   ____  _   _
//  / ___/ _ \|  _ \| ____| / ___|| | | |
// | |  | | | | | | |  _|   \___ \| | | |
// | |__| |_| | |_| | |___   ___) | |_| |
//  \____\___/|____/|_____| |____/ \___/
//
//                  GEOMETRY TEMPLATE (CSES Supported)

#include <bits/stdc++.h>
using namespace std;
#define fastio ios::sync_with_stdio(0); cin.tie(0); cout.tie(0);

typedef long long ll;
typedef long double ld;
const ld EPS = 1e-9;

// Función segura para comparar punto flotante
int cmp(ld x, ld y = 0) {
    return (x <= y + EPS) ? (x + EPS < y ? -1 : 0) : 1;
}

// ====================================================================
// 1. ESTRUCTURAS BÁSICAS (Point)
// ====================================================================
template<typename T>
struct Point {
    T x, y;
    Point() : x(0), y(0) {}
    Point(T x, T y) : x(x), y(y) {}

    Point operator+(const Point& p) const { return Point(x + p.x, y + p.y); }
    Point operator-(const Point& p) const { return Point(x - p.x, y - p.y); }
    Point operator*(T c) const { return Point(x * c, y * c); }
    Point operator/(T c) const { return Point(x / c, y / c); }

    bool operator==(const Point& p) const { return x == p.x && y == p.y; }
    bool operator!=(const Point& p) const { return !(*this == p); }
    // Ordenamiento lexicográfico (Sweep Line, Convex Hull)
    bool operator<(const Point& p) const { return x < p.x || (x == p.x && y < p.y); }
};

template<typename T>
T dot(Point<T> a, Point<T> b) { return a.x * b.x + a.y * b.y; }

template<typename T>
T cross(Point<T> a, Point<T> b) { return a.x * b.y - a.y * b.x; }

// Producto cruz entre vectores AP y BP.
// >0 antihorario (izquierda), <0 horario (derecha), 0 colineales.
template<typename T>
T cross(Point<T> p, Point<T> a, Point<T> b) { return cross(a - p, b - p); }

template<typename T>
T norm2(Point<T> p) { return dot(p, p); }

template<typename T>
long double norm(Point<T> p) { return sqrtl(norm2(p)); }

template<typename T>
int sgn(T val) { return (T(0) < val) - (val < T(0)); }

// ====================================================================
// 2. INTERSECCIÓN DE SEGMENTOS Y PUNTO-SEGMENTO
// ====================================================================
template<typename T>
int orient(Point<T> a, Point<T> b, Point<T> c) { return sgn(cross(a, b, c)); }

template<typename T>
bool on_segment(Point<T> p, Point<T> a, Point<T> b) {
    return orient(a, b, p) == 0 &&
           min(a.x, b.x) <= p.x && p.x <= max(a.x, b.x) &&
           min(a.y, b.y) <= p.y && p.y <= max(a.y, b.y);
}

// Verifica si los segmentos AB y CD se intersectan (incluyendo colinealidad)
template<typename T>
bool segment_intersect(Point<T> a, Point<T> b, Point<T> c, Point<T> d) {
    int o1 = orient(a, b, c), o2 = orient(a, b, d);
    int o3 = orient(c, d, a), o4 = orient(c, d, b);

    // Intersección propia
    if (o1 != o2 && o3 != o4) return true;

    // Casos colineales y extremos
    if (o1 == 0 && on_segment(c, a, b)) return true;
    if (o2 == 0 && on_segment(d, a, b)) return true;
    if (o3 == 0 && on_segment(a, c, d)) return true;
    if (o4 == 0 && on_segment(b, c, d)) return true;

    return false;
}

// ====================================================================
// 3. POLÍGONOS (Área, Puntos en el Retículo, Punto dentro de Polígono)
// ====================================================================
// Retorna el DOBLE del área (Shoelace formula). Útil para evitar float.
template<typename T>
T polygon_area_2(const vector<Point<T>>& p) {
    T area = 0;
    int n = p.size();
    for (int i = 0; i < n; i++) {
        area += cross(p[i], p[(i + 1) % n]);
    }
    return abs(area);
}

// Pick's Theorem: Área = Interiores + Borde/2 - 1
// -> 2*Área = 2*Interiores + Borde - 2
// Esta función calcula los puntos del límite (borde)
long long boundary_points(const vector<Point<long long>>& poly) {
    long long b = 0;
    int n = poly.size();
    for (int i = 0; i < n; i++) {
        Point<long long> p1 = poly[i], p2 = poly[(i + 1) % n];
        b += std::gcd(abs(p1.x - p2.x), abs(p1.y - p2.y));
    }
    return b;
}

// Verifica si el punto p está dentro de un polígono.
// Retorna: 0 = afuera, 1 = en el borde, 2 = adentro.
template<typename T>
int point_in_polygon(const vector<Point<T>>& poly, Point<T> p) {
    int n = poly.size();
    int hits = 0;
    for (int i = 0; i < n; i++) {
        Point<T> a = poly[i], b = poly[(i + 1) % n];
        if (on_segment(p, a, b)) return 1;
        if (a.y > b.y) swap(a, b);
        if (p.y >= a.y && p.y < b.y && cross(p, a, b) > 0) hits++;
    }
    return (hits % 2 != 0) ? 2 : 0;
}

// ====================================================================
// 4. CONVEX HULL (Monotone Chain)
// ====================================================================
// Retorna la envoltura convexa del conjunto de puntos (en orden antihorario)
template<typename T>
vector<Point<T>> convex_hull(vector<Point<T>> pts) {
    int n = pts.size(), k = 0;
    if (n <= 2) return pts;

    vector<Point<T>> h(2 * n);
    sort(pts.begin(), pts.end());

    // Mitad inferior
    for (int i = 0; i < n; i++) {
        // Cambiar <= a < si se quieren incluir los puntos colineales en el borde
        while (k >= 2 && cross(h[k - 2], h[k - 1], pts[i]) <= 0) k--;
        h[k++] = pts[i];
    }
    // Mitad superior
    for (int i = n - 2, t = k + 1; i >= 0; i--) {
        while (k >= t && cross(h[k - 2], h[k - 1], pts[i]) <= 0) k--;
        h[k++] = pts[i];
    }

    h.resize(k - 1);
    return h;
}

// ====================================================================
// 5. DISTANCIA MÍNIMA (CLOSEST PAIR) - O(N log N)
// ====================================================================
// Retorna la distancia euclidiana AL CUADRADO entre los 2 puntos más cercanos.
long long closest_pair(vector<Point<long long>> pts) {
    int n = pts.size();
    if (n < 2) return 8e18; // Evitar RTE

    sort(pts.begin(), pts.end());
    set<pair<long long, long long>> active; // {y, x}
    long long best_dist2 = 8e18; // Infinito seguro
    int left = 0;

    for (int i = 0; i < n; i++) {
        long long d = ceil(sqrt(best_dist2));
        while (left < i && pts[i].x - pts[left].x >= d) {
            active.erase({pts[left].y, pts[left].x});
            left++;
        }
        auto it1 = active.lower_bound({pts[i].y - d, -8e18});
        auto it2 = active.upper_bound({pts[i].y + d, 8e18});
        for (auto it = it1; it != it2; ++it) {
            long long dy = pts[i].y - it->first;
            long long dx = pts[i].x - it->second;
            best_dist2 = min(best_dist2, dx * dx + dy * dy);
        }
        active.insert({pts[i].y, pts[i].x});
    }
    return best_dist2;
}

// ====================================================================
// 6. ESTRUCTURAS SWEEP LINE (1D y 2D)
// ====================================================================
struct Event1D {
    long long x;
    int type; // +1 Inicio, -1 Fin
    int id;

    bool operator<(const Event1D& o) const {
        if (x != o.x) return x < o.x;
        return type > o.type; // Ajustar (type > o.type) = procesa Entradas primero
    }
};

struct Event2D {
    long long x;
    int type; // +1 Borde Izquierdo, -1 Borde Derecho
    long long y1, y2;

    bool operator<(const Event2D& o) const {
        if (x != o.x) return x < o.x;
        return type > o.type;
    }
};


// ====================================================================
// 7. ORDENAMIENTO POLAR (O(N log N) puramente entero)
// ====================================================================
// Ordena los puntos en sentido antihorario alrededor de 'center'.
// Empieza en el eje X positivo. Cero floats involucrados.
template<typename T>
void polar_sort(vector<Point<T>>& pts, Point<T> center = Point<T>(0, 0)) {
    auto half = [](Point<T> p) {
        return p.y > 0 || (p.y == 0 && p.x > 0) ? 1 : -1;
    };

    sort(pts.begin(), pts.end(), [&](Point<T> a, Point<T> b) {
        Point<T> va = a - center, vb = b - center;
        int ha = half(va), hb = half(vb);
        if (ha != hb) return ha > hb; // ¡CAMBIO AQUÍ! '>' pone la mitad superior (1) primero
        T cr = cross(va, vb);
        if (cr != 0) return cr > 0;
        return norm2(va) < norm2(vb);
    });
}

// ====================================================================
// 8. LÍNEAS E INTERSECCIONES EXACTAS (Requiere 'ld')
// ====================================================================
template<typename T>
ld dist_to_line(Point<T> p, Point<T> a, Point<T> b) {
    return abs((ld)cross(p, a, b)) / norm(a - b);
}

template<typename T>
ld dist_to_segment(Point<T> p, Point<T> a, Point<T> b) {
    if (a == b) return norm(p - a);
    if (dot(p - a, b - a) <= 0) return norm(p - a);
    if (dot(p - b, a - b) <= 0) return norm(p - b);
    return dist_to_line(p, a, b);
}

// Intersección exacta de dos líneas infinitas AB y CD
// Verifica antes que no sean paralelas: cross(b-a, d-c) != 0
Point<ld> line_intersect(Point<ld> a, Point<ld> b, Point<ld> c, Point<ld> d) {
    ld cr1 = cross(b - a, c - a);
    ld cr2 = cross(b - a, d - a);
    return c + (d - c) * (cr1 / (cr1 - cr2));
}

// ====================================================================
// 9. ROTATING CALIPERS (Diámetro de un Polígono) - O(N)
// ====================================================================
// Devuelve la máxima distancia al CUADRADO entre cualquier par de puntos.
// 'poly' DEBE ser un Convex Hull en orden antihorario.
template<typename T>
T rotating_calipers(const vector<Point<T>>& poly) {
    int n = poly.size();
    if (n < 2) return 0;
    if (n == 2) return norm2(poly[0] - poly[1]);

    T max_dist2 = 0;
    int j = 1;

    for (int i = 0; i < n; i++) {
        int nxt_i = (i + 1) % n;
        while (abs(cross(poly[j], poly[i], poly[nxt_i])) <
               abs(cross(poly[(j + 1) % n], poly[i], poly[nxt_i]))) {
            j = (j + 1) % n;
        }
        max_dist2 = max({max_dist2, norm2(poly[i] - poly[j]), norm2(poly[nxt_i] - poly[j])});
    }
    return max_dist2;
}

// ====================================================================
// 10. CÍRCULOS E INTERSECCIONES
// ====================================================================
struct Circle {
    Point<ld> c;
    ld r;
    Circle() : c(Point<ld>(0,0)), r(0) {}
    Circle(Point<ld> c, ld r) : c(c), r(r) {}
};

// Intersección de Línea infinita AB con Círculo
vector<Point<ld>> circle_line_intersect(Point<ld> a, Point<ld> b, Circle circ) {
    Point<ld> dir = b - a;
    ld len2 = norm2(dir);
    if (cmp(len2) == 0) return {}; // A y B son el mismo punto

    ld t = dot(circ.c - a, dir) / len2;
    Point<ld> proj = a + dir * t; // Proyección del centro a la línea

    ld dist2 = norm2(circ.c - proj);
    if (cmp(dist2, circ.r * circ.r) > 0) return {}; // No intersecta
    if (cmp(dist2, circ.r * circ.r) == 0) return {proj}; // Tangente (1 punto)

    ld offset = sqrtl(max((ld)0.0, circ.r * circ.r - dist2) / len2);
    return {proj - dir * offset, proj + dir * offset}; // Secante (2 puntos)
}

// Intersección de dos Círculos
vector<Point<ld>> circle_circle_intersect(Circle c1, Circle c2) {
    Point<ld> dir = c2.c - c1.c;
    ld d = norm(dir);

    // Separados, uno dentro de otro, o concéntricos
    if (cmp(d, c1.r + c2.r) > 0 || cmp(d, abs(c1.r - c2.r)) < 0 || cmp(d) == 0) return {};

    ld a = (c1.r * c1.r - c2.r * c2.r + d * d) / (2 * d);
    ld h = sqrtl(max((ld)0.0, c1.r * c1.r - a * a));

    Point<ld> p2 = c1.c + dir * (a / d);
    Point<ld> perp(-dir.y, dir.x);
    perp = perp * (h / d);

    if (cmp(h) == 0) return {p2}; // Tangentes exterior/interiormente
    return {p2 + perp, p2 - perp}; // Se cortan en 2 puntos
}



// Nota para Union de Rectángulos / Area: requiere un Segment Tree
// que maneje Range Sum (Active Length) con la compresión de coordenadas Y.

int main() {
    fastio

    // ==========================================
    // EJEMPLOS MÁGICOS DE GEOMETRÍA
    // ==========================================

    // 1. ÁREA DE UN POLÍGONO (Shoelace)
    // cout << "--- 1. AREA DE POLIGONO ---\n";
    // vector<Point<ll>> poly = {{0,0}, {4,0}, {4,3}, {0,3}}; // Rectángulo 4x3
    // El método devuelve el DOBLE del área para evitar floats
    // cout << "Area real: " << polygon_area_2(poly) / 2.0 << "\n\n"; // 12

    // 2. PUNTO EN POLÍGONO (Point in Polygon)
    // cout << "--- 2. PUNTO EN POLIGONO ---\n";
    // Point<ll> p_in(2, 1), p_out(5, 5), p_edge(4, 1);
    // 0 = Afuera, 1 = Borde, 2 = Adentro
    // cout << "Adentro: " << point_in_polygon(poly, p_in) << "\n";
    // cout << "Borde: " << point_in_polygon(poly, p_edge) << "\n";
    // cout << "Afuera: " << point_in_polygon(poly, p_out) << "\n\n";

    // 3. CONVEX HULL & ROTATING CALIPERS
    // cout << "--- 3. CONVEX HULL & MAX DISTANCE ---\n";
    // vector<Point<ll>> scatter = {{0,0}, {1,1}, {2,2}, {3,1}, {4,0}, {2,-2}};
    // vector<Point<ll>> hull = convex_hull(scatter);
    // cout << "Puntos en el Hull: " << hull.size() << "\n";
    // Calipers requiere la envolvente convexa limpia
    // ll max_dist_sq = rotating_calipers(hull);
    // cout << "Distancia maxima al cuadrado: " << max_dist_sq << "\n\n";

    // 4. DISTANCIA MÍNIMA (Closest Pair)
    // cout << "--- 4. CLOSEST PAIR ---\n";
    // vector<Point<ll>> cities = {{0,0}, {100,100}, {1000,1000}, {102, 102}};
    // Debe detectar que (100,100) y (102,102) son los más cercanos
    // cout << "Distancia minima al cuadrado: " << closest_pair(cities) << "\n\n";

    // 5. ORDENAMIENTO POLAR
    // cout << "--- 5. POLAR SORT ---\n";
    // vector<Point<ll>> angles = {{0,-1}, {-1,0}, {0,1}, {1,0}};
    // polar_sort(angles, Point<ll>(0,0));
    // cout << "Orden Antihorario: ";
    // for(auto p : angles) cout << "(" << p.x << "," << p.y << ") ";
    // Imprime: (1,0) (0,1) (-1,0) (0,-1)
    // cout << "\n\n";

    // 6. INTERSECCIÓN EXACTA DE LÍNEAS (Requiere 'ld')
    // cout << "--- 6. LINE INTERSECT ---\n";
    // Point<ld> l1_a(0,0), l1_b(10,10); // Línea y = x
    // Point<ld> l2_a(0,10), l2_b(10,0); // Línea y = -x + 10
    // Point<ld> inters = line_intersect(l1_a, l1_b, l2_a, l2_b);
    // cout << "Cruce en: (" << inters.x << ", " << inters.y << ")\n\n"; // (5, 5)

    // 7. CÍRCULOS
    // cout << "--- 7. CIRCLE INTERSECTIONS ---\n";
    // Circle c1(Point<ld>(0,0), 5); // Círculo radio 5 origen
    // Circle c2(Point<ld>(8,0), 5); // Círculo radio 5 en X=8

    // vector<Point<ld>> c_inter = circle_circle_intersect(c1, c2);
    // cout << "Puntos de corte entre circulos: " << c_inter.size() << "\n";
    // for(auto p : c_inter) cout << "Cortan en: (" << p.x << ", " << p.y << ")\n";

    // Línea horizontal y=3 cortando el círculo c1
    // Point<ld> la(0, 3), lb(10, 3);
    // vector<Point<ld>> l_inter = circle_line_intersect(la, lb, c1);
    // cout << "La linea corta el circulo en " << l_inter.size() << " puntos.\n";

    return 0;
}
