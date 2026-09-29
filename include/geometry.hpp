/**
 * Uso: cmp(a, b);
 * Comparacion de numeros de punto flotante con tolerancia $EPS$; devuelve -1, 0 o 1.
 * Complejidad: $O(1)$.
 */
int cmp(ld x, ld y = 0) { return x <= y + EPS ? (x + EPS < y ? -1 : 0) : 1; }

/**
 * Uso: Point<ll> p(x, y); p+q, p-q, p*c, p/c, p==q; orden lexicografico con <
 * Punto 2D generico con operadores aritmeticos, igualdad y orden lexicografico.
 * Complejidad: $O(1)$ por operacion.
 */
template<typename T> struct Point {
  T x, y;
  Point(T x = 0, T y = 0) : x(x), y(y) {}
  Point operator+(const Point& p) const { return {x + p.x, y + p.y}; }
  Point operator-(const Point& p) const { return {x - p.x, y - p.y}; }
  Point operator*(T c) const { return {x * c, y * c}; }
  Point operator/(T c) const { return {x / c, y / c}; }
  bool operator==(const Point& p) const { return x == p.x && y == p.y; }
  bool operator<(const Point& p) const { return x < p.x || (x == p.x && y < p.y); } // lexicografico
};
/**
 * Uso: dot(a, b), cross(a, b), cross(p, a, b), norm2(p), norm(p), sgn(v);
 * Operaciones vectoriales basicas: producto punto, producto cruz, norma al cuadrado, norma y signo.
 * Complejidad: $O(1)$ por operacion.
 */
template<typename T> T dot(Point<T> a, Point<T> b) { return a.x * b.x + a.y * b.y; }
template<typename T> T cross(Point<T> a, Point<T> b) { return a.x * b.y - a.y * b.x; }
template<typename T> T cross(Point<T> p, Point<T> a, Point<T> b) { return cross(a - p, b - p); }
template<typename T> T norm2(Point<T> p) { return dot(p, p); }
template<typename T> ld norm(Point<T> p) { return sqrtl(norm2(p)); }
template<typename T> int sgn(T v) { return (T(0) < v) - (v < T(0)); }
/**
 * Uso: orient(a, b, c), on_segment(p, a, b), segment_intersect(a, b, c, d);
 * Orientacion de tres puntos (horaria, colineal o antihoraria), pertenencia de un punto a un segmento e interseccion de segmentos.
 * Complejidad: $O(1)$ por operacion.
 */
template<typename T> int orient(Point<T> a, Point<T> b, Point<T> c) { return sgn(cross(a, b, c)); }

template<typename T> bool on_segment(Point<T> p, Point<T> a, Point<T> b) {
  return orient(a, b, p) == 0 && min(a.x, b.x) <= p.x && p.x <= max(a.x, b.x) &&
         min(a.y, b.y) <= p.y && p.y <= max(a.y, b.y);
}
template<typename T> bool segment_intersect(Point<T> a, Point<T> b, Point<T> c, Point<T> d) {
  int o1 = orient(a, b, c), o2 = orient(a, b, d), o3 = orient(c, d, a), o4 = orient(c, d, b);
  if (o1 != o2 && o3 != o4) return true;
  if (o1 == 0 && on_segment(c, a, b)) return true;
  if (o2 == 0 && on_segment(d, a, b)) return true;
  if (o3 == 0 && on_segment(a, c, d)) return true;
  if (o4 == 0 && on_segment(b, c, d)) return true;
  return false;
}

/**
 * Uso: polygon_area_2(poly) / 2 = area, boundary_points(poly), point_in_polygon(poly, p);
 * Area doble de un poligono (para evitar floats), puntos del reticulo en el borde (teorema de Pick) y clasificacion de un punto (0 fuera, 1 borde, 2 dentro).
 * Complejidad: $O(N)$ por operacion con $N$ vertices.
 */
template<typename T> T polygon_area_2(const vector<Point<T>>& p) {
  T area = 0;
  for (int i = 0; i < p.size(); i++) area += cross(p[i], p[(i + 1) % p.size()]);
  return abs(area);
}
// Puntos del reticulo en el borde (Pick: 2*A = 2*I + B - 2)
long long boundary_points(const vector<Point<ll>>& poly) {
  ll b = 0;
  for (int i = 0; i < poly.size(); i++) {
    auto p1 = poly[i], p2 = poly[(i + 1) % poly.size()];
    b += gcd(abs(p1.x - p2.x), abs(p1.y - p2.y));
  }
  return b;
}
// Retorna 0 afuera, 1 en el borde, 2 adentro
template<typename T> int point_in_polygon(const vector<Point<T>>& poly, Point<T> p) {
  int hits = 0;
  for (int i = 0; i < poly.size(); i++) {
    auto a = poly[i], b = poly[(i + 1) % poly.size()];
    if (on_segment(p, a, b)) return 1;
    if (a.y > b.y) swap(a, b);
    if (p.y >= a.y && p.y < b.y && cross(p, a, b) > 0) hits++;
  }
  return hits & 1 ? 2 : 0;
}

/**
 * Uso: auto hull = convex_hull(puntos);
 * Cascara convexa (monotonic chain) en sentido antihorario, sin puntos colineales en el borde.
 * Complejidad: $O(N \log N)$.
 */
template<typename T> vector<Point<T>> convex_hull(vector<Point<T>> pts) {
  int n = pts.size(), k = 0;
  if (n <= 2) return pts;
  vector<Point<T>> h(2 * n);
  sort(all(pts));
  for (int i = 0; i < n; i++) {
    while (k >= 2 && cross(h[k - 2], h[k - 1], pts[i]) <= 0) k--;
    h[k++] = pts[i];
  }
  for (int i = n - 2, t = k + 1; i >= 0; i--) {
    while (k >= t && cross(h[k - 2], h[k - 1], pts[i]) <= 0) k--;
    h[k++] = pts[i];
  }
  h.resize(k - 1);
  return h;
}

/**
 * Uso: closest_pair(puntos);
 * Distancia al cuadrado entre el par de puntos mas cercanos.
 * Complejidad: $O(N \log N)$.
 */
long long closest_pair(vector<Point<ll>> pts) {
  int n = pts.size();
  if (n < 2) return 8e18;
  sort(all(pts));
  set<pair<ll, ll>> active; // {y, x}
  ll best = 8e18;
  int left = 0;
  for (int i = 0; i < n; i++) {
    ll d = ceil(sqrtl(best));
    while (left < i && pts[i].x - pts[left].x >= d) active.erase({pts[left].y, pts[left].x}), left++;
    auto lo = active.lower_bound({pts[i].y - d, -8e18});
    auto hi = active.upper_bound({pts[i].y + d, 8e18});
    for (auto it = lo; it != hi; ++it) {
      ll dy = pts[i].y - it->fi, dx = pts[i].x - it->se;
      best = min(best, dx * dx + dy * dy);
    }
    active.insert({pts[i].y, pts[i].x});
  }
  return best;
}

/**
 * Uso: Event1D / Event2D para ordenar eventos de sweep line por coordenada x.
 * Eventos para sweep line en 1D y 2D; type +1 entrada/inicio, -1 salida/fin.
 * Complejidad: $O(1)$ por evento, $O(K \log K)$ al ordenarlos con $K$ eventos.
 */
struct Event1D {
  ll x;
  int type, id;
  bool operator<(const Event1D& o) const { return x != o.x ? x < o.x : type > o.type; }
};
struct Event2D {
  ll x, y1, y2;
  int type; // +1 borde izq, -1 borde der
  bool operator<(const Event2D& o) const { return x != o.x ? x < o.x : type > o.type; }
};

/**
 * Uso: polar_sort(puntos, centro);
 * Ordena puntos por angulo polar alrededor de un centro, en sentido antihorario desde el eje $+x$.
 * Complejidad: $O(N \log N)$.
 */
template<typename T> void polar_sort(vector<Point<T>>& pts, Point<T> center = {}) {
  auto half = [](Point<T> p) { return p.y > 0 || (p.y == 0 && p.x > 0) ? 1 : -1; };
  sort(all(pts), [&](Point<T> a, Point<T> b) {
    auto va = a - center, vb = b - center;
    int ha = half(va), hb = half(vb);
    if (ha != hb) return ha > hb;
    T cr = cross(va, vb);
    if (cr != 0) return cr > 0;
    return norm2(va) < norm2(vb);
  });
}

/**
 * Uso: dist_to_line(p, a, b), dist_to_segment(p, a, b), line_intersect(a, b, c, d);
 * Distancias de un punto a una linea o segmento e interseccion exacta de dos lineas infinitas no paralelas.
 * Complejidad: $O(1)$ por operacion.
 */
template<typename T> ld dist_to_line(Point<T> p, Point<T> a, Point<T> b) {
  return abs((ld)cross(p, a, b)) / norm(a - b);
}
template<typename T> ld dist_to_segment(Point<T> p, Point<T> a, Point<T> b) {
  if (a == b) return norm(p - a);
  if (dot(p - a, b - a) <= 0) return norm(p - a);
  if (dot(p - b, a - b) <= 0) return norm(p - b);
  return dist_to_line(p, a, b);
}
// Interseccion exacta de lineas infinitas AB y CD (no paralelas)
Point<ld> line_intersect(Point<ld> a, Point<ld> b, Point<ld> c, Point<ld> d) {
  ld cr1 = cross(b - a, c - a), cr2 = cross(b - a, d - a);
  return c + (d - c) * (cr1 / (cr1 - cr2));
}

/**
 * Uso: rotating_calipers(hull);
 * Diametro al cuadrado de un poligono convexo mediante calipers rotativos.
 * Complejidad: $O(N)$.
 */
template<typename T> T rotating_calipers(const vector<Point<T>>& poly) {
  int n = poly.size();
  if (n < 2) return 0;
  if (n == 2) return norm2(poly[0] - poly[1]);
  T best = 0;
  int j = 1;
  for (int i = 0; i < n; i++) {
    int ni = (i + 1) % n;
    while (abs(cross(poly[j], poly[i], poly[ni])) < abs(cross(poly[(j + 1) % n], poly[i], poly[ni]))) j = (j + 1) % n;
    best = max({best, norm2(poly[i] - poly[j]), norm2(poly[ni] - poly[j])});
  }
  return best;
}

/**
 * Uso: Circle c(centro, radio); circle_line_intersect(a, b, c); circle_circle_intersect(c1, c2);
 * Circulo e intersecciones con lineas y otros circulos; devuelven 0, 1 o 2 puntos.
 * Complejidad: $O(1)$ por operacion.
 */
struct Circle {
  Point<ld> c;
  ld r;
  Circle(Point<ld> c = {}, ld r = 0) : c(c), r(r) {}
};
vector<Point<ld>> circle_line_intersect(Point<ld> a, Point<ld> b, Circle circ) {
  auto dir = b - a;
  ld len2 = norm2(dir);
  if (cmp(len2) == 0) return {};
  ld t = dot(circ.c - a, dir) / len2;
  auto proj = a + dir * t;
  ld d2 = norm2(circ.c - proj);
  if (cmp(d2, circ.r * circ.r) > 0) return {};
  if (cmp(d2, circ.r * circ.r) == 0) return {proj};
  ld off = sqrtl(max((ld)0, circ.r * circ.r - d2) / len2);
  return {proj - dir * off, proj + dir * off};
}
vector<Point<ld>> circle_circle_intersect(Circle c1, Circle c2) {
  auto dir = c2.c - c1.c;
  ld d = norm(dir);
  if (cmp(d, c1.r + c2.r) > 0 || cmp(d, abs(c1.r - c2.r)) < 0 || cmp(d) == 0) return {};
  ld a = (c1.r * c1.r - c2.r * c2.r + d * d) / (2 * d);
  ld h = sqrtl(max((ld)0, c1.r * c1.r - a * a));
  auto p2 = c1.c + dir * (a / d);
  auto perp = Point<ld>(-dir.y, dir.x) * (h / d);
  if (cmp(h) == 0) return {p2};
  return {p2 + perp, p2 - perp};
}

/**
 * Uso: vector<Halfplane> H; H.pb(Halfplane(a, b)); halfplane_intersection(H);
 * Interseccion de semiplanos (valida a la izquierda de cada recta $a \to b$); devuelve el poligono convexo resultante.
 * Complejidad: $O(N \log N)$.
 */
struct Halfplane {
  Point<ld> p, pq; // p = punto, pq = vector direccional
  ld angle;
  Halfplane() {}
  Halfplane(Point<ld> a, Point<ld> b) : p(a), pq(b - a), angle(atan2l(pq.y, pq.x)) {}
  bool out(Point<ld> r) const { return cross(pq, r - p) < -EPS; } // r afuera (lado derecho estricto)
  bool operator<(const Halfplane& e) const {
    if (abs(angle - e.angle) > EPS) return angle < e.angle;
    return cross(pq, e.p - p) > EPS;
  }
};
Point<ld> intersect(Halfplane s, Halfplane t) {
  ld alpha = cross(t.p - s.p, t.pq) / cross(s.pq, t.pq);
  return s.p + s.pq * alpha;
}
vector<Point<ld>> halfplane_intersection(vector<Halfplane>& H) {
  sort(all(H));
  int n = H.size(), k = 0;
  vector<Halfplane> U(n);
  for (int i = 0; i < n; i++) {
    if (i > 0 && abs(H[i].angle - H[i - 1].angle) < EPS) continue;
    U[k++] = H[i];
  }
  H = U; H.resize(k); n = k;
  deque<Halfplane> dq;
  dq.push_back(H[0]);
  if (n > 1) dq.push_back(H[1]);
  for (int i = 2; i < n; i++) {
    while (dq.size() >= 2 && H[i].out(intersect(dq.back(), dq[dq.size() - 2]))) dq.pop_back();
    while (dq.size() >= 2 && H[i].out(intersect(dq.front(), dq[1]))) dq.pop_front();
    dq.push_back(H[i]);
  }
  while (dq.size() >= 3 && dq.front().out(intersect(dq.back(), dq[dq.size() - 2]))) dq.pop_back();
  while (dq.size() >= 3 && dq.back().out(intersect(dq.front(), dq[1]))) dq.pop_front();
  vector<Point<ld>> poly;
  if (dq.size() < 3) return poly;
  for (int i = 0; i < dq.size(); i++) poly.pb(intersect(dq[i], dq[(i + 1) % dq.size()]));
  vector<Point<ld>> clean; // quita vertices duplicados y colineales
  for (int i = 0; i < poly.size(); i++) {
    auto cur = poly[i], nxt = poly[(i + 1) % poly.size()];
    if (hypot(cur.x - nxt.x, cur.y - nxt.y) < EPS) continue;
    if (!clean.empty()) {
      auto prv = clean.back();
      if (abs(cross(cur - prv, nxt - cur)) < EPS) continue;
    }
    clean.pb(cur);
  }
  return clean;
}

/**
 * Uso: Point3D p(x, y, z); dist_point_plane(p, a, n); dist_point_line(p, a, b); line_plane_intersect(L1, L2, P, N);
 * Punto 3D con producto punto y cruz, distancias punto-plano y punto-linea e interseccion linea-plano.
 * Complejidad: $O(1)$ por operacion.
 */
struct Point3D {
  double x, y, z;
  Point3D(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
  Point3D operator+(const Point3D& p) const { return {x + p.x, y + p.y, z + p.z}; }
  Point3D operator-(const Point3D& p) const { return {x - p.x, y - p.y, z - p.z}; }
  Point3D operator*(double c) const { return {x * c, y * c, z * c}; }
  double dot(const Point3D& p) const { return x * p.x + y * p.y + z * p.z; }
  Point3D cross(const Point3D& p) const { return {y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x}; }
  double length() const { return sqrt(dot(*this)); }
};
double dist_point_plane(Point3D p, Point3D a, Point3D n) { return abs((p - a).dot(n)) / n.length(); }
double dist_point_line(Point3D p, Point3D a, Point3D b) { return (p - a).cross(b - a).length() / (b - a).length(); }
Point3D line_plane_intersect(Point3D L1, Point3D L2, Point3D P, Point3D N) {
  auto u = L2 - L1;
  double d = N.dot(u);
  if (abs(d) > EPS) return L1 + u * ((P - L1).dot(N) / d);
  return {1e9, 1e9, 1e9}; // paralelos
}