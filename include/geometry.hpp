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
 * Uso: dot(a, b), cross(a, b), cross(p, a, b), norm2(p);
 * Operaciones que devuelven el mismo tipo $T$ del punto: producto punto, producto cruz (escalar) y norma al cuadrado.
 * Complejidad: $O(1)$.
 */
template<typename T> T dot(Point<T> a, Point<T> b) { return a.x * b.x + a.y * b.y; }
template<typename T> T cross(Point<T> a, Point<T> b) { return a.x * b.y - a.y * b.x; }
template<typename T> T cross(Point<T> p, Point<T> a, Point<T> b) { return cross(a - p, b - p); }
template<typename T> T norm2(Point<T> p) { return dot(p, p); }
/**
 * Uso: norm(p);
 * Norma (longitud) de un punto o vector; devuelve long double.
 * Complejidad: $O(1)$.
 */
template<typename T> ld norm(Point<T> p) { return sqrtl(norm2(p)); }
/**
 * Uso: sgn(v);
 * Signo de un valor: devuelve 1 si $v > 0$, -1 si $v < 0$ y 0 si $v = 0$.
 * Complejidad: $O(1)$.
 */
template<typename T> int sgn(T v) { return (T(0) < v) - (v < T(0)); }
/**
 * Uso: orient(a, b, c);
 * Orientacion del giro $a \to b \to c$: devuelve +1 si es antihorario (izquierda), -1 si es horario (derecha) y 0 si son colineales.
 * Complejidad: $O(1)$.
 */
template<typename T> int orient(Point<T> a, Point<T> b, Point<T> c) { return sgn(cross(a, b, c)); }

/**
 * Uso: on_segment(p, a, b);
 * Devuelve true si el punto $p$ pertenece al segmento $[a, b]$ (incluyendo los extremos).
 * Complejidad: $O(1)$.
 */
template<typename T> bool on_segment(Point<T> p, Point<T> a, Point<T> b) {
  return orient(a, b, p) == 0 && min(a.x, b.x) <= p.x && p.x <= max(a.x, b.x) &&
         min(a.y, b.y) <= p.y && p.y <= max(a.y, b.y);
}
/**
 * Uso: segment_intersect(a, b, c, d);
 * Devuelve true si los segmentos $[a, b]$ y $[c, d]$ se intersectan (incluye tocarse en un punto o superponerse).
 * Complejidad: $O(1)$.
 */
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
 * Uso: polygon_area_2(poly);
 * Devuelve el doble del area del poligono (para evitar floats); el area real es $area / 2$. Vertices en orden horario o antihorario.
 * Complejidad: $O(N)$.
 */
template<typename T> T polygon_area_2(const vector<Point<T>>& p) {
  T area = 0;
  for (int i = 0; i < p.size(); i++) area += cross(p[i], p[(i + 1) % p.size()]);
  return abs(area);
}
/**
 * Uso: boundary_points(poly);
 * Devuelve (long long) la cantidad de puntos del reticulo sobre el borde de un poligono (teorema de Pick: $2A = 2I + B - 2$).
 * Complejidad: $O(N)$.
 */
long long boundary_points(const vector<Point<ll>>& poly) {
  ll b = 0;
  for (int i = 0; i < poly.size(); i++) {
    auto p1 = poly[i], p2 = poly[(i + 1) % poly.size()];
    b += gcd(abs(p1.x - p2.x), abs(p1.y - p2.y));
  }
  return b;
}
/**
 * Uso: point_in_polygon(poly, p);
 * Clasifica un punto respecto a un poligono: devuelve 0 fuera, 1 en el borde, 2 dentro.
 * Complejidad: $O(N)$.
 */
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
 * Cascara convexa (monotonic chain) en sentido antihorario, preservando puntos colineales en el borde; devuelve vector<Point<T>>.
 * Complejidad: $O(N \log N)$.
 */
template<typename T> vector<Point<T>> convex_hull(vector<Point<T>> pts) {
  int n = pts.size(), k = 0;
  if (n <= 2) return pts;
  vector<Point<T>> h(2 * n);
  sort(all(pts));
  for (int i = 0; i < n; i++) {
    while (k >= 2 && cross(h[k - 2], h[k - 1], pts[i]) < 0) k--;
    h[k++] = pts[i];
  }
  for (int i = n - 2, t = k + 1; i >= 0; i--) {
    while (k >= t && cross(h[k - 2], h[k - 1], pts[i]) < 0) k--;
    h[k++] = pts[i];
  }
  h.resize(k - 1);
  return h;
}

/**
 * Uso: closest_pair(puntos);
 * Devuelve (long long) la distancia al cuadrado entre el par de puntos mas cercanos.
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
 * Uso: Event1D ev = {x, type, id};
 * Evento de sweep line 1D: type +1 marca el inicio (entrada) de un rango en $x$ y type -1 su fin (salida); id identifica al intervalo.
 * Al ordenar por $x$ se da prioridad a las entradas: en empate de $x$, primero el de mayor type (+1 antes que -1).
 * Complejidad: $O(1)$ por evento, $O(K \log K)$ al ordenarlos con $K$ eventos.
 */
struct Event1D {
  ll x;
  int type, id;
  bool operator<(const Event1D& o) const { return x != o.x ? x < o.x : type > o.type; }
};
/**
 * Uso: Event2D ev = {x, y1, y2, type};
 * Evento de sweep line 2D para rectangulos: cubre el rango vertical $[y1, y2]$ (ambos extremos inclusivos); type +1 borde izquierdo (entrada) y -1 borde derecho (salida).
 * Al ordenar por $x$ se da prioridad a las entradas: en empate de $x$, primero el de mayor type (+1 antes que -1).
 * Complejidad: $O(1)$ por evento, $O(K \log K)$ al ordenarlos con $K$ eventos.
 */
struct Event2D {
  ll x, y1, y2;
  int type; // +1 borde izq, -1 borde der
  bool operator<(const Event2D& o) const { return x != o.x ? x < o.x : type > o.type; }
};

/**
 * Uso: polar_sort(puntos, centro);
 * Ordena los puntos en su lugar por angulo polar alrededor de un centro, en sentido antihorario desde el eje $+x$.
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
 * Uso: dist_to_line(p, a, b);
 * Devuelve la distancia (long double) del punto $p$ a la linea infinita que pasa por $a$ y $b$.
 * Complejidad: $O(1)$.
 */
template<typename T> ld dist_to_line(Point<T> p, Point<T> a, Point<T> b) {
  return abs((ld)cross(p, a, b)) / norm(a - b);
}
/**
 * Uso: dist_to_segment(p, a, b);
 * Devuelve la distancia (long double) del punto $p$ al segmento $[a, b]$.
 * Complejidad: $O(1)$.
 */
template<typename T> ld dist_to_segment(Point<T> p, Point<T> a, Point<T> b) {
  if (a == b) return norm(p - a);
  if (dot(p - a, b - a) <= 0) return norm(p - a);
  if (dot(p - b, a - b) <= 0) return norm(p - b);
  return dist_to_line(p, a, b);
}
/**
 * Uso: line_intersect(a, b, c, d);
 * Devuelve el punto de interseccion (Point<long double>) de las lineas infinitas $AB$ y $CD$; no es valido si son paralelas.
 * Complejidad: $O(1)$.
 */
Point<ld> line_intersect(Point<ld> a, Point<ld> b, Point<ld> c, Point<ld> d) {
  ld cr1 = cross(b - a, c - a), cr2 = cross(b - a, d - a);
  return c + (d - c) * (cr1 / (cr1 - cr2));
}

/**
 * Uso: rotating_calipers(hull);
 * Devuelve el diametro al cuadrado de un poligono convexo mediante calipers rotativos (mismo tipo $T$ del punto).
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
 * Uso: Circle c(centro, radio);
 * Circulo con centro Point<long double> y radio long double.
 * Complejidad: $O(1)$.
 */
struct Circle {
  Point<ld> c;
  ld r;
  Circle(Point<ld> c = {}, ld r = 0) : c(c), r(r) {}
};
/**
 * Uso: circle_line_intersect(a, b, circ);
 * Devuelve vector<Point<long double>> con los 0, 1 o 2 puntos donde la linea $AB$ corta al circulo.
 * Complejidad: $O(1)$.
 */
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
/**
 * Uso: circle_circle_intersect(c1, c2);
 * Devuelve vector<Point<long double>> con los 0, 1 o 2 puntos donde se cortan dos circulos (vacio si son concentricos o no se tocan).
 * Complejidad: $O(1)$.
 */
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
 * Uso: Halfplane hp(a, b);
 * Un semiplano se crea con dos puntos: el plano valido es el que queda a la izquierda de la recta orientada $a \to b$ (producto cruz positivo).
 * Complejidad: $O(1)$ por semiplano.
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
/**
 * Uso: auto p = intersect(hp1, hp2);
 * Devuelve el punto de corte (Point<long double>) entre las rectas de dos semiplanos.
 * Complejidad: $O(1)$.
 */
Point<ld> intersect(Halfplane s, Halfplane t) {
  ld alpha = cross(t.p - s.p, t.pq) / cross(s.pq, t.pq);
  return s.p + s.pq * alpha;
}
/**
 * Uso: auto poly = halfplane_intersection(H);
 * Interseccion de un conjunto de semiplanos (cada uno valido a la izquierda de su recta); devuelve el poligono convexo resultante (vacio si la region es ilimitada o no existe).
 * Complejidad: $O(N \log N)$.
 */
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
 * Uso: Point3D p(x, y, z); p.dot(q); p.cross(q); p.length();
 * Punto 3D con operadores aritmeticos, producto punto, producto cruz y norma.
 * Complejidad: $O(1)$.
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
/**
 * Uso: dist_point_plane(p, a, n), dist_point_line(p, a, b);
 * Distancia (double) de un punto a un plano (dado por el punto $a$ y la normal $n$) y a una linea (por los puntos $a$ y $b$).
 * Complejidad: $O(1)$.
 */
double dist_point_plane(Point3D p, Point3D a, Point3D n) { return abs((p - a).dot(n)) / n.length(); }
double dist_point_line(Point3D p, Point3D a, Point3D b) { return (p - a).cross(b - a).length() / (b - a).length(); }
/**
 * Uso: line_plane_intersect(L1, L2, P, N);
 * Punto de interseccion (Point3D) de la linea $L1 \to L2$ con el plano definido por el punto $P$ y la normal $N$; devuelve $(10^9, 10^9, 10^9)$ si son paralelos.
 * Complejidad: $O(1)$.
 */
Point3D line_plane_intersect(Point3D L1, Point3D L2, Point3D P, Point3D N) {
  auto u = L2 - L1;
  double d = N.dot(u);
  if (abs(d) > EPS) return L1 + u * ((P - L1).dot(N) / d);
  return {1e9, 1e9, 1e9}; // paralelos
}