#include "test_base.hpp"
#include "geometry.hpp"
using namespace std;
int failures = 0;
#define CHECK(cond) do { if (!(cond)) { cerr << "  FAIL L" << __LINE__ << ": " #cond "\n"; failures++; } } while (0)
#define NEAR(a, b) CHECK(fabs((a) - (b)) < 1e-6)

int main() {
  // operaciones de Point
  Point<ll> a(1, 2), b(3, 4);
  CHECK((a + b) == Point<ll>(4, 6));
  CHECK(dot(a, b) == 11);
  CHECK(cross(a, b) == -2);
  CHECK(orient(Point<ll>(0, 0), Point<ll>(1, 0), Point<ll>(0, 1)) == 1); // ccw
  CHECK(norm2(a) == 5);

  // interseccion de segmentos
  CHECK(segment_intersect(Point<ll>(0, 0), Point<ll>(4, 0), Point<ll>(2, -1), Point<ll>(2, 1)));
  CHECK(!segment_intersect(Point<ll>(0, 0), Point<ll>(2, 0), Point<ll>(0, 1), Point<ll>(2, 1)));
  CHECK(segment_intersect(Point<ll>(0, 0), Point<ll>(2, 0), Point<ll>(1, 0), Point<ll>(3, 0))); // colineal solapado
  CHECK(on_segment(Point<ll>(1, 0), Point<ll>(0, 0), Point<ll>(2, 0)));

  // area (doble) y borde del rectangulo 4x3
  vector<Point<ll>> rect = {{0, 0}, {4, 0}, {4, 3}, {0, 3}};
  CHECK(polygon_area_2(rect) == 24);
  CHECK(boundary_points(rect) == 14);
  CHECK(point_in_polygon(rect, Point<ll>(2, 1)) == 2);  // adentro
  CHECK(point_in_polygon(rect, Point<ll>(4, 1)) == 1);  // borde
  CHECK(point_in_polygon(rect, Point<ll>(5, 5)) == 0);  // afuera
  CHECK(point_in_polygon(rect, Point<ll>(0, 0)) == 1);  // vertice

  // convex hull
  auto hull = convex_hull(vector<Point<ll>>{{0, 0}, {2, 0}, {2, 2}, {0, 2}, {1, 1}});
  CHECK(hull.size() == 4);
  auto hull2 = convex_hull(vector<Point<ll>>{{0, 0}, {1, 1}}); // 2 puntos
  CHECK(hull2.size() == 2);

  // closest pair (distancia al cuadrado)
  CHECK(closest_pair({Point<ll>(0, 0), Point<ll>(100, 100), Point<ll>(1000, 1000), Point<ll>(102, 102)}) == 8);
  CHECK(closest_pair({Point<ll>(1, 1), Point<ll>(4, 5)}) == 25);

  // rotating calipers sobre el rectangulo
  CHECK(rotating_calipers(rect) == 25);

  // interseccion de lineas infinitas (y = x vs y = -x + 10)
  auto inter = line_intersect(Point<ld>(0, 0), Point<ld>(10, 10), Point<ld>(0, 10), Point<ld>(10, 0));
  NEAR(inter.x, 5); NEAR(inter.y, 5);

  // distancias a linea/segmento
  NEAR(dist_to_segment(Point<ll>(0, 0), Point<ll>(3, 0), Point<ll>(3, 3)), 3);
  NEAR(dist_to_segment(Point<ll>(4, 0), Point<ll>(3, 0), Point<ll>(3, 3)), 1);

  // circulo vs linea horizontal y=3 (r=5) -> x = +-4
  auto cl = circle_line_intersect(Point<ld>(0, 3), Point<ld>(10, 3), Circle(Point<ld>(0, 0), 5));
  CHECK(cl.size() == 2);
  NEAR(fabs(cl[0].x), 4); NEAR(fabs(cl[1].x), 4);
  auto cl2 = circle_line_intersect(Point<ld>(0, 7), Point<ld>(10, 7), Circle(Point<ld>(0, 0), 5));
  CHECK(cl2.empty());

  // circulo vs circulo (0,0,5) y (8,0,5) -> (4, +-3)
  auto cc = circle_circle_intersect(Circle(Point<ld>(0, 0), 5), Circle(Point<ld>(8, 0), 5));
  CHECK(cc.size() == 2);
  NEAR(cc[0].x, 4); NEAR(fabs(cc[0].y), 3);

  // ordenamiento polar antihorario
  vector<Point<ll>> pts = {{0, -1}, {-1, 0}, {0, 1}, {1, 0}};
  polar_sort(pts);
  CHECK(pts[0] == Point<ll>(1, 0) && pts[1] == Point<ll>(0, 1) && pts[2] == Point<ll>(-1, 0) && pts[3] == Point<ll>(0, -1));

  // interseccion de semiplanos: cuadrado 10x10
  vector<Halfplane> H = {
    Halfplane(Point<ld>(0, 0), Point<ld>(10, 0)),
    Halfplane(Point<ld>(10, 0), Point<ld>(10, 10)),
    Halfplane(Point<ld>(10, 10), Point<ld>(0, 10)),
    Halfplane(Point<ld>(0, 10), Point<ld>(0, 0)),
  };
  auto sq = halfplane_intersection(H);
  CHECK(sq.size() == 4);

  // 3D
  NEAR(dist_point_plane(Point3D(3, 4, 5), Point3D(0, 0, 0), Point3D(0, 1, 0)), 4);
  NEAR(dist_point_line(Point3D(1, 2, 3), Point3D(0, 0, 0), Point3D(1, 0, 0)), sqrt(13.0));

  if (failures) { cout << "test_geometry: " << failures << " fallos\n"; return 1; }
  cout << "test_geometry: OK\n";
  return 0;
}