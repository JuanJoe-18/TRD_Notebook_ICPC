#include "test_base.hpp"
#include "extra.hpp"
using namespace std;
int failures = 0;
#define CHECK(cond) do { if (!(cond)) { cerr << "  FAIL L" << __LINE__ << ": " #cond "\n"; failures++; } } while (0)

int main() {
  // Meet in the middle
  CHECK(MeetInTheMiddle::subset_sum({10, 20, 15, 5, 30}, 35));
  CHECK(MeetInTheMiddle::subset_sum({10, 20, 15, 5, 30}, 80));
  CHECK(!MeetInTheMiddle::subset_sum({10, 20, 15, 5, 30}, 12));
  CHECK(MeetInTheMiddle::subset_sum({7}, 7));
  CHECK(MeetInTheMiddle::subset_sum({}, 0));

  // Zobrist hashing: anagramas
  ZobristHash zh;
  vi A = {1, 2, 3, 2, 1}, B = {3, 1, 1, 2, 2};
  auto pa = zh.build_pref(A), pb = zh.build_pref(B);
  CHECK(zh.query_range(pa, 0, 4) == zh.query_range(pb, 0, 4));
  CHECK(zh.query_range(pa, 0, 0) == zh.query_range(pa, 0, 0));
  CHECK(zh.query_range(pa, 1, 3) == zh.query_range(pa, 1, 3));

  // Bitset reachability (cierre transitivo)
  BitsetReachability<2500> br;
  br.build(5, {{1, 2}, {2, 3}, {4, 5}});
  CHECK(br.can_reach(1, 3));
  CHECK(br.can_reach(1, 1));
  CHECK(!br.can_reach(1, 5));
  CHECK(!br.can_reach(3, 1));
  CHECK(br.can_reach(4, 5));

  // Compresion de coordenadas
  {
    auto [vals, ids] = compress_coords(vi{5, 3, 3, 1, 5, 2});
    CHECK(vals == (vi{1, 2, 3, 5}));
    CHECK(ids == (vi{3, 2, 2, 0, 3, 1}));
  }
  {
    auto [vals, ids] = compress_coords(vi{1, 2, 3});
    CHECK(vals == (vi{1, 2, 3}) && ids == (vi{0, 1, 2}));
  }
  {
    auto [vals, ids] = compress_coords(vi{-3, 10, -3, 0});
    CHECK(vals == (vi{-3, 0, 10}) && ids == (vi{0, 2, 0, 1}));
  }
  {
    auto [vals, ids] = compress_coords(vi{7});
    CHECK(vals == (vi{7}) && ids == (vi{0}));
  }
  {
    auto [vals, ids] = compress_coords(vll{});
    CHECK(vals.empty() && ids.empty());
  }
  // robustez aleatoria vs brute force
  {
    mt19937 rnd(42);
    for (int it = 0; it < 20; it++) {
      int n = rnd() % 15;
      vll a(n);
      for (ll& x : a) x = (ll)(rnd() % 20) - 10;
      auto [vals, ids] = compress_coords(a);
      vll exp_vals = a;
      sort(all(exp_vals));
      exp_vals.erase(unique(all(exp_vals)), exp_vals.end());
      CHECK(vals == exp_vals);
      bool ok = true;
      for (int i = 0; i < n; i++) {
        int exp_id = lower_bound(all(exp_vals), a[i]) - exp_vals.begin();
        ok &= (ids[i] == exp_id);
      }
      CHECK(ok);
    }
  }

  if (failures) { cout << "test_extra: " << failures << " fallos\n"; return 1; }
  cout << "test_extra: OK\n";
  return 0;
}