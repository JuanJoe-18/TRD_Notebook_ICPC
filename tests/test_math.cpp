#include "test_base.hpp"
#include "math.hpp"
using namespace std;
int failures = 0;
#define CHECK(cond) do { if (!(cond)) { cerr << "  FAIL L" << __LINE__ << ": " #cond "\n"; failures++; } } while (0)

int main() {
  // modpow / modinv / crt
  CHECK(MathAlgo::modpow(2, 10) == 1024);
  CHECK(MathAlgo::modpow(3, 5, 13) == 9);
  CHECK(MathAlgo::modinv(3, 7) == 5);
  CHECK(MathAlgo::modinv(2, 5) == 3);
  auto [x, lcm] = MathAlgo::crt(2, 3, 3, 5);
  CHECK(x == 8 && lcm == 15);
  auto [x2, lcm2] = MathAlgo::crt(1, 2, 3, 4); // x = 3 mod 4
  CHECK(x2 == 3 && lcm2 == 4);
  CHECK(MathAlgo::crt(1, 2, 2, 4).fi == -1); // 1 mod 2 y 2 mod 4: inconsistente

  // Sieve
  MathAlgo::Sieve sv(100);
  auto f = sv.factorize(100);
  CHECK(f == (vi{2, 2, 5, 5}));
  CHECK(sv.primes.size() == 25);
  CHECK(sv.mu[4] == 0 && sv.mu[6] == 1);
  CHECK(sv.phi[10] == 4 && sv.phi[1] == 1);

  // Combinatorics
  MathAlgo::Combinatorics comb(20);
  CHECK(comb.nCr(10, 3) == 120);
  CHECK(comb.nCr(5, 0) == 1 && comb.nCr(5, 6) == 0);
  CHECK(comb.catalan(5) == 42);

  // Gauss: x + y = 5, x - y = 1 -> x=3, y=2
  vector<vector<double>> A = {{1, 1, 5}, {1, -1, 1}};
  vector<double> ans;
  int sols = MathAlgo::gauss(A, ans);
  CHECK(sols == 1);
  CHECK(fabs(ans[0] - 3) < 1e-6 && fabs(ans[1] - 2) < 1e-6);

  // Matrix
  MathAlgo::Matrix A1(2, 2), A2(2, 2);
  A1.mat = {{1, 2}, {3, 4}};
  A2.mat = {{5, 6}, {7, 8}};
  auto P = A1 * A2;
  CHECK(P.mat[0][0] == 19 && P.mat[0][1] == 22 && P.mat[1][0] == 43 && P.mat[1][1] == 50);
  auto I = MathAlgo::Matrix::identity(3);
  CHECK(I.power(5).mat[0][0] == 1);

  // FFT multiply
  auto poly = MathAlgo::multiply({1, 2, 3}, {4, 5, 6});
  CHECK(poly == (vi{4, 13, 28, 27, 18}));

  // berlekamp-massey sobre Fibonacci
  auto C = MathAlgo::berlekamp_massey({0, 1, 1, 2, 3, 5, 8, 13});
  CHECK(C.size() == 3 && C[0] == 1 && C[1] == MOD - 1 && C[2] == MOD - 1);

  // PollardRho
  CHECK(MathAlgo::PollardRho::is_prime(2) && MathAlgo::PollardRho::is_prime(97));
  CHECK(!MathAlgo::PollardRho::is_prime(100));
  CHECK(MathAlgo::PollardRho::is_prime(1000000007));
  auto fac = MathAlgo::PollardRho::factorize(60);
  CHECK(fac[2] == 2 && fac[3] == 1 && fac[5] == 1);

  // Lucas (P primo ~1e6, n,r normales)
  MathAlgo::Lucas lucas(1000003);
  CHECK(lucas.solve(10, 3) == 120);
  CHECK(lucas.solve(5, 2) == 10);

  // sumatorias
  CHECK(MathAlgo::sum_n(10) == 55);
  CHECK(MathAlgo::sum_squares(10) == 385);
  CHECK(MathAlgo::geom_sum(1, 2, 10) == 1023);

  // nim
  CHECK(!MathAlgo::nim_game({1, 2, 3}));
  CHECK(MathAlgo::nim_game({1, 2}));

  // FWHT: convoluciones XOR / OR / AND vs brute force
  {
    using namespace MathAlgo;
    // casos exactos pequenos
    CHECK(xor_convolution({1, 2}, {3, 4}) == (vll{11, 10}));
    CHECK(or_convolution({1, 2}, {3, 4}) == (vll{3, 18}));
    CHECK(and_convolution({1, 2}, {3, 4}) == (vll{13, 8}));
    // brute force aleatorio (valores mod MOD)
    mt19937 rnd(123);
    for (int it = 0; it < 30; it++) {
      int n = 1 << (1 + rnd() % 4);
      vll a(n), b(n);
      for (int i = 0; i < n; i++) a[i] = rnd() % MOD, b[i] = rnd() % MOD;
      vll ex_xor(n, 0), ex_or(n, 0), ex_and(n, 0);
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
          ex_xor[i ^ j] = (ex_xor[i ^ j] + a[i] * b[j]) % MOD;
          ex_or[i | j] = (ex_or[i | j] + a[i] * b[j]) % MOD;
          ex_and[i & j] = (ex_and[i & j] + a[i] * b[j]) % MOD;
        }
      CHECK(xor_convolution(a, b) == ex_xor);
      CHECK(or_convolution(a, b) == ex_or);
      CHECK(and_convolution(a, b) == ex_and);
    }
  }

  // floor_sum vs brute force
  {
    CHECK(MathAlgo::floor_sum(4, 5, 3, 1) == 3);
    CHECK(MathAlgo::floor_sum(5, 7, 2, 3) == 3);
    mt19937 rnd(456);
    for (int it = 0; it < 30; it++) {
      ll n = 1 + rnd() % 50, m = 1 + rnd() % 100, a = rnd() % 100, b = rnd() % 100;
      ll brute = 0;
      for (ll i = 0; i < n; i++) brute += (a * i + b) / m;
      CHECK(MathAlgo::floor_sum(n, m, a, b) == brute);
    }
  }

  if (failures) { cout << "test_math: " << failures << " fallos\n"; return 1; }
  cout << "test_math: OK\n";
  return 0;
}