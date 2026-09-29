#include "test_base.hpp"
#include "dp.hpp"
using namespace std;
int failures = 0;
#define CHECK(cond) do { if (!(cond)) { cerr << "  FAIL L" << __LINE__ << ": " #cond "\n"; failures++; } } while (0)

int main() {
  // Li Chao (minimizar)
  LiChaoTree cht;
  cht.add(-2, 10);
  cht.add(1, 2);
  CHECK(cht.query(3) == 4);  // -2*3+10 = 4 vs 1*3+2 = 5
  CHECK(cht.query(10) == -10); // -20+10 = -10 vs 12
  cht.add(0, 0);
  CHECK(cht.query(0) == 0);

  // LIS + reconstruccion
  auto [len, seq] = lis({10, 22, 9, 33, 21, 50, 41, 60});
  CHECK(len == 5);
  CHECK(seq.size() == 5);
  CHECK(is_sorted(all(seq)));
  CHECK(lis({}).fi == 0);
  CHECK(lis({7}).fi == 1);

  // bitset knapsack
  CHECK(bitset_knapsack({2, 5, 8}, 13));
  CHECK(bitset_knapsack({2, 5, 8}, 7)); // 2 + 5
  CHECK(!bitset_knapsack({2, 5, 8}, 4));
  CHECK(bitset_knapsack({5}, 5));
  CHECK(bitset_knapsack({3}, 0));

  // SOS DP: dp[3] = suma de todos los subsets de {0,1,2,3}
  vll dp = {1, 2, 3, 4};
  sos_dp(dp, 2);
  CHECK(dp[3] == 10 && dp[1] == 3 && dp[0] == 1);

  // Knuth: cost(i,j)=1, n=4 -> dp[0][3]=2
  auto cost = [](int i, int j) { return 1; };
  auto kn = knuth_dp(4, cost);
  CHECK(kn[0][3] == 2 && kn[0][2] == 1);

  // DnC DP: particionar, cost 0, prev[0]=0
  vll prev(4, LINF), curr(4, 0);
  prev[0] = 0;
  auto cost0 = [](int l, int r) { return 0; };
  dnc_dp(1, 3, 0, 3, prev, curr, cost0);
  CHECK(curr[1] == 0 && curr[2] == 0 && curr[3] == 0);

  // Digit DP: todos los numeros son validos -> n+1
  CHECK(DigitDP::solve(0) == 1);
  CHECK(DigitDP::solve(9) == 10);
  CHECK(DigitDP::solve(99) == 100);

  // CHT con deque (min, pendientes y x crecientes)
  {
    CHTDeque cht;
    cht.add(1, 0);
    cht.add(2, -1);
    cht.add(3, -5);
    // brute force
    vector<array<ll, 2>> lines = {{1, 0}, {2, -1}, {3, -5}};
    for (ll x : {0LL, 1LL, 2LL, 3LL, 10LL}) {
      ll best = LINF;
      for (auto& [m, c] : lines) best = min(best, m * x + c);
      CHECK(cht.query(x) == best);
      CHECK(cht.query_bs(x) == best);
    }
    // max: insertar (-m, -c) y negar
    CHTDeque cht_max;
    for (auto& [m, c] : lines) cht_max.add(-m, -c);
    for (ll x : {0LL, 1LL, 2LL, 5LL}) {
      ll best = -LINF;
      for (auto& [m, c] : lines) best = max(best, m * x + c);
      CHECK(-cht_max.query(x) == best);
      CHECK(-cht_max.query_bs(x) == best);
    }
  }

  // WQS: particionar en exactamente K segmentos minimizando suma de (suma)^2
  {
    vll arr = {1, 2, 3, 4, 1, 2};
    int n = arr.size();
    vll pref(n + 1, 0);
    for (int i = 0; i < n; i++) pref[i + 1] = pref[i] + arr[i];
    auto check = [&](ll lam) -> pair<ll, int> {
      vll dp(n + 1, LINF);
      vi cnt(n + 1, 0);
      dp[0] = 0;
      for (int i = 1; i <= n; i++) {
        for (int j = 0; j < i; j++) {
          ll seg = pref[i] - pref[j];
          ll cand = dp[j] + seg * seg + lam;
          if (cand < dp[i] || (cand == dp[i] && cnt[j] + 1 > cnt[i])) {
            dp[i] = cand;
            cnt[i] = cnt[j] + 1;
          }
        }
      }
      return {dp[n], cnt[n]};
    };
    // brute force: particionar en exactamente k segmentos
    auto brute = [&](int k) {
      ll best = LINF;
      function<void(int, int, ll)> rec = [&](int pos, int segs, ll cost) {
        if (segs == 0) {
          if (pos == n) best = min(best, cost);
          return;
        }
        for (int end = pos + 1; end <= n - (segs - 1); end++) {
          ll s = pref[end] - pref[pos];
          rec(end, segs - 1, cost + s * s);
        }
      };
      rec(0, k, 0);
      return best;
    };
    for (int k = 1; k <= n; k++) {
      CHECK(WQS::solve(check, k) == brute(k));
    }
  }

  if (failures) { cout << "test_dp: " << failures << " fallos\n"; return 1; }
  cout << "test_dp: OK\n";
  return 0;
}