// #pragma GCC optimize("O3,unroll-loops")
// #pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")

#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace std;
using namespace __gnu_pbds;

#define fastio ios::sync_with_stdio(false); cin.tie(0); cout.tie(0)
#define all(v) v.begin(), v.end()
#define rall(v) v.rbegin(), v.rend()
#define pb push_back
#define fi first
#define se second

typedef long long ll;
typedef long double ld;
typedef pair<int,int> pii;
typedef vector<int> vi;
typedef vector<ll> vll;

const int INF = 1e9;
const ll LINF = 1e18;
const int MOD = 1e9+7;
const double EPS = 1e-9;

ll rand_val(ll a, ll b) {
  static mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
  return uniform_int_distribution<ll>(a, b)(rng);
}

void solve(){}

int main() {
  #ifdef ICESI
    freopen("input.txt", "r", stdin);
  #endif
  fastio;
  int t = 1;
  // cin >> t;
  while (t--) {
    solve();
  }
  return 0;
}
