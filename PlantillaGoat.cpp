/**
 * @file PlantillaGoat.cpp
 * @brief Plantilla general para problemas de programacion competitiva.
 * @details Utilidades de E/S, matematicas, grafos, rangos y strings.
 * @note Elimina las secciones que no necesites antes de enviar la solucion.
 */
//   ____ ___  ____  _____   ____  _   _
//  / ___/ _ \|  _ \| ____| / ___|| | | |
// | |  | | | | | | |  _|   \___ \| | | |
// | |__| |_| | |_| | |___   ___) | |_| |
//  \____\___/|____/|_____| |____/ \___/
//
//              COMPETITIVE PROGRAMMING TEMPLATE

// #pragma GCC optimize("O3,unroll-loops")
// #pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")

#include <bits/stdc++.h>
using namespace std;

// ================================
// Entrada/Salida rapida
// ================================
#define fastio ios::sync_with_stdio(false); cin.tie(0); cout.tie(0)

// ================================
// Atajos
// ================================
#define all(v) v.begin(), v.end()
#define rall(v) v.rbegin(), v.rend()
#define pb push_back
#define fi first
#define se second

typedef long long ll;
typedef pair<int,int> pii;
typedef vector<int> vi;
typedef vector<ll> vll;

// ================================
// Constantes globales
// ================================
const int INF = 1e9;
const ll LINF = 1e18;
const int MOD = 1e9+7;  // cambiar según el problema
const double EPS = 1e-9;


// ================================
// MAIN
// ================================
int main(){
    fastio;
    int t=1;
    // cin >> t; // descomentar si hay multiples casos
    while(t--){
        // ---------------------------
        // Aquí resuelves el problema
        // ---------------------------

        int n; cin >> n;
        vector<int> a(n);
        for(int i=0;i<n;i++) cin >> a[i];

        // ejemplo: suma
        ll sum = accumulate(all(a),0LL);
        cout << sum << "\n";
    }
    return 0;
}
