#include "test_base.hpp"
#include "string.hpp"
using namespace std;
int failures = 0;
#define CHECK(cond) do { if (!(cond)) { cerr << "  FAIL L" << __LINE__ << ": " #cond "\n"; failures++; } } while (0)

int main() {
  // utilidades
  string s = "  Hola Mundo  ";
  trim(s);
  CHECK(s == "Hola Mundo");
  to_lower(s);
  CHECK(s == "hola mundo");
  to_upper(s);
  CHECK(s == "HOLA MUNDO");
  string r = "a,b,c";
  replace_all(r, ",", "-");
  CHECK(r == "a-b-c");
  CHECK(split("a,b,c", ',') == (vector<string>{"a", "b", "c"}));
  CHECK(join(split("a,b,c", ','), "+") == "a+b+c");

  // KMP / Z
  CHECK(prefix_function("aabaaab") == (vi{0, 1, 0, 1, 2, 2, 3}));
  CHECK(prefix_function("aaa") == (vi{0, 1, 2}));
  auto z = z_function("aaaaa");
  CHECK(z[1] == 4 && z[4] == 1);

  // eliminar ocurrencias
  CHECK(remove_all_occurrences("abcabc", "abc") == "");
  CHECK(remove_all_occurrences("hello", "zz") == "hello");

  // rabin-karp
  CHECK(rabin_karp("abracadabra", "abra") == (vi{0, 7}));
  CHECK(rabin_karp("aaaa", "aa") == (vi{0, 1, 2}));

  // manacher: max radio = 4 para "abba"
  auto p = manacher("abba");
  CHECK(*max_element(all(p)) == 4);
  auto p2 = manacher("abcba");
  CHECK(*max_element(all(p2)) == 5);

  // StringHash
  StringHash sh("hello world");
  CHECK(sh.get_hash(0, 4) == sh.get_hash(0, 4));
  CHECK(sh.get_hash(0, 2) != sh.get_hash(1, 3));

  // Trie
  Trie tr;
  tr.insert("hello");
  tr.insert("help");
  CHECK(tr.search("hello") && tr.search("help"));
  CHECK(!tr.search("hel") && !tr.search("hell"));

  // Aho-Corasick
  AhoCorasick ac;
  vector<string> pats = {"a", "ab", "bab", "bc", "bca", "c", "caa"};
  for (auto& pstr : pats) ac.insert(pstr);
  ac.build();
  auto occ = ac.count_occurrences("abccab");
  CHECK(occ == (vi{2, 2, 0, 1, 0, 2, 0}));
  auto first = ac.first_occurrences("abccab");
  CHECK(first == (vi{0, 0, -1, 1, -1, 2, -1}));

  // Suffix automaton
  SuffixAutomaton sam("abc");
  CHECK(sam.distinct_substrings() == 6);
  SuffixAutomaton sam2("abab");
  CHECK(sam2.count_occurrences("ab") == 2);
  CHECK(sam.first_occurrence("bc") == 1);
  CHECK(sam.first_occurrence("zz") == -1);

  // Suffix array
  SuffixArray sa("ababba");
  CHECK(sa.p[0] == 5); // sufijo lexicograficamente menor: "a$"
  CHECK(sa.get_lcp(0, 2) == 2); // "ababba$" vs "abba$"
  CHECK(sa.get_lcp(1, 3) == 1); // "babba$" vs "bba$"

  // Palindromic tree
  {
    PalindromicTree pt("abacaba");
    CHECK(pt.count_palindromes() == 7); // a,b,c,aba,aca,bacab,abacaba
    pt.compute_occurrences();
    // encontrar nodos por longitud y primer posicion para verificar ocurrencias
    // "a" aparece 4 veces (pos 0,2,4,6)
    bool found_a = false;
    for (auto& nd : pt.t) if (nd.len == 1 && nd.occ == 4) found_a = true;
    CHECK(found_a);
    bool found_aba = false;
    for (auto& nd : pt.t) if (nd.len == 3 && nd.occ == 2) found_aba = true;
    CHECK(found_aba);
    bool found_full = false;
    for (auto& nd : pt.t) if (nd.len == 7 && nd.occ == 1 && nd.first_pos == 0) found_full = true;
    CHECK(found_full);
  }
  {
    PalindromicTree pt("aaaa");
    CHECK(pt.count_palindromes() == 4); // a, aa, aaa, aaaa
    pt.compute_occurrences();
    bool found_aa = false;
    for (auto& nd : pt.t) if (nd.len == 2 && nd.occ == 3) found_aa = true;
    CHECK(found_aa);
    bool found_aaaa = false;
    for (auto& nd : pt.t) if (nd.len == 4 && nd.occ == 1 && nd.first_pos == 0) found_aaaa = true;
    CHECK(found_aaaa);
  }
  {
    PalindromicTree pt("abccba"); // palindromos: a, b, c, cc, bccb, abccba
    CHECK(pt.count_palindromes() == 6);
    pt.compute_occurrences();
    bool found_long = false;
    for (auto& nd : pt.t) if (nd.len == 6 && nd.occ == 1 && nd.first_pos == 0) found_long = true;
    CHECK(found_long);
    bool found_cc = false;
    for (auto& nd : pt.t) if (nd.len == 2 && nd.occ == 1 && nd.first_pos == 2) found_cc = true;
    CHECK(found_cc);
  }

  if (failures) { cout << "test_string: " << failures << " fallos\n"; return 1; }
  cout << "test_string: OK\n";
  return 0;
}