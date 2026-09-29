#include "test_base.hpp"
#include "bits.hpp"
using namespace std;
int failures = 0;
#define CHECK(cond) do { if (!(cond)) { cerr << "  FAIL L" << __LINE__ << ": " #cond "\n"; failures++; } } while (0)

int main() {
  // utilidades O(1)
  CHECK(check_bit(5, 0) && !check_bit(5, 1) && check_bit(5, 2));
  CHECK(set_bit(0, 3) == 8);
  CHECK(clear_bit(7, 1) == 5);
  CHECK(toggle_bit(5, 1) == 7);
  CHECK(lsb(12) == 4);
  CHECK(count_bits_ll(255) == 8);
  CHECK(highest_bit(18) == 4 && highest_bit(0) == -1);
  CHECK(ctz(12) == 2);

  // gospers hack: mascaras de 2 bits en 4
  auto masks = gospers_hack(2, 4);
  CHECK(masks == (vll{3, 5, 6, 9, 10, 12}));

  // XorBasis
  XorBasis b;
  b.insert(1); b.insert(2); b.insert(3); // 3 = 1^2, no independiente
  CHECK(b.sz == 2);
  CHECK(b.can_form(3));
  CHECK(!b.can_form(5));
  CHECK(b.get_min() == 1);
  b.insert(4); b.insert(8);
  CHECK(b.get_max() == 15);
  CHECK(b.get_max(5) == 15); // 5 ^ 2 ^ 8 = 15

  // BitTrie: max XOR contra el arreglo
  BitTrie trie;
  trie.insert(3); trie.insert(10); trie.insert(5);
  CHECK(trie.max_xor(6) == 12); // 6 ^ 10
  CHECK(trie.max_xor(0) == 10);
  trie.insert(0);
  CHECK(trie.max_xor(5) == 15); // 5 ^ 10
  trie.insert(5, -1); // borrar
  CHECK(trie.max_xor(5) == 15); // sigue 5 ^ 10
  trie.insert(10, -1);
  CHECK(trie.max_xor(5) == 6); // 5 ^ 3

  if (failures) { cout << "test_bits: " << failures << " fallos\n"; return 1; }
  cout << "test_bits: OK\n";
  return 0;
}