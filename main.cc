#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cassert>

using namespace std;

typedef unsigned char byte;
typedef vector<byte> bytes;
typedef uint32_t block;

/* === Block cipher === */

typedef byte sbox[16];

const sbox S    = {7, 3, 6, 1, 13, 9, 10, 11, 2, 12, 0, 4, 5, 15, 8, 14};
const sbox Sinv = {10, 3, 8, 1, 11, 12, 2, 0, 14, 5, 6, 7, 9, 4, 15, 13};

block apply_subst(const sbox s, block b) {
  block res = 0;
  for (int i = 0; i < 8; i++) { // From right to left
    block shifted = b >> 4*i;
    byte cut = shifted & 0xF;
    block img = s[cut];
    res = res | (img << 4*i);
  }
  return res;
}

block rotr(block x, byte n)
{
  return (x>>n) | (x<<(32-n)); // compiled into a rorl instruction
}

block B32_turn(block K, block b) {
  return rotr(apply_subst(S, b), 2) xor K;
}

block B32_inv_turn(block K, block b) {
  return apply_subst(Sinv, rotr((b xor K), 30));
}

struct B32_Cipher {
  block K0;
  block K1;
  block K2;
};


block B32_Encode(B32_Cipher c, block b) {
  block xored = b xor c.K0;
  block round1 = B32_turn(c.K1, xored);
  return B32_turn(c.K2, round1);
}

block B32_Decode(B32_Cipher c, block b) {
  block round2 = B32_inv_turn(c.K2, b);
  block round1 = B32_inv_turn(c.K1, round2);
  return round1 xor c.K0;
}

/* === Matrix of linear approximations === */
int compute_L_of(const sbox s, byte a, byte b) {
  int res = 0;
  for (int x = 0; x < 0xF; x++) {
    if ((a ^ x) == (b ^ s[x]))
      res += 1;
  }
  return res;
}

typedef vector<vector<double>> matrix;

matrix populate_linapp_matrix(const sbox s) {
  matrix m;
  for (int i = 0; i < 0xF; i++) {
    m.push_back(vector<double>());
    for (int j = 0; j < 0xF; j++) {
      double mij = (double)compute_L_of(s, i, j) / 16.0;
      m[i].push_back(mij);
    }
  }
  return m;
}

void print_matrix(matrix m) {
  for (auto x : m) {
    for (auto y : x) {
      cout << fixed << setprecision(4);
      cout << y << " ";
    }
    cout << endl;
  }
}

/* ==== main ==== */

int main() {
  block k = 0u | 1u | (1u << 31);
  B32_Cipher test_key = {
    k,
    ~0u,
    ~k
  };

  cout << bitset<32>(test_key.K0) << endl;
  cout << bitset<32>(test_key.K1) << endl;
  cout << bitset<32>(test_key.K2) << endl;
  cout << bitset<32>(B32_Encode(test_key, 0)) << endl;
  cout << bitset<32>(B32_Decode(test_key, B32_Encode(test_key, 0))) << endl;

  print_matrix(populate_linapp_matrix(S));

  return 0;
}
