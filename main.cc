#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cassert>
#include <random>
#include <time.h>

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
  return rotr(apply_subst(S, b), 2) ^ K;
}

block B32_inv_turn(block K, block b) {
  return apply_subst(Sinv, rotr((b ^ K), 30));
}

struct B32_Cipher {
  block K0;
  block K1;
  block K2;
};


block B32_Encode(B32_Cipher c, block b) {
  block xored = b ^ c.K0;
  block round1 = B32_turn(c.K1, xored);
  return B32_turn(c.K2, round1);
}

block B32_Decode(B32_Cipher c, block b) {
  block round2 = B32_inv_turn(c.K2, b);
  block round1 = B32_inv_turn(c.K1, round2);
  return round1 ^ c.K0;
}

/* === Matrix of linear approximations === */
int strange_op(block a, block b) {
  int res = 0;
  for (int i = 0; i < 32; i++) {
    byte a1 = (a & (1 << i)) >> i;
    byte a2 = (b & (1 << i)) >> i;
    res ^= a1 & a2;
  }
  return res;
}

int compute_L_of(const sbox s, byte a, byte b) {
  int res = 0;
  for (int x = 0; x < 16; x++) {
    if (strange_op(a, x) == strange_op(b, s[x]))
      res += 1;
  }
  return res;
}

typedef vector<vector<double>> matrix;

matrix populate_linapp_matrix(const sbox s) {
  matrix m;
  for (int i = 0; i < 16; i++) {
    m.push_back(vector<double>());
    for (int j = 0; j < 16; j++) {
      double mij = (double)compute_L_of(s, i, j); // / 16.0;
      if (mij == 14 || mij == 2) cout << mij << " " << i << " " << j << endl;
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

/* ==== Question 4 ==== */
typedef std::mt19937 Rng;
Rng rng;
std::uniform_int_distribution<block> block_dist;

block random_block() {
  return block_dist(rng);
}

void pb(block b) {
  cout << bitset<32>(b) << endl;
}

bool check_rel(block K0, block K1, block msg, block A, block B) {
  bool x = strange_op(rotr(B, 2), B32_turn(K1, msg ^ K0));
  // cout << (x == strange_op(B, 0xF0000000 & apply_subst(S, msg ^ K0) ^ (rotr(K1, 30)))) << endl;
  return x == strange_op(A, msg);
}

double question4_aux(byte a, byte b) {
  assert (a < 16 && b < 16);
  block K0 = random_block();
  block K1 = random_block();
  const int N = 2000;
  block A = a << 28;
  block B = b << 28;
  int total = 0;
  for (int i = 0; i < N; i++) {
    block msg = random_block();
    if (check_rel(K0, K1, msg, A, B)) {
      total += 1;
    }
  }
  return double(total) / double(N);
}

double question4(byte a, byte b) {
  const int N = 100;
  double total = 0;
  for (int i = 0; i < N; i++) {
      auto t = question4_aux(a, b);
      if (t > .5)
	  t = 1.0 - t;
      total += t;
  }
  return total / double(N);
}

bool check_rel2(block K0, block K1, block msg,block B) {
  bool x = strange_op(rotr(B, 2), B32_turn(K1, msg ^ K0));
  // cout << (x == strange_op(B, 0xF0000000 & apply_subst(S, msg ^ K0) ^ (rotr(K1, 30)))) << endl;
  return  strange_op(B, apply_subst(S, msg ^ K0) ^ K1)== strange_op(B, apply_subst(S, msg));
}

double experiment(byte a, byte b) {
    assert (a < 16 && b < 16);
    block K0 = random_block();
    block K1 = random_block();
    const int N = 2000;
    block B = b << 28;
    int total = 0;
    for (int i = 0; i < N; i++) {
	block msg = random_block();
	if (check_rel2(K0, K1, msg, B)) {
	    total += 1;
	}
    }
    return double(total) / double(N);
}

/* ==== main ==== */

int main() {
  rng.seed(time(0));

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

  cout << endl;

  sbox SS = {12,5,6,11,9,0,10,13,3,14,15,8,4,7,1,2};
  matrix m = populate_linapp_matrix(S);
  print_matrix(m);

  cout << endl;

  /*for (int i = 0; i < 16; i++)
      for (int j = 0; j < 16; j++)
	  cout << i << "," << j << ": " << question4(i, j) << endl;
	  cout << question4(5, 1) << endl; */
  cout << question4(1, 5) << endl;

  int x = 12;
  cout << strange_op(11, 9) << " " << strange_op(11<<x, 9<<x);

  return 0;
}
