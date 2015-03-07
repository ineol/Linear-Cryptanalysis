#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cassert>
#include <random>
#include <time.h>
#include <algorithm>
#include <iterator>
#include <unordered_map>

using namespace std;

typedef unsigned char byte;
typedef vector<byte> bytes;
typedef uint32_t block;

/* === Block cipher === */

typedef byte sbox[16];

const sbox S = { 7, 3, 6, 1, 13, 9, 10, 11, 2, 12, 0, 4, 5, 15, 8, 14 };
const sbox Sinv = { 10, 3, 8, 1, 11, 12, 2, 0, 14, 5, 6, 7, 9, 4, 15, 13 };

block apply_subst(const sbox s, block b)
{
    block res = 0;
    for (int i = 0; i < 8; i++) { // From right to left
        block shifted = b >> 4 * i;
        byte cut = shifted & 0xF;
        block img = s[cut];
        res = res | (img << 4 * i);
    }
    return res;
}

block rotr(block x, byte n)
{
    return (x >> n) | (x << (32 - n)); // compiled into a rorl instruction
}

block B32_turn(block K, block b) { return rotr(apply_subst(S, b), 2) ^ K; }

block B32_inv_turn(block K, block b)
{
    return apply_subst(Sinv, rotr((b ^ K), 30));
}

struct B32_Cipher {
    block K0;
    block K1;
    block K2;
};

block B32_Encode(B32_Cipher c, block b)
{
    block xored = b ^ c.K0;
    block round1 = B32_turn(c.K1, xored);
    return B32_turn(c.K2, round1);
}

block B32_Decode(B32_Cipher c, block b)
{
    block round2 = B32_inv_turn(c.K2, b);
    block round1 = B32_inv_turn(c.K1, round2);
    return round1 ^ c.K0;
}

/* === Matrix of linear approximations === */
int strange_op(block a, block b)
{
    int res = 0;
    for (int i = 0; i < 32; i++) {
        byte a1 = (a & (1 << i)) >> i;
        byte a2 = (b & (1 << i)) >> i;
        res ^= a1 & a2;
    }
    return res;
}

int compute_L_of(const sbox s, byte a, byte b)
{
    int res = 0;
    for (int x = 0; x < 16; x++) {
        if (strange_op(a, x) == strange_op(b, s[x]))
            res += 1;
    }
    return res;
}

typedef vector<vector<double> > matrix;
typedef vector<vector<int> > imatrix;

matrix populate_linapp_matrix(const sbox s)
{
    matrix m;
    for (int i = 0; i < 16; i++) {
        m.push_back(vector<double>());
        for (int j = 0; j < 16; j++) {
            double mij = (double)compute_L_of(s, i, j) / 16.0;
            m[i].push_back(mij);
        }
    }
    return m;
}

imatrix populate_L(const sbox s)
{
    imatrix m;
    for (int i = 0; i < 16; i++) {
        m.push_back(vector<int>());
        for (int j = 0; j < 16; j++) {
            int mij = compute_L_of(s, i, j);
            m[i].push_back(mij);
        }
    }
    return m;
}

void print_matrix(matrix m)
{
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

block random_block() { return block_dist(rng); }

void pb(block b) { cout << bitset<32>(b) << endl; }

bool check_rel(block K0, block K1, block msg, block A, block B)
{
    bool x = strange_op(rotr(B, 2), B32_turn(K1, msg ^ K0));
    // cout << (x == strange_op(B, 0xF0000000 & apply_subst(S, msg ^ K0) ^
    // (rotr(K1, 30)))) << endl;
    return x == strange_op(A, msg);
}

double question4_aux(byte a, byte b)
{
    assert(a < 16 && b < 16);
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

double question4(byte a, byte b)
{
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

bool check_rel2(block K0, block K1, block msg, block B)
{
    return strange_op(B, apply_subst(S, msg ^ K0) ^ K1) == strange_op(B, apply_subst(S, msg));
}

double experiment(byte a, byte b)
{
    assert(a < 16 && b < 16);
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

/* ==== Question 5 ==== */

/* Returns 0 if the active box is the one bellow the position of a
 *         1 ---------------------------- after ------------------
 *        -1 if there are two
 * a must be a 4-bit integer
 */
int active_sbox(byte a)
{
    assert(a <= 0xF);
    if ((a & 0xC) == 0) { // 1100
        return 1;
    }
    else if ((a & 0x3) == 0) { // 0011
        return 0;
    }
    else {
        return -1;
    }
}

typedef pair<byte, byte> couple;

/* Finds the couples (a, b) of question 3 */
vector<couple> interresting_couples(const imatrix L)
{
    vector<couple> res;
    int N = L.size();
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            auto Lij = L[i][j];
            if (Lij == 14 || Lij == 2) {
                res.push_back(couple(i, j));
            }
        }
    }
    return res;
}

template <class X, class Y>
basic_ostream<X, Y>& operator<<(basic_ostream<X, Y>& out, const couple& c)
{
    out << "(" << (int)c.first << ", " << (int)c.second << ")";
    return out;
}

ostream& operator<<(ostream& out, const vector<couple> cs)
{
    out << "[";
    for (auto p : cs) {
        out << p << ", ";
    }
    out << "]";
    return out;
}

void print_linear_eq(block X, char var)
{
    bool first = true;
    for (int i = 0; i < 32; i++) {
        block bit = (X & (1 << i)) >> i;
        assert(bit == 1 || bit == 0); // duh
        if (bit == 1) {
            if (!first) {
                cout << " + ";
            }
            else {
                first = false;
            }
            cout << var << i;
        }
    }
}

void question6(const imatrix& L)
{
    for (auto& p : interresting_couples(L)) {
        byte a = p.first;
        byte b = p.second;
        block A = a << 28;
        block B = b << 28;

        cout << "===== " << p << " =====" << endl;
        cout << "P(B) = ";
        pb(rotr(B, 2));
        print_linear_eq(A, 'M');
        cout << " = ";
        print_linear_eq(rotr(B, 2), 'X');
        cout << "\n";
        int boxes = active_sbox(b);
        cout << "Active S-boxes: ";
        if (boxes == 0) {
            cout << "0\n";
        }
        else if (boxes == 1) {
            cout << "1\n";
        }
        else if (boxes == -1) {
            cout << "0 and 1\n";
        }
        else {
            assert(false);
        }
    }
}

/* === Question 7 === */
#include "known_ciphertexts_test.cc" // Booooooohhhh

/* c : cipher, k2: key guess
 */
block x1_guess(block c, block k2)
{
    block x2 = c ^ k2;
    block before_S = rotr(x2, 30);
    return apply_subst(Sinv, before_S);
}

double linear_rel_stats(byte a, byte b, block k2, int shift)
{
    block A = a << (28 - 4 * shift);
    block B = b << (28 - 4 * shift);

    int N = Plaintext.size();
    double total = 0;
    for (int i = 0; i < N; i++) {
        block x1 = x1_guess(Ciphertext[i], k2);
        if (strange_op(A, Plaintext[i]) == strange_op(rotr(B, 2), x1)) {
            total += 1;
        }
    }
    double proba = total / double(N);
    return (proba > .5) ? 1.0 - proba : proba;
}

typedef vector<double> bvec;

bvec populate_map_1_box(byte a, byte b, int shift)
{
    bvec res(16);
    assert(active_sbox(b) == 0); // Always the case
    for (int i = 0; i <= 0xF; i++) {
        block guess = i << (28 - 4 * shift);
        guess = rotr(guess, 2);
        res[i] = linear_rel_stats(a, b, guess, shift);
    }
    return res;
}

bvec average(bvec m1, bvec m2)
{
    bvec m(m1.size());
    for (int i = 0; i < m.size(); i++) {
        m[i] = (m1[i] + m2[i]) / 2.0;
    }
    return m;
}

byte find_subkey_1_box(vector<pair<byte, byte> > as, int shift)
{
    bvec res = populate_map_1_box(as[0].first, as[0].second, shift);
    for (int i = 1; i < as.size(); i++) {
        res = average(res, populate_map_1_box(as[i].first, as[i].second, shift));
    }

    return distance(res.begin(), min_element(res.begin(), res.end()));
}

block find_key_1_block(vector<pair<byte, byte> > as)
{
    block res = 0;
    for (int i = 0; i < 8; i++) {
        res |= (find_subkey_1_box(as, i) << (28 - 4 * i));
    }
    return rotr(res, 2);
}

/* ==== Question 9 ==== */

typedef bitset<32> bskey;

const vector<byte> k0_of_k = { 17, 31, 0, 0, 18, 7, 20, 18, 8, 1, 27, 27, 2, 4, 11, 20, 25, 13, 17, 10, 24, 9, 29, 15, 21, 18, 28, 20, 4, 5, 24, 15 };
const vector<byte> k1_of_k = { 15, 2, 5, 0, 13, 31, 5, 10, 18, 2, 3, 14, 14, 0, 11, 1, 20, 15, 14, 27, 6, 11, 19, 3, 6, 20, 14, 2, 28, 11, 5, 8 };
const vector<byte> k2_of_k = { 4, 24, 23, 12, 22, 21, 31, 15, 29, 1, 0, 26, 17, 24, 16, 5, 31, 0, 20, 21, 26, 30, 15, 11, 16, 23, 18, 30, 30, 19, 28, 23 };

/* Note: mk = master key = K */

block create_key_from_mk(bskey K, vector<byte> r)
{
    bskey res;
    for (int i = 0; i < r.size(); i++) {
        res[i] = K[r[i]];
    }
    unsigned long x = res.to_ulong();
    assert(x < 1UL << 32);
    return (block)x;
}

unordered_map<byte, byte> inverse_key_schedule(vector<byte> r)
{
    unordered_map<byte, byte> m;
    for (int i = 0; i < r.size(); i++) {
        m[r[i]] = i;
    }
    return m;
}

bskey create_mk_template_from_k2(block K2)
{
    bskey mk;
    bskey k2(K2);
    for (auto p : inverse_key_schedule(k2_of_k)) {
        mk[p.first] = k2[p.second];
    }
    return mk;
}

/* The bits of the master key we need to brute force are : */
const vector<byte> to_bruteforce = { 2, 3, 6, 7, 8, 9, 10, 13, 14, 25, 27 };

block complete_mk_with(bskey mk, block seed)
{
    assert(seed < (1 << to_bruteforce.size()));
    for (int i = 0; i < to_bruteforce.size(); i++) {
        byte ith_bit = (seed & (1 << i)) >> i;
        assert(ith_bit == 0 || ith_bit == 1);
        mk[to_bruteforce[i]] = ith_bit;
    }

    unsigned long x = mk.to_ulong();
    assert(x < 1UL << 32);
    return (block)x;
}

block find_mk(block K2)
{
    bskey mkt = create_mk_template_from_k2(K2);
    for (int i = 0; i < (1 << to_bruteforce.size()); i++) {
        block guess = complete_mk_with(mkt, i);

        B32_Cipher c = {
            create_key_from_mk(guess, k0_of_k),
            create_key_from_mk(guess, k1_of_k),
            create_key_from_mk(guess, k2_of_k)
        };
        assert(c.K2 == K2);
        for (int j = 0; j < Ciphertext.size(); j++) {
            if (Ciphertext[j] != B32_Encode(c, Plaintext[j])) {
                goto next_guess; // Take that Dijkstra
            }
        }
        return guess;

    next_guess:
        ;
    }
    assert(false && "Could not find master key!");
}

/* ==== main ==== */

int main()
{
    assert(Plaintext.size() == Ciphertext.size());

    rng.seed(time(0));

    block k = 0u | 1u | (1u << 31);
    B32_Cipher test_key = { k, ~0u, ~k };

    cout << bitset<32>(test_key.K0) << endl;
    cout << bitset<32>(test_key.K1) << endl;
    cout << bitset<32>(test_key.K2) << endl;
    cout << bitset<32>(B32_Encode(test_key, 0)) << endl;
    cout << bitset<32>(B32_Decode(test_key, B32_Encode(test_key, 0))) << endl;

    cout << endl;

    matrix m = populate_linapp_matrix(S);
    print_matrix(m);

    cout << endl;

    cout << question4(1, 5) << endl;

    imatrix L = populate_L(S);
    auto couples = interresting_couples(L);
    cout << couples << endl;

    question6(L);

    block MK = random_block();
    TEST_MK = MK;

    TEST_CIPHER = {
        create_key_from_mk(MK, k0_of_k),
        create_key_from_mk(MK, k1_of_k),
        create_key_from_mk(MK, k2_of_k)
    };
    init_test();

    byte guess = find_subkey_1_box({ { 4, 8 }, { 9, 4 }, { 13, 12 } }, 0);
    cout << (int)guess << endl;
    pb(guess);
    block k2 = find_key_1_block({ { 4, 8 }, { 9, 4 }, { 13, 12 } });
    pb(k2);
    pb(TEST_CIPHER.K2);
    pb(MK);
    pb(find_mk(TEST_CIPHER.K2));

    return 0;
}
