
vector<block> Plaintext;
vector<block> Ciphertext;

B32_Cipher TEST_CIPHER;
block TEST_MK;

void init_test() {
    B32_Cipher c = TEST_CIPHER;
    
    const int N = 1000;

    for (int i = 0; i < N; i++) {
	block m = random_block();
	block cph = B32_Encode(c, m);
	Plaintext.push_back(m);
	Ciphertext.push_back(cph);
	assert (B32_Decode(c, cph) == m);
    }
}
