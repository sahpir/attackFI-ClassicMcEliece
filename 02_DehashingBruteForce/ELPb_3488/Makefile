bruteForce_3488: Makefile main.c
	gcc -O3 -march=native -mtune=native -Wall -DSYS_N=3488 -Dcrypto_kem_CIPHERTEXTBYTES=128 -DSYND_BYTES=96 -o bruteForce_3488 main.c keccak/SimpleFIPS202.c keccak/KeccakSponge.c keccak/KeccakP-1600-opt64.c

clean: 
	rm -f *.o bruteForce_3488
