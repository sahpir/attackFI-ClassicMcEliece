bruteForce_6688: Makefile main.c
	gcc -O3 -march=native -mtune=native -Wall -DSYS_N=6688 -Dcrypto_kem_CIPHERTEXTBYTES=240 -DSYND_BYTES=208 -o bruteForce_6688 main.c keccak/SimpleFIPS202.c keccak/KeccakSponge.c keccak/KeccakP-1600-opt64.c

clean: 
	rm -f *.o bruteForce_6688
