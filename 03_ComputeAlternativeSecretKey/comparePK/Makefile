#output.txt: comparePK
#	./run

comparePK_3488: Makefile main.c gf_3488.c pk_gen_mod.c root.c
	gcc -O3 -march=native -mtune=native -Wall -DGFBITS=12 -DSYS_N=3488 -DSYS_T=64 -Dcrypto_kem_PUBLICKEYBYTES=261120 -Dcrypto_kem_SECRETKEYBYTES=6492 -Dcrypto_kem_CIPHERTEXTBYTES=128 "-DCRYPTO_NAMESPACE(x)=x" "-D_CRYPTO_NAMESPACE(x)=_##x" -o comparePK_3488 main.c gf_3488.c pk_gen_mod.c root.c 
comparePK_4608: Makefile main.c gf_4608.c pk_gen_mod.c root.c
	gcc -O3 -march=native -mtune=native -Wall -DGFBITS=13 -DSYS_N=4608 -DSYS_T=96 -Dcrypto_kem_PUBLICKEYBYTES=524160 -Dcrypto_kem_SECRETKEYBYTES=13608 -Dcrypto_kem_CIPHERTEXTBYTES=188 "-DCRYPTO_NAMESPACE(x)=x" "-D_CRYPTO_NAMESPACE(x)=_##x" -o comparePK_4608 main.c gf_4608.c pk_gen_mod.c root.c
comparePK_6688: Makefile main.c gf_6688.c pk_gen_mod.c root.c
	gcc -O3 -march=native -mtune=native -Wall -DGFBITS=13 -DSYS_N=6688 -DSYS_T=128 -Dcrypto_kem_PUBLICKEYBYTES=1044992 -Dcrypto_kem_SECRETKEYBYTES=13932 -Dcrypto_kem_CIPHERTEXTBYTES=240 "-DCRYPTO_NAMESPACE(x)=x" "-D_CRYPTO_NAMESPACE(x)=_##x" -o comparePK_6688 main.c gf_6688.c pk_gen_mod.c root.c
comparePK_6960: Makefile main.c gf_6960.c pk_gen_mod.c root.c
	gcc -O3 -march=native -mtune=native -Wall -DGFBITS=13 -DSYS_N=6960 -DSYS_T=119 -Dcrypto_kem_PUBLICKEYBYTES=1047319 -Dcrypto_kem_SECRETKEYBYTES=13948 -Dcrypto_kem_CIPHERTEXTBYTES=226 "-DCRYPTO_NAMESPACE(x)=x" "-D_CRYPTO_NAMESPACE(x)=_##x" -o comparePK_6960 main.c gf_6960.c pk_gen_mod.c root.c
comparePK_8192: Makefile main.c gf_8192.c pk_gen_mod.c root.c
	gcc -O3 -march=native -mtune=native -Wall -DGFBITS=13 -DSYS_N=8192 -DSYS_T=128 -Dcrypto_kem_PUBLICKEYBYTES=1357824 -Dcrypto_kem_SECRETKEYBYTES=14120 -Dcrypto_kem_CIPHERTEXTBYTES=240 "-DCRYPTO_NAMESPACE(x)=x" "-D_CRYPTO_NAMESPACE(x)=_##x" -o comparePK_8192 main.c gf_8192.c pk_gen_mod.c root.c

clean: 
	rm -f *.o comparePK_3488 comparePK_4096 comparePK_6688 comparePK_6960 comparePK_8192
