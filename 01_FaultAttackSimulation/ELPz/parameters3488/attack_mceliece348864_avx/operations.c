#include "operations.h"

#include "controlbits.h"
#include "randombytes.h"
#include "crypto_hash.h"
#include "encrypt_modified.h"
#include "decrypt_faulty.h"
#include "params.h"
#include "sk_gen.h"
#include "pk_gen.h"
#include "util.h"

#include <stdint.h>
#include <string.h>

int crypto_kem_enc_modified(
       unsigned char *c,
       unsigned char *key,
       const unsigned char *pk,
	   unsigned char *two_e
)
{
	//unsigned char two_e[ 1 + SYS_N/8 ] = {2};
	unsigned char *e = two_e + 1;
	unsigned char one_ec[ 1 + SYS_N/8 + (SYND_BYTES + 32) ] = {1};

	//

	encrypt_modified(c, pk, e);

	crypto_hash_32b(c + SYND_BYTES, two_e, sizeof(two_e)); 

	memcpy(one_ec + 1, e, SYS_N/8);
	memcpy(one_ec + 1 + SYS_N/8, c, SYND_BYTES + 32);

	crypto_hash_32b(key, one_ec, sizeof(one_ec));

	return 0;
}

int crypto_kem_dec_faulty(
       unsigned char *key,
       const unsigned char *c,
       const unsigned char *sk,
	   const unsigned char fault_position, // specifies where the Eroor Locator Polynomial need to be attacked
	   short unsigned int  fault_value
)
{
	int i;

	unsigned char ret_confirm = 0;
	unsigned char ret_decrypt = 0;

	uint16_t m;

	unsigned char conf[32];
	unsigned char two_e[ 1 + SYS_N/8 ] = {2};
	unsigned char *e = two_e + 1;
	unsigned char preimage[ 1 + SYS_N/8 + (SYND_BYTES + 32) ];
	unsigned char *x = preimage;
	const unsigned char *s = sk + 40 + IRR_BYTES + COND_BYTES;

	//

	ret_decrypt = decrypt_faulty(e, sk + 40, c, fault_position,fault_value);

	crypto_hash_32b(conf, two_e, sizeof(two_e)); 

	for (i = 0; i < 32; i++) 
		ret_confirm |= conf[i] ^ c[SYND_BYTES + i];

	m = 0; // Fault here
	//m = ret_decrypt | ret_confirm;
	m -= 1;
	m >>= 8;

	*x++ = m & 1;
	for (i = 0; i < SYS_N/8; i++) 
		*x++ = (~m & s[i]) | (m & e[i]);

	for (i = 0; i < SYND_BYTES + 32; i++) 
		*x++ = c[i];

	crypto_hash_32b(key, preimage, sizeof(preimage)); 

	return 0;
}

int crypto_kem_keypair
(
       unsigned char *pk,
       unsigned char *sk 
)
{
	int i;
	unsigned char seed[ 33 ] = {64};
	unsigned char r[ SYS_N/8 + (1 << GFBITS)*sizeof(uint32_t) + SYS_T*2 + 32 ];
	unsigned char *rp, *skp;

	gf f[ SYS_T ]; // element in GF(2^mt)
	gf irr[ SYS_T ]; // Goppa polynomial
	uint32_t perm[ 1 << GFBITS ]; // random permutation as 32-bit integers
	int16_t pi[ 1 << GFBITS ]; // random permutation
	
	// added to output support elements in file
// 	gf L[ SYS_N ]; // support

	randombytes(seed+1, 32);

	while (1)
	{
		rp = &r[ sizeof(r)-32 ];
		skp = sk;

		// expanding and updating the seed

		shake(r, sizeof(r), seed, 33);
		memcpy(skp, seed+1, 32);
		skp += 32 + 8;
		memcpy(seed+1, &r[ sizeof(r)-32 ], 32);

		// generating irreducible polynomial

		rp -= sizeof(f); 

		for (i = 0; i < SYS_T; i++) 
			f[i] = load_gf(rp + i*2); 

		if (genpoly_gen(irr, f)) 
			continue;
		
// 		printf("\nirr=[");
// 		for (i = 0; i < SYS_T; i++)
// 			printf("%.4x,",irr[i]);
// 		printf("%.4x]",1);
// 		

		for (i = 0; i < SYS_T; i++)
			store_gf(skp + i*2, irr[i]);

		skp += IRR_BYTES;

		// generating permutation

		rp -= sizeof(perm);

		for (i = 0; i < (1 << GFBITS); i++) 
			perm[i] = load4(rp + i*4); 

		if (pk_gen(pk, skp - IRR_BYTES, perm, pi))
			continue;
		
		// sp: Extract support elements from pi
// 		for (i = 0; i < SYS_N;         i++) L[i] = bitrev(pi[i]);
	
// 		printf("\nL=[");
// 		for (i = 0; i < SYS_N; i++)
// 			printf("%.4x,",L[i]);
// 		printf("]");

		controlbitsfrompermutation(skp, pi, GFBITS, 1 << GFBITS);
		skp += COND_BYTES;

		// storing the random string s

		rp -= SYS_N/8;
		memcpy(skp, rp, SYS_N/8);

		// storing positions of the 32 pivots

		store8(sk + 32, 0xFFFFFFFF);

		break;
	}

	return 0;
}

