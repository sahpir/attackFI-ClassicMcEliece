#ifndef crypto_kem_H
#define crypto_kem_H

#include "crypto_kem_mceliece6688128.h"

#define crypto_kem_keypair crypto_kem_mceliece6688128_keypair
#define crypto_kem_enc_modified crypto_kem_mceliece6688128_enc
#define crypto_kem_dec_faulty crypto_kem_mceliece6688128_dec
#define crypto_kem_PUBLICKEYBYTES crypto_kem_mceliece6688128_PUBLICKEYBYTES
#define crypto_kem_SECRETKEYBYTES crypto_kem_mceliece6688128_SECRETKEYBYTES
#define crypto_kem_BYTES crypto_kem_mceliece6688128_BYTES
#define crypto_kem_CIPHERTEXTBYTES crypto_kem_mceliece6688128_CIPHERTEXTBYTES
#define crypto_kem_PRIMITIVE "mceliece6688128"

#endif
