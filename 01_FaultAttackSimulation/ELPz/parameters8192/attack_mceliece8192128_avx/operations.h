#ifndef OPERATIONS_H
#define OPERATIONS_H

#include "crypto_kem.h"

int crypto_kem_enc_modified(
       unsigned char *c,
       unsigned char *key,
       const unsigned char *pk,
       unsigned char *two_e
);

int crypto_kem_dec_faulty(
       unsigned char *key,
       const unsigned char *c,
       const unsigned char *sk,
       const unsigned char fault_position, 
       short unsigned int fault_value
);

int crypto_kem_keypair
(
       unsigned char *pk,
       unsigned char *sk 
);

#endif

