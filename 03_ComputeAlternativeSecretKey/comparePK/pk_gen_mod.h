/*
  This file is for public-key generation
*/

#ifndef PK_GEN_H
#define PK_GEN_H
#define pk_gen_mod CRYPTO_NAMESPACE(pk_gen_mod)

#include "gf.h"

int pk_gen_mod(unsigned char *, gf*,gf*);

#endif

