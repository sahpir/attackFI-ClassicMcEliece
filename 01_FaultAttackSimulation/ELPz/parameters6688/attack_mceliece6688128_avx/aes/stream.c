#include <string.h>
#include "aes.h"
#include "crypto_stream_aes256ctr.h"

int crypto_stream(
  unsigned char *out,
  unsigned long long outlen,
  const unsigned char *n,
  const unsigned char *k
)
{
  unsigned char temp[outlen];
  memset(temp, 0, outlen);
  return crypto_stream_xor(out, temp, outlen, n, k);
}

int crypto_stream_xor(
  unsigned char *out,
  const unsigned char *in,
  unsigned long long inlen,
  const unsigned char *n,
  const unsigned char *k
)
{
  struct AES_ctx ctx;
  //EVP_CIPHER_CTX *x;
  //int ok;
  //int outl = 0;

  //x = EVP_CIPHER_CTX_new();
  //if (!x) return -111;
    memcpy(out,in,inlen); //copy 128 bits of plaintext to buffer for encryption
    AES_init_ctx_iv(&ctx, k, n);
    AES_CTR_xcrypt_buffer(&ctx, out, inlen);
  
  //ok = EVP_EncryptInit_ex(x,EVP_aes_256_ctr(),0,k,n);
  //if (ok == 1) ok = EVP_CIPHER_CTX_set_padding(x, 0);
  //if (ok == 1) ok = EVP_EncryptUpdate(x, out, &outl, in, inlen);
  //if (ok == 1) ok = EVP_EncryptFinal_ex(x, out, &outl);

  //EVP_CIPHER_CTX_free(x);
  return 1 ;//ok == 1 ? 0 : -111;
}
