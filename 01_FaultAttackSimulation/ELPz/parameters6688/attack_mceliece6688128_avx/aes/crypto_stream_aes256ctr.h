#ifndef crypto_stream_aes256ctr_H
#define crypto_stream_aes256ctr_H

#define crypto_stream_aes256ctr crypto_stream
#define crypto_stream_aes256ctr_xor crypto_stream_xor


#ifdef __cplusplus
#include <string>
extern std::string crypto_stream_aes256ctr(size_t,const std::string &,const std::string &) __attribute__((visibility("default")));
extern std::string crypto_stream_aes256ctr_xor(const std::string &,const std::string &,const std::string &) __attribute__((visibility("default")));
extern "C" {
#endif

extern int crypto_stream_aes256ctr(unsigned char *,unsigned long long,const unsigned char *,const unsigned char *) __attribute__((visibility("default")));
extern int crypto_stream_aes256ctr_xor(unsigned char *,const unsigned char *,unsigned long long,const unsigned char *,const unsigned char *) __attribute__((visibility("default")));
#ifdef __cplusplus
}
#endif

#endif
