#ifndef NTT_H
#define NTT_H

#include <stdint.h>
#include "params.h"

#define zetas KYBER_NAMESPACE(zetas)
extern const int16_t zetas[128];

#define ntt KYBER_NAMESPACE(ntt)
void ntt(int16_t r[KYBER_N]);

#define invntt KYBER_NAMESPACE(invntt)
void invntt(int16_t r[KYBER_N]);

#define basemul KYBER_NAMESPACE(basemul)
void basemul(int16_t r[2], const int16_t a[2], const int16_t b[2], int16_t zeta);

//#define fqmul KYBER_NAMESPACE(fqmul)
//static int16_t fqmul(int16_t a, int16_t b);

#endif
