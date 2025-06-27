#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "kem.h"
#include "kex.h"
#include "params.h"
#include "indcpa.h"
#include "poly.h"
#include "polyvec.h"
#include "cpucycles.h"
#include "ntt.h"
#include "randombytes.h"
#include "symmetric.h"
#include "indcpa.h"
#include <time.h>

/*
static 
void VectorVectorMul(poly *mp, polyvec *b, polyvec *skpv)
{
  polyvec_ntt(b);
  polyvec_basemul_acc_montgomery(mp, skpv, b);
  poly_invntt_tomont(mp);
}

static 
void MatrixVectorMul(polyvec at[KYBER_K], polyvec *sp, polyvec *b)
{
  polyvec_ntt(sp);
  // matrix-vector multiplication
  for(int i=0;i<KYBER_K;i++)
    polyvec_basemul_acc_montgomery(&b->vec[i], &at[i], sp);

  polyvec_invntt_tomont(b);
}
*/


// micro second 
#define NTESTS 100000

#define TIME(s,name) clock_gettime(CLOCK_MONOTONIC_RAW, &s);
// Result is nanosecond per call 
#define  CALC(start, stop) \
  ((double) ((stop.tv_sec - start.tv_sec) * 1000000000 + (stop.tv_nsec - start.tv_nsec))) / NTESTS;


static 
void VectorVectorMul(poly *mp, polyvec *b, polyvec *skpv)
{
  polyvec_ntt(b);
  polyvec_basemul_acc_montgomery(mp, skpv, b);
  poly_invntt_tomont(mp);
}

static 
void MatrixVectorMul(polyvec at[KYBER_K], polyvec *sp, polyvec *b)
{
  polyvec_ntt(sp);
  // matrix-vector multiplication
  for(int i=0;i<KYBER_K;i++)
    polyvec_basemul_acc_montgomery(&b->vec[i], &at[i], sp);

  polyvec_invntt_tomont(b);
}

uint8_t seed[KYBER_SYMBYTES] = {0};

int main()
{
  unsigned int i;
  unsigned char pk[CRYPTO_PUBLICKEYBYTES] = {0};
  unsigned char sk[CRYPTO_SECRETKEYBYTES] = {0};
  unsigned char ct[CRYPTO_CIPHERTEXTBYTES] = {0};
  unsigned char key[CRYPTO_BYTES] = {0};
  // unsigned char kexsenda[KEX_AKE_SENDABYTES] = {0};
  // unsigned char kexsendb[KEX_AKE_SENDBBYTES] = {0};
  // unsigned char kexkey[KEX_SSBYTES] = {0};
  uint8_t ss[32]={0};
  uint8_t kr[64]={0};
  uint8_t c[KYBER_INDCPA_BYTES];
 
  unsigned char msg[KYBER_INDCPA_MSGBYTES] = {0};
  polyvec matrix[KYBER_K];
  poly ap;
  polyvec sp, b, skpv, e;
  uint8_t buf[2*KYBER_SYMBYTES]={0};
  const uint8_t *publicseed = buf;
  // poly bp;

  struct timespec start, stop;
  long ns;


  TIME(start, "gen_matrix");
  for(i=0;i<NTESTS;i++) {
    gen_matrix(matrix, seed, 0);
  }
  TIME(stop, "gen_matrix");
  ns = CALC(start, stop);
  printf("gen_matrix: %lld\n", ns);

  TIME(start, "poly_getnoise_eta1");
  for(i=0;i<NTESTS;i++) {
    poly_getnoise_eta1(&ap, seed, 1);
  }
  TIME(stop, "poly_getnoise_eta1");
  ns = CALC(start, stop);
  printf("poly_getnoise_eta1: %lld\n", ns);

  TIME(start, "poly_getnoise_eta2");
  for(i=0;i<NTESTS;i++) {
    poly_getnoise_eta2(&ap, seed, 0);
  }
  TIME(stop, "poly_getnoise_eta2");
  ns = CALC(start, stop);
  printf("poly_getnoise_eta2: %lld\n", ns);

  TIME(start, "poly_tomsg");
  for(i=0;i<NTESTS;i++) {
    poly_tomsg(msg, &ap);
  }
  TIME(stop, "poly_tomsg");
  ns = CALC(start, stop);
  printf("poly_tomsg: %lld\n", ns);

  TIME(start, "poly_frommsg");
  for(i=0;i<NTESTS;i++) {
    poly_frommsg(&ap, msg);
  }
  TIME(stop, "poly_frommsg");
  ns = CALC(start, stop);
  printf("poly_frommsg: %lld\n", ns);


  TIME(start, "ref_ntt");
  for(i=0;i<NTESTS;i++) {
    ntt(ap.coeffs);
  }
  TIME(stop, "ref_ntt");
  ns = CALC(start, stop);
  printf("ref_ntt: %lld\n", ns);

  TIME(start, "ref_invntt");
  for(i=0;i<NTESTS;i++) {
    invntt(ap.coeffs);
  }
  TIME(stop, "ref_invntt");
  ns = CALC(start, stop);
  printf("ref_invntt: %lld\n", ns);
  /********************************/
  printf("*************Key pair Gen Results***********************\n");
  TIME(start, "stime");
  for(i=0; i<NTESTS; i++){
    randombytes(buf, KYBER_SYMBYTES);
  }
  TIME(stop, "etime");
  ns=CALC(start, stop);
  printf("random_bytes: %lld\n", ns);

  TIME(start, "stime");
  for(i=0; i<NTESTS; i++){
    hash_g(buf, buf, KYBER_SYMBYTES);
  }
  TIME(stop, "etime");
  ns=CALC(start, stop);
  printf("hash_g SHA3-512: %lld\n", ns);

  TIME(start, "stime");
  for(i=0; i<NTESTS; i++){
    polyvec_ntt(&skpv);
  }
  TIME(stop, "etime");
  ns=CALC(start, stop);
  printf("polyvec_ntt: %lld\n", ns);

TIME(start, "stime");
  for(i=0; i<NTESTS; i++){
    polyvec_basemul_acc_montgomery(&ap, &matrix[0], &skpv);
  
  }
  TIME(stop, "etime");
  ns=CALC(start, stop);
  printf("polyvec_basemul_acc_montgomery: %lld\n", ns);

  TIME(start, "stime");
  for(i=0; i<NTESTS; i++){
    poly_tomont(&ap);
  }
  TIME(stop, "etime");
  ns=CALC(start, stop);
  printf("poly_tomont: %lld\n", ns);

  TIME(start, "stime");
  for(i=0; i<NTESTS; i++){
    polyvec_add(&skpv, &skpv, &e);
  }
  TIME(stop, "etime");
  ns=CALC(start, stop);
  printf("polyvec_add: %lld\n", ns);

  TIME(start, "stime");
  for(i=0; i<NTESTS; i++){
    polyvec_reduce(&skpv);
  }
  TIME(stop, "etime");
  ns=CALC(start, stop);
  printf("polyvec_reduce: %lld\n", ns);

  TIME(start, "stime");
  for(i=0; i<NTESTS; i++){
    pack_sk(sk, &skpv);
  }
  TIME(stop, "etime");
  ns=CALC(start, stop);
  printf("pack_sk: %lld\n", ns);

  TIME(start, "stime");
  for(i=0; i<NTESTS; i++){
    pack_pk(sk, &skpv, publicseed);
  }
  TIME(stop, "etime");
  ns=CALC(start, stop);
  printf("pack_pk: %lld\n", ns);

  /*******************************/

  TIME(start, "crypto_kem_keypair");
  for(i=0;i<NTESTS;i++) {
    crypto_kem_keypair(pk, sk);
  }
  TIME(stop, "crypto_kem_keypair");
  ns = CALC(start, stop);
  printf("crypto_kem_keypair_overall: %lld\n", ns);

printf("*******************Encryption Results*************************\n");
  TIME(start, "crypto_kem_enc");
  for(i=0;i<NTESTS;i++) {
    crypto_kem_enc(ct, key, pk);
  }
  TIME(stop, "crypto_kem_enc");
  ns = CALC(start, stop);
  printf("crypto_kem_enc: %lld\n", ns);
  
  TIME(start, "sha3-256");
  for(i=0;i<NTESTS;i++) {
    hash_h(buf, buf,KYBER_SYMBYTES);
  }
  TIME(stop, "sha3-256");
  ns = CALC(start, stop);
  printf("SHA3_256: %lld\n", ns);


  TIME(start, "poly_from_msg");
  for(i=0;i<NTESTS;i++) {
    poly_frommsg(&ap, msg);
  }
  TIME(stop, "poly_from_msg");
  ns = CALC(start, stop);
  printf("poly_from_msg: %lld\n", ns);


    TIME(start, "polyvec_invntt_tomont");
  for(i=0;i<NTESTS;i++) {
    polyvec_invntt_tomont(&sp);
  }
  TIME(stop, "polyvec_invntt_tomont");
  ns = CALC(start, stop);
  printf("polyvec_invntt_tomont: %lld\n", ns);

  TIME(start, "poly_invntt_tomont");
  for(i=0;i<NTESTS;i++) {
    poly_invntt_tomont(&ap);
  }
  TIME(stop, "poly_invntt_tomont");
  ns = CALC(start, stop);
  printf("poly_invntt_tomont: %lld\n", ns);


  TIME(start, "poly_add");
  for(i=0;i<NTESTS;i++) {
    poly_add(&ap, &ap, &ap);
  }
  TIME(stop, "poly_add");
  ns = CALC(start, stop);
  printf("poly_add: %lld\n", ns);


   TIME(start, "poly_reduce");
  for(i=0;i<NTESTS;i++) {
    poly_reduce(&ap);
  }
  TIME(stop, "poly_reduce");
  ns = CALC(start, stop);
  printf("poly_reduce: %lld\n", ns);


  TIME(start, "pack_ciphertext");
  for(i=0;i<NTESTS;i++) {
    pack_ciphertext(c, &sp, &ap);
  }
  TIME(stop, "pack_ciphertext");
  ns = CALC(start, stop);
  printf("pack_ciphertext: %lld\n", ns);


   TIME(start, "kdf");
  for(i=0;i<NTESTS;i++) {
    kdf(ss,kr, 2*KYBER_SYMBYTES);
  }
  TIME(stop, "kdf");
  ns = CALC(start, stop);
  printf("KDF: %lld\n", ns);

  printf("*******************Decryption Results****************************\n");

  TIME(start, "crypto_kem_dec");
  for(i=0;i<NTESTS;i++) {
    crypto_kem_dec(key, ct, sk);
  }
  TIME(stop, "crypto_kem_dec");
  ns = CALC(start, stop);
  printf("crypto_kem_dec: %lld\n", ns);



  TIME(start, "unpack_ciphertext");
  for(i=0;i<NTESTS;i++) {
    unpack_ciphertext(&sp, &ap, c);
  }
  TIME(stop, "unpack_ciphertext");
  ns = CALC(start, stop);
  printf("unpack_ciphertext: %lld\n", ns);

  TIME(start, "unpack_sk");
  for(i=0;i<NTESTS;i++) {
    unpack_sk(&skpv, sk);
  }
  TIME(stop, "unpack_sk");
  ns = CALC(start, stop);
  printf("unpack_sk: %lld\n", ns);


   TIME(start, "poly_sub");
  for(i=0;i<NTESTS;i++) {
    poly_sub(&ap, &ap, &ap);
  }
  TIME(stop, "poly_sub");
  ns = CALC(start, stop);
  printf("poly_sub: %lld\n", ns);


   TIME(start, "poly_tomsg");
  for(i=0;i<NTESTS;i++) {
    poly_tomsg(msg, &ap);
  }
  TIME(stop, "poly_tomsg");
  ns = CALC(start, stop);
  printf("poly_tomsg: %lld\n", ns);

  printf("**********************************************************\n");

  TIME(start, "VectorVectorMul");
  for(i=0;i<NTESTS;i++) {
    VectorVectorMul(&ap, &sp, &b);
  }
  TIME(stop, "VectorVectorMul");
  ns = CALC(start, stop);
  printf("VectorVectorMul: %lld\n", ns);

  
  TIME(start, "stime");
  for(i=0;i<NTESTS;i++) {
    MatrixVectorMul(matrix, &sp, &b);
  }
  TIME(stop, "etime");
  ns = CALC(start, stop);
  printf("MatrixVectorMul: %lld\n", ns);

  return 0;
}
