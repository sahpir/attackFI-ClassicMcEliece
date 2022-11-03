/*
   Compare Public Keys 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "params.h"
#include "pk_gen_mod.h"
#include "gf.h"

#define KAT_FILE_OPEN_ERROR -1
#define KAT_CRYPTO_FAILURE -4
// Parameters for n=3488 t=64 m =12
//#define crypto_kem_PUBLICKEYBYTES 261120
//#define crypto_kem_SECRETKEYBYTES 6492
//#define crypto_kem_CIPHERTEXTBYTES 128
#define crypto_kem_BYTES 32


int
main(int argc, char **argv)
{
    FILE                *fp_pk, *fp_rec_L,*fp_rec_g;
    int                 ret_val;
    int i;
    unsigned char *pk = 0;
    unsigned char *rec_pk = 0;
    
    //fp_pk = fopen("061550_mceliece3488_pk.txt", "r");
    fp_pk = fopen(argv[1], "r");
    if (!fp_pk)
        return KAT_FILE_OPEN_ERROR;
    //fp_rec_g = fopen("output_sage_g.txt", "r");
    fp_rec_g = fopen(argv[2], "r");
    if (!fp_rec_g)
        return KAT_FILE_OPEN_ERROR;
    //fp_rec_L = fopen("output_sage_a.txt", "r");
    fp_rec_L = fopen(argv[3], "r");
    if (!fp_rec_L)
        return KAT_FILE_OPEN_ERROR;
    

    if (!pk) pk = malloc(crypto_kem_PUBLICKEYBYTES);
    if (!pk) abort();
    if (!rec_pk) rec_pk = malloc(crypto_kem_PUBLICKEYBYTES);
    if (!rec_pk) abort();
   
    gf g[ SYS_T+1 ]; // storage for recovered Goppa polynomial 
    gf L[ SYS_N ]; // recovered support elements
    
    printf("PK:\n");
    for(i=0;i<crypto_kem_PUBLICKEYBYTES;i++){
        fscanf(fp_pk, "%hhx,",&pk[i]);
        printf("%.2X",pk[i]);
    }
    printf("\ng:\n");
    for(i=0;i<SYS_T+1;i++){
        fscanf(fp_rec_g, "%hx,",&g[i]);
        printf("%.4x ",g[i]);
    }
    printf("\nL:\n");
    for(i=0;i<SYS_N;i++){
        fscanf(fp_rec_L, "%hx,",&L[i]);
        printf("%.4x ",L[i]);
    }
    if((ret_val=pk_gen_mod(rec_pk,g,L))!=0){
        printf("Cannot calculate public key! Returned <%d>\n", ret_val);
        return KAT_CRYPTO_FAILURE; 
    }
    printf("\nrec_pk:\n");
    for(i=0;i<crypto_kem_PUBLICKEYBYTES;i++){
        printf("%.2X",rec_pk[i]);
    }
    // Compare public key with recovered public key
    for(i=0;i<crypto_kem_PUBLICKEYBYTES;i++){
        if(pk[i]!=rec_pk[i]){
            printf("Public Keys not identical! First position <%d>\n",i);
            return KAT_CRYPTO_FAILURE;
        }
    }
    printf("\nSUCCESS: PUBLIC KEYS IDENTICAL \n");
        
    
    fclose(fp_pk);
    fclose(fp_rec_g);
    fclose(fp_rec_L);

    return 0;
}
