#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nist/rng.h"
#include "crypto_kem.h"
#include "params.h"
#include "crypto_hash.h"
#include "decrypt_faulty.h"
#include "util.h"
#include "gf.h" // only used to check equations


#define FILE_OPEN_ERROR -1
#define COMPUTATION_ERROR -4
#define KEY_RECOVERY_ERROR -5


static void gen_e_modified(unsigned char *e, unsigned int e_weight, uint16_t *ind);
int search_zero_element(uint16_t*faulty_indizes, unsigned char* faulty_hash, unsigned char* ct);
int get_indizes_weight1(uint16_t* faulty_indizes, unsigned char * output_faulty_e, uint16_t zero_position);
void fprintBstr(FILE *fp, char *S, unsigned char *A, unsigned long long L);

unsigned char entropy_input[48];
unsigned char seed[KATNUM][48];



int main(){
    
    FILE *fp_rsp, *fp_sage_lin,*fp_sage_quad,*fp_simulation,*fp_pk, *fp_hashes_lin, *fp_hashes_quad;
    //FILE *fp_support; // used to check equations
    int ret_val;
    int i;
    unsigned char *ct = 0; // storage for ciphertext ct=(c0,c1)=(eH^T,Hash(2,e))
    unsigned char *ss = 0; // storage for session key ss=K=Hash(1,e,(ct))
    unsigned char *faulty_hash = 0; //storage for faulty hash output
    unsigned char *pk = 0; //stores public key
    unsigned char *sk = 0; //stores secret key
    
    unsigned int sum_injections_constant=0;
    unsigned int sum_equations_constant=0;
    unsigned int sum_injections_linear=0;
    unsigned int sum_equations_linear=0;
    
    unsigned char fault_degree;
    unsigned char fault_position;
    uint16_t fault_value=0;
    
    unsigned char two_e[ 1 + SYS_N/8 ] = {2}; // storage for (2,e)
    unsigned char *e = two_e + 1; // pointer to plaintext
    unsigned int e_weight; // weight of plaintext
    uint16_t indizes[2]; // indizes of plaintext
    unsigned char output_faulty_e[SYS_N/8]; // storage for faulty output
    uint16_t faulty_indizes[SYS_T]; // indizes of faulty output
    uint16_t faulty_weight; // weight of faulty output
    
    uint16_t zero_position=-1; // That's the position of the zero element in the support set (if available), -1 if not exists
    uint64_t found[105]; // works only if SYS_N=6688; stores which different support elements are already found
    unsigned int nr_found_variables=0; // counts how many different support elements are already found
    unsigned char found_all=0; // Indicates if all N support elements are found or not
    unsigned int nr_found_equations=0; // counts how many equations are gathered in total for the polynomial system
    unsigned int nr_found_constant_equations; // counts how many quations for constant injection are gathered
    unsigned int nr_loops_constant=0;  // counts how often the fault_value had to be changed for constant injection
    unsigned char linear_injection=1; // specifies if injections on the linear coefficient of ELP should be made or not
    unsigned int nr_injections_constant=0, nr_injections_linear=0; // counts how many injections were needed to gather the equations
    
    // Used to mark which variables are already found
    for(i=0;i<105;i++) // works only if SYS_N=6688
        found[i]=0;
    

    for (i=0; i<48; i++)
        entropy_input[i] = i;
    randombytes_init(entropy_input, NULL, 256);

    for (i=0; i<KATNUM; i++)
        randombytes(seed[i], 48);

    fp_rsp = fdopen(9, "w");
    if (!fp_rsp)
        return FILE_OPEN_ERROR;

    fprintf(fp_rsp, "# kem/%s\n\n", crypto_kem_PRIMITIVE);
    
    fp_simulation = fopen("simulation_results.txt", "w");
    if (!fp_simulation)
        return FILE_OPEN_ERROR;

    // Start of simulation
    for (int num=0; num<KATNUM; num++){
        printf("\nNUM: %d\n",num);
        fprintf(fp_simulation,"NUM: %d \n", num);
        
        if (!ct) ct = malloc(crypto_kem_CIPHERTEXTBYTES);
        if (!ct) abort();
        if (!ss) ss = malloc(crypto_kem_BYTES);
        if (!ss) abort();
        if (!faulty_hash) faulty_hash = malloc(crypto_kem_BYTES);
        if (!faulty_hash) abort();
        if (!pk) pk = malloc(crypto_kem_PUBLICKEYBYTES);
        if (!pk) abort();
        if (!sk) sk = malloc(crypto_kem_SECRETKEYBYTES);
        if (!sk) abort();

    
        fp_sage_lin = fopen("input_sage_lin.sage", "w"); // File storing the linear equations
        if (!fp_sage_lin)
            return FILE_OPEN_ERROR;
        fp_sage_quad = fopen("input_sage_quad.sage", "w"); // File storing the quadratic equations
        if (!fp_sage_quad)
            return FILE_OPEN_ERROR;
        fp_pk = fopen("input_sage_pk_hex.txt", "w"); // File containing public key in hex format
        if (!fp_pk)
            return FILE_OPEN_ERROR;
        fp_hashes_lin = fopen("input_hashes_lin.txt", "w"); // File containing the hash output for constant injections/linear equations
        if (!fp_hashes_lin)
            return FILE_OPEN_ERROR;
        fp_hashes_quad = fopen("input_hashes_quad.txt", "w"); // File containing the hash output for linear injections/quadratic equations
        if (!fp_hashes_quad)
            return FILE_OPEN_ERROR;
        
        randombytes_init(seed[num], NULL, 256);

        fprintf(fp_rsp, "count = %d\n", num);
        fprintBstr(fp_rsp, "seed = ", seed[num], 48);
        
        if ( (ret_val = crypto_kem_keypair(pk, sk)) != 0) { // Generates public/secret keypair
            fprintf(stderr, "crypto_kem_keypair returned <%d>\n", ret_val);
            return COMPUTATION_ERROR;
        }
        
        fprintBstr(fp_rsp, "pk = ", pk, crypto_kem_PUBLICKEYBYTES);
        fprintBstr(fp_rsp, "sk = ", sk, crypto_kem_SECRETKEYBYTES);
        
        for(i=0;i<crypto_kem_PUBLICKEYBYTES;i++){  //output in hex
            fprintf(fp_pk,"%.2X,",pk[i]);
        }
        printf("\n");
        fclose(fp_pk);
    
    
        /******** CHECK EQUATIONS WITH SECRET KEY *************/
        /*fp_support = fopen("64335BF2_support.txt", "r");
        if (!fp_support)
            return FILE_OPEN_ERROR;
        
        // To check equations
        gf g[ SYS_T+1 ]; // storage for Goppa polynomial
        gf L[ SYS_N ]; // support element

        g[ SYS_T ] = 1;
        for (i = 0; i < SYS_T; i++) { 
            g[i] = load_gf(sk+2*i);
        }
        printf("\ng in kem_attack: \n");
        for (i = 0; i < SYS_T+1; i++)
            printf("%.4x ",g[i]);
        
        for(i=0;i<SYS_N;i++){
            fscanf(fp_support, "%hx,",&L[i]);
        }
        fclose(fp_support);
        
       printf("\nL in support: \n");
        for (i = 0; i < SYS_N; i++)
            printf("%.4x ",L[i]);*/
        /******** END: CHECK EQUATIONS WITH SECRET KEY *************/
        
        
        /******* SEARCH FOR ZERO ELEMENT ***************/
        
        // set e to zero
        memset(e, 0x00, SYS_N/8);
        
        if ( (ret_val = crypto_kem_enc_modified(ct, ss, pk,two_e)) != 0) { // Encrypt the error vector e, ss=K=Hash(1,e,C), ct=C=(c_0,Hash(2,e))
                fprintf(stderr, "crypto_kem_enc returned <%d>\n", ret_val);
                return COMPUTATION_ERROR;
            }
        fault_position=-1;
        if ( (ret_val = crypto_kem_dec_faulty(faulty_hash, ct, sk, fault_position, fault_value)) != 0) {
            fprintf(stderr, "crypto_kem_dec returned <%d>\n", ret_val);
            return COMPUTATION_ERROR;
        }
        
        //If no zero element then faulty_hash=Hash(1,0,(0,Hash(2,0)))
        faulty_weight=search_zero_element(faulty_indizes,faulty_hash,ct);
        if(faulty_weight!=0){
            zero_position=faulty_indizes[0];
        }
            
        if(zero_position!=((uint16_t)-1)){
#if DEBUG_LEVEL >0
            printf("Zero-element at position: %d \n", zero_position);
#endif
            fprintf(fp_simulation,"Zero-element at position: %d\n",zero_position);
            fprintf(fp_sage_lin, "%d\n",zero_position);
            // Set found variable for zero position
            found[zero_position/64]=found[zero_position/64]|((uint64_t)1<<(zero_position%64));
        }else
            fprintf(fp_simulation,"No zero-element in support set.\n");
        
        // Set ciphertext to zero
        memset(ct, 0x00, crypto_kem_CIPHERTEXTBYTES);
        // set e to zero
        memset(e, 0x00, SYS_N/8);
        
        /******* END: SEARCH FOR ZERO ELEMENT ***************/
        
        /************* CONSTANT INJECTIONS - wt(e)=2 ***************/
        // Run through possible error vectors of weight 2 of the form (n-1,0), (0,1),...,(n-2,n-1)
        fault_degree=0; // Constant injection
        e_weight=2;
        fault_position=SYS_T-e_weight+fault_degree; // Calculated using fault_degree and e_weight, gives where to inject in locator
        fault_value=0b1111111111111; // sets coeficcient to zero
        unsigned int minimum_nr_equations=SYS_N;
        unsigned int additional_fi=-1;
        int k1=1,k2=1;
        
        while(k1<SYS_N && (nr_found_variables < SYS_N || nr_found_equations < minimum_nr_equations || additional_fi<750)){ // needed if all variables are searched
            while(k2 <SYS_N && (nr_found_variables < SYS_N || nr_found_equations < minimum_nr_equations || additional_fi<750)){
                indizes[0]=(SYS_N+k2+k1)%SYS_N;
                indizes[1]=k2;
                gen_e_modified(e,e_weight, indizes); // construct e with e_weight at indizes positions
                
                if ( (ret_val = crypto_kem_enc_modified(ct, ss, pk,two_e)) != 0) { // Encrypt the error vector e, ss=K=Hash(1,e,C), ct=C=(C0,Hash(2,e))
                    fprintf(stderr, "crypto_kem_enc returned <%d>\n", ret_val);
                    return COMPUTATION_ERROR;
                }
                
                // Run the faulty decryption and output a faulty hash ss1
                if ( (ret_val = crypto_kem_dec_faulty(faulty_hash, ct, sk, fault_position, fault_value)) != 0) {
                    fprintf(stderr, "crypto_kem_dec returned <%d>\n", ret_val);
                    return COMPUTATION_ERROR;
                }
                nr_injections_constant++;
                
                // Write Hashes to file in format ct faulty_hash index1 index2
                for (i=0; i<crypto_kem_CIPHERTEXTBYTES; i++ )
                    fprintf(fp_hashes_lin,"%02X",ct[i]);
                fprintf(fp_hashes_lin," ");
                for (i=0; i<crypto_kem_BYTES; i++ )
                    fprintf(fp_hashes_lin,"%02X",faulty_hash[i]);
                fprintf(fp_hashes_lin," ");
                fprintf(fp_hashes_lin,"%04d %04d\n",indizes[0], indizes[1]);
                
                
                // for simulation purposes use this to verify ouput of parallelized brute force for finding e' 
                decrypt_faulty(output_faulty_e, sk +40, ct, fault_position,fault_value); //sk=(s,g(x),alpha_1,...,alpha_n), s has n-bits
                faulty_weight=get_indizes_weight1(faulty_indizes, output_faulty_e, zero_position); //  check if faulty e' has weight 1 (without zero position) 
                    
#if DEBUG_LEVEL > 2
                printf("indizes: %d, %d\n", indizes[0], indizes[1]);
                printf("faulty_indizes (%d,%d)\n\n",faulty_indizes[0],faulty_indizes[1]);
#endif
                if(faulty_weight==1 && (indizes[0]!=faulty_indizes[0]&& indizes[1]!= faulty_indizes[0]) && !(indizes[0]==faulty_indizes[1] && indizes[1]==faulty_indizes[0])){
                    // Check equations
//                  if((L[indizes[0]]^L[indizes[1]]^L[faulty_indizes[0]]^L[faulty_indizes[1]]) !=0)
//                      printf("WRONG EQ");
                            
                    fprintf(fp_sage_lin,"%d %d %d\n",indizes[0],indizes[1],faulty_indizes[0]);           
                    
                    // Sets the bit of found support elements in order to check if all are found
                    found[indizes[0]/64]=found[indizes[0]/64]|((uint64_t)1<<(indizes[0]%64));
                    found[indizes[1]/64]=found[indizes[1]/64]|((uint64_t)1<<(indizes[1]%64));
                    found[faulty_indizes[0]/64]=found[faulty_indizes[0]/64]|((uint64_t)1<<(faulty_indizes[0]%64));
            
                    nr_found_equations++;
                }

                faulty_indizes[0]=-1; faulty_indizes[1]=-1;
                
                // checks if equations contain all 6688 support eleements
                found_all=1;                
                nr_found_variables=0;
                for(i=0;i<104;i++){
                    // counts how many variables are already found
                    uint64_t result= found[i];
                    for (int j = 0; j < 64; j++){
                        nr_found_variables += (result & 1);
                        result = (result >> 1);
                    }
                    //checks if all variables are found
                    if(found[i]!=0xFFFFFFFFFFFFFFFF)
                        found_all=0;
                }
                uint64_t result= found[104];
                for (int j = 0; j < 32; j++){
                    nr_found_variables += (result & 1);
                    result = (result >> 1); 
                }
                //Checks last 32 variables
                if(found[104]!=0x00000000ffffffff)
                        found_all=0;
                
#if DEBUG_LEVEL >2
                printf("nr_variables: %d ", nr_found_variables);
                printf("nr_equations: %d\n", nr_found_equations);
#endif                
                if(nr_found_variables >= SYS_N && nr_found_equations >= minimum_nr_equations){
                    additional_fi++;
                }
                // Set ciphertext to zero
                memset(ct, 0x00, crypto_kem_CIPHERTEXTBYTES);
                // set two_e to zero
                memset(e, 0x00, SYS_N/8);
                // Set output faulty e to zero
                memset(output_faulty_e, 0x00, SYS_N/8);
                
                k2++; 
                
                
            }//end while k2<SYS_N
#if DEBUG_LEVEL >1            
            printf("nr_found_variables: %d, nr_found_equations: %d, found_all: %d, nr_loops_constant: %d\n",nr_found_variables,nr_found_equations, found_all, nr_loops_constant);
#endif            
            k2=k1+1;
            k1++;
            
            nr_loops_constant++;
            
#if DEBUG_LEVEL >2            
            // prints how many variables are already found
            printf("\n nr %d\n", nr_found_variables);
            for(i=0;i<105;i++)  
                printf("%016lx ", found[i]);
#endif
            
        }//end while k1

        nr_found_constant_equations=nr_found_equations;
        
        /************* END: CONSTANT INJECTIONS - wt(2) ***************/
        
        /************* LINEAR INJECTIONS (Quadratic Equations) ***************/
        if(linear_injection){
#if DEBUG_LEVEL >1
        printf("nr_constant_equations %d\n", nr_found_constant_equations);
#endif
        fault_degree=1;
        e_weight=2;
        fault_position=SYS_T-e_weight+fault_degree; // Calculated using fault_degree and e_weight, gives where to inject in locator
        fault_value=0b1111111111111; // sets coefficient to zero
        unsigned int nr_requested_lineq=100;
        
        int k1=1,k2=1;
                
                
        // Does linear injections until nr_requested_lineq equations are found
        while(k1<SYS_N && (nr_found_equations-nr_found_constant_equations<nr_requested_lineq)){
            while(k2 <SYS_N && (nr_found_equations-nr_found_constant_equations<nr_requested_lineq)){
                indizes[0]=(SYS_N+k2+k1)%SYS_N;
                indizes[1]=k2;
                
                gen_e_modified(e,e_weight, indizes); // construct e with e_weight at indizes positions
                
                if ( (ret_val = crypto_kem_enc_modified(ct, ss, pk,two_e)) != 0) { // Encrypt the error vector e, ss=K=Hash(1,e,C), ct=C=(C0,Hash(2,e))
                    fprintf(stderr, "crypto_kem_enc returned <%d>\n", ret_val);
                    return COMPUTATION_ERROR;
                }
                
                // Run the faulty decryption and output a faulty hash ss1
                if ( (ret_val = crypto_kem_dec_faulty(faulty_hash, ct, sk, fault_position, fault_value)) != 0) {
                    fprintf(stderr, "crypto_kem_dec returned <%d>\n", ret_val);
                    return COMPUTATION_ERROR;
                }
                nr_injections_linear++;
                
                // Write Hashes to file in format ct faulty_hash index1 index2
                for (i=0; i<crypto_kem_CIPHERTEXTBYTES; i++ )
                    fprintf(fp_hashes_quad,"%02X",ct[i]);
                fprintf(fp_hashes_quad," ");
                for (i=0; i<crypto_kem_BYTES; i++ )
                    fprintf(fp_hashes_quad,"%02X",faulty_hash[i]);
                fprintf(fp_hashes_quad," ");
                fprintf(fp_hashes_quad,"%04d %04d\n",indizes[0], indizes[1]);
                
                // for simulation purposes use this to verify ouput of parallelized brute force for finding e' 
                decrypt_faulty(output_faulty_e, sk +40, ct, fault_position,fault_value);
                faulty_weight=get_indizes_weight1(faulty_indizes, output_faulty_e, zero_position);                   
                    
                if(faulty_weight==1 && (indizes[0]!=faulty_indizes[0] || indizes[1]!= faulty_indizes[0]) && !(indizes[0]==faulty_indizes[1] && indizes[1]==faulty_indizes[0])){      
                        
                    // Check equations
//                  if((gf_mul(L[indizes[0]],L[indizes[1]])^gf_mul(L[faulty_indizes[0]],L[faulty_indizes[1]])) !=0){
//                     printf("WRONG EQ\n");
//                  }
                    fprintf(fp_sage_quad,"%d %d %d %d\n",indizes[0],indizes[1],faulty_indizes[0],faulty_indizes[0]);

                    nr_found_equations++;
                }

                faulty_indizes[0]=-1; faulty_indizes[1]=-1;
                
                
                // Set ciphertext to zero
                memset(ct, 0x00, crypto_kem_CIPHERTEXTBYTES);
                // set two_e to zero
                memset(e, 0x00, SYS_N/8);
                // Set output faulty e to zero
                memset(output_faulty_e, 0x00, SYS_N/8);
                
                k2++;
#if DEBUG_LEVEL >2               
                printf("nr_equations %d\n", nr_found_equations);
#endif                
            }//end while k2<SYS_N
#if DEBUG_LEVEL >1            
            printf("nr_found_variables: %d, nr_found_equations: %d, found_all: %d\n",nr_found_variables,nr_found_equations, found_all);
#endif            
            k2=k1+1;
            k1++;
            
        }// end k1<SYS_N
        }//end if linear injection
        
        /************* END: LINEAR INJECTIONS ***************/
        

        fprintf(fp_simulation,"nr_loops_constant: %d\n", nr_loops_constant);
        fprintf(fp_simulation,"nr_injections_constant: %d ", nr_injections_constant);
        fprintf(fp_simulation,"nr_injections_linear: %d \n", nr_injections_linear);
        fprintf(fp_simulation,"nr_constant_equations: %d\n", nr_found_constant_equations);
        fprintf(fp_simulation,"nr_linear_equations %d\n", nr_found_equations-nr_found_constant_equations);
        fprintf(fp_simulation,"All equations: %d\n", nr_found_equations);
        fprintf(fp_simulation,"\n\n");
#if DEBUG_LEVEL >0        
        printf("\nnr_loops_constant: %d\n", nr_loops_constant);
        printf("nr_injections_constant: %d, ", nr_injections_constant);
        printf("nr_injections_linear %d, ", nr_injections_linear);
        printf("\nnr_constant_equations %d,", nr_found_constant_equations);
        printf("nr_linear_equations %d\n", nr_found_equations-nr_found_constant_equations);
        printf("All equations: %d\n", nr_found_equations);
#endif        
        sum_injections_constant+=nr_injections_constant;
        sum_injections_linear+=nr_injections_linear;
        sum_equations_constant+=nr_found_constant_equations;
        sum_equations_linear+=nr_found_equations-nr_found_constant_equations;

        
        fclose(fp_sage_lin);
        fclose(fp_sage_quad);
        
        if ( (ret_val = system("./runKeyRecovery.sh")) != 0) { // Encrypt the error vector e, ss=K=Hash(1,e,C), ct=C=(C0,Hash(2,e))
            fprintf(stderr, "Failure of key recovery <%d> in cryptosystem num %d\n", ret_val, num);
            return KEY_RECOVERY_ERROR;
        } 
        
    }
    
    fprintf(fp_simulation,"\nsum_injections_constant: %d\n",sum_injections_constant);
    fprintf(fp_simulation,"sum_injections_linear: %d\n",sum_injections_linear);
    fprintf(fp_simulation,"sum_equations_constant: %d\n",sum_equations_constant);
    fprintf(fp_simulation,"sum_equations_linear: %d\n",sum_equations_linear);
    
    fprintf(fp_simulation,"\nmean injections_constant: %d\n",sum_injections_constant/KATNUM);
    fprintf(fp_simulation,"mean injections_linear: %d\n",sum_injections_linear/KATNUM);
    fprintf(fp_simulation,"mean equations_constant: %d\n",sum_equations_constant/KATNUM);
    fprintf(fp_simulation,"mean equations_linear: %d\n",sum_equations_linear/KATNUM);
    
    
    fclose(fp_simulation);
    fclose(fp_rsp);
    return 0;
}

void
fprintBstr(FILE *fp, char *S, unsigned char *A, unsigned long long L)
{
	unsigned long long i;

	fprintf(fp, "%s", S);

	for ( i=0; i<L; i++ )
		fprintf(fp, "%02X", A[i]);

	if ( L == 0 )
		fprintf(fp, "00");

	fprintf(fp, "\n");
}


/*
 * Determines the non-zero positions of the fauly error vector and removes the position of the zero element (if exists)
 * Input: faulty error vector, position of zero element
 * Output: positions/indizes of faulty error vector if Hamming weight is 1
 * Returns: Hamming weight of error vector with removed position of zero element
 */
int get_indizes_weight1(uint16_t* faulty_indizes, unsigned char * output_faulty_e, uint16_t zero_position){
    int k, found=0;
    for (k = 0;k < SYS_N;++k){
        if (output_faulty_e[k/8] & (1 << (k&7))){
            faulty_indizes[found]=k;
            found++;  
        }
    }
    // remove zero_position
    if(found==2 && faulty_indizes[0]==zero_position){
        faulty_indizes[0]=faulty_indizes[1];
        faulty_indizes[1]=-1;
        found--;
    }
    if(found==2 && faulty_indizes[1]==zero_position){
        faulty_indizes[1]=-1;
        found--;
    }
    //remove all equations that contain the zero_position, because it messes up the polynomial system
    if(faulty_indizes[0]==zero_position){
        found --;
    }
    
    return found;
}

/*
 * Determines the position of the zero element in the secret key (if exists)
 * which can be found from error vector if plaintext has smaller weight than SYS_T
 * Input: hash output, ciphertext
 * Output: zero_element (if exists)
 * Return: 0 if zero_element could be found, 1 if secret key does not contain the zero element
 */
int search_zero_element(uint16_t*faulty_indizes,unsigned char* faulty_hash, unsigned char* ct){
    uint16_t ret=1;
    int i=0;
    unsigned char one_ec[ 1 + SYS_N/8 + (SYND_BYTES + 32) ] = {1};
    unsigned char computed_hash[32];
    unsigned char try_e[SYS_N/8];
    uint16_t ind[1];
    unsigned char exists_zero_position=1;
  
#if DEBUG_LEVEL >2    
    printf("Faulty hash: \n");
    for(int j=0;j<32;j++)
        printf("%.2X",faulty_hash[j]);
    printf("\n");
#endif    
    
    memset(try_e, 0x00, SYS_N/8);
    memcpy(one_ec + 1, try_e, SYS_N/8);  // Concatenates (1,e')=(1,0)
    memcpy(one_ec + 1 + SYS_N/8, ct, SYND_BYTES + 32); // Concatenates (1,e',ct)
    crypto_hash_32b(computed_hash, one_ec, sizeof(one_ec)); // = Hash(1,e',C)
    
#if DEBUG_LEVEL >2     
    printf("computed hash: \n");
    for(int j=0;j<32;j++)
        printf("%.2X",computed_hash[j]);
    printf("\n");
#endif
    
    if(!memcmp(faulty_hash,computed_hash,crypto_kem_BYTES)){ // Compare the computed_hash with the faulty_hash
        exists_zero_position=0;
        ret=0;
    }

    while(exists_zero_position && ret!=0 && i<SYS_N){
        ind[0]=i;
        gen_e_modified(try_e,1,ind); // Generates different e' with wt(e')=1
                
        memcpy(one_ec + 1, try_e, SYS_N/8);  // Concatenates (1,e')
        memcpy(one_ec + 1 + SYS_N/8, ct, SYND_BYTES + 32); // Concatenates (1,e',ct)
        crypto_hash_32b(computed_hash, one_ec, sizeof(one_ec)); // = Hash(1,e',C)

        if(!memcmp(faulty_hash,computed_hash,crypto_kem_BYTES)){ // Compare the computed_hash with the faulty_hash
            faulty_indizes[0]=ind[0];
            ret=1;
        }

        i++;
    }
    return ret; // outputs 0 if no faulty_indizes could be found, otherwise 1
}


static inline uint32_t same_mask(uint16_t x, uint16_t y)
{
	uint32_t mask;

	mask = x ^ y;
	mask -= 1;
	mask >>= 31;
	mask = -mask;

	return mask;
}

/* Generates an error vector e of Hamming weight e_weight with specified non-zero positions
 * Input: desired Hamming weight e_weight, non-zero positions ind
 * Output: e, an error vector of weight e_weight 
 */
static void gen_e_modified(unsigned char *e, unsigned int e_weight, uint16_t *ind)
{
	int i, j;
	unsigned char mask;	
	unsigned char val[ e_weight];

	for (j = 0; j < e_weight; j++)
		val[j] = 1 << (ind[j] & 7);

	for (i = 0; i < SYS_N/8; i++) 
	{
		e[i] = 0;

		for (j = 0; j < e_weight; j++)
		{
			mask = same_mask(i, (ind[j] >> 3));

			e[i] |= val[j] & mask;
		}
	}
}
