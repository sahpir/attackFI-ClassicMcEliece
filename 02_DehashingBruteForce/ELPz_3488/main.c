#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include "keccak/SimpleFIPS202.h"
#include <time.h>

#define crypto_hash_32b(out,in,inlen) \
  SHAKE256(out,32,in,inlen)

#define NR_THREADS 64

#define FILE_OPEN_ERROR -1
#define THREAD_ERROR -2
#define INPUT_ERROR -3
#define COMPUTATION_ERROR -4
// #define crypto_kem_CIPHERTEXTBYTES
// #define SYND_BYTES
// #define SYS_N


struct Element{
    unsigned char ct[crypto_kem_CIPHERTEXTBYTES]; // ct=(c0,c1)=(eH^T,Hash(2,e))
    unsigned char faulty_hash[32]; // output of FI on ELP
    uint16_t indizes[2]; // indizes of e used to generate ct and inject FI
};

struct Arguments{
    char out_eq_type; // If out_eq_type=1 output for linear equations, if out_eq_type=2 output for quadratic equations
    char *filename_input;
    char *filename_output_thread;
    unsigned int thread;
    uint32_t beginReadingAtPosition;
    uint32_t read_nr_lines;
    uint16_t zero_position;
    unsigned char exists_zero_position;
};

static inline uint32_t same_mask(uint16_t x, uint16_t y)
{
	uint32_t mask;

	mask = x ^ y;
	mask -= 1;
	mask >>= 31;
	mask = -mask;

	return mask;
}

/* output: e, an error vector of weight e_weight */
static void gen_e_modified(unsigned char *e, unsigned char e_weight, uint16_t *ind)
{
	int i, j;//, eq;
	unsigned char mask;	
	unsigned char val[ e_weight];
    
    if(e_weight<=2){ // checks if e_weight is not bigger than ind size of 2
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
}

void *findErrorVector(void *_arguments){
    struct Arguments *parameters=(struct Arguments *) _arguments;
    FILE *fp_input, *fp_output;
    //char filename_output_thread[strlen(parameters->filename_output_thread)];
    uint16_t ret=0;
    uint16_t i=0;
    uint32_t lines_processed=0;
    unsigned char one_ec[ 1 + SYS_N/8 + (SYND_BYTES + 32) ] = {1}; // used to concatenate and calculate Hash
    unsigned char try_e[SYS_N/8]; // ind are indizes of try_e
    unsigned char computed_hash[32];
    uint16_t ind[2]; // try indizes of e'
    uint16_t faulty_indizes[3];
    struct Element element;
    
    
    
    // Assumption that input file has not changed on disk during start of program
    // Read from File the faulty_hash and indizes from those lines in file given by parent process (split input file for all threads equally)
//     printf("filename_input: %s\n", parameters->filename_input);
    fp_input = fopen(parameters->filename_input, "r");
    if (!fp_input)
        return  (void*)FILE_OPEN_ERROR;
    
    fseek(fp_input, parameters->beginReadingAtPosition, SEEK_SET); // Sets the pointer to read from beginReadingAtLine
    
     // Open outputfile with write access
    fp_output = fopen(parameters->filename_output_thread, "w");
    if (!fp_output)
            pthread_exit((void*)FILE_OPEN_ERROR);
    
    while (lines_processed < parameters->read_nr_lines){
//         printf("T%d: %ld\n",parameters->thread, ftell(fp_input));
        for(i=0;i<crypto_kem_CIPHERTEXTBYTES;i++){ // Read crypto_kem_CIPHERTEXTBYTES 
            fscanf(fp_input, "%2hhx",&element.ct[i]);
//             printf("%.2X",element.ct[i]);
        }
//         printf("\n");
        for(i=0;i<32;i++){ // Read 32 Bytes for Hash
            fscanf(fp_input, "%2hhx",&element.faulty_hash[i]);
//             printf("%.2X",element.faulty_hash[i]);
        }
//         printf("\n");
        fscanf(fp_input, "%4hd %4hd",&element.indizes[0], &element.indizes[1]);
//         printf("%d %d",element.indizes[0], element.indizes[1]);
//         printf("\n");
        
        
        
        
        // Do Calculations that needs ct, faulty_hash, zero_position, exists_zero_position
        // Search for wt(e')=1 +exists_zero_position
        
        if(parameters->exists_zero_position){ // zero element exists in support
            ind[1]=parameters->zero_position;
        }
        
        i=0;ret=0;
        while (ret==0 && i<SYS_N){
            ind[0]=i;
            gen_e_modified(try_e,1+parameters->exists_zero_position,ind); // Generates different e' with wt(e')=1(+1 for zero position)
            memcpy(one_ec + 1, try_e, SYS_N/8);  // Concatenates (1,e')
            memcpy(one_ec + 1 + SYS_N/8, &element.ct, SYND_BYTES + 32); // Concatenates (1,e',ct)
            crypto_hash_32b(computed_hash, one_ec, sizeof(one_ec)); // = Hash(1,e',C)
                
//             printf("computed_hash <%d>: \n",i);
//             for(a=0;a<32;a++)
//                 printf("%.2X",computed_hash[a]);
//             printf("\n");
            
            if(!memcmp(element.faulty_hash,computed_hash,32)){ // Compare the computed_hash with the faulty_hash
                faulty_indizes[0]=ind[0];
                ret=-1;
            }
            i++;
        }
        
        
//         printf("indizes: %d, %d\n", element.indizes[0], element.indizes[1]);
//         printf("faulty_indizes (%d,%d)\n\n",faulty_indizes[0],faulty_indizes[1]);
           
        // Write to outputfile depending on linear equations or quadratic equations specified in out_eq_type
        if(faulty_indizes[0]!=(uint16_t)-1 && faulty_indizes[0]!=parameters->zero_position && element.indizes[0]!=faulty_indizes[0] && element.indizes[1]!=faulty_indizes[0]){
            if(parameters->out_eq_type==1){
//             printf("Written nr of characters: %d\n",fprintf(fp_output,"%d %d %d %d\n",element.indizes[0],element.indizes[1],faulty_indizes[0],faulty_indizes[1]));
                fprintf(fp_output,"%d %d %d\n",element.indizes[0],element.indizes[1],faulty_indizes[0]);
            }else if (parameters->out_eq_type==2){
                fprintf(fp_output,"%d %d %d %d\n",element.indizes[0],element.indizes[1],faulty_indizes[0],faulty_indizes[0]);
            }
        }
        
        faulty_indizes[0]=-1;
    
        
        lines_processed++;
    }
 
    fclose(fp_input);
    fclose(fp_output);
    //free(parameters);
    pthread_exit(0);
    
}

    
    // Preliminary : searched_e_weight=1 for all computations, zero-position is the same for all computations
    // Einlesen von input File (consists of ciphertext ct, faulty_hash, indizes) 
    // Calculate Hashes of Hash(1,e',ct) of different e' with wt(e')=1 and compare with faulty_hash
    // Output: polynomial system of equations (e.g. 3487 0 862 1291 per line)
    // Compile with -pthread
    // argv reihenfolge: program name, "lin" or "quad", filename_input, filename_output, (zero_position)

int main(int argc, char **argv){

    FILE *fp_input, *fp_output, *fp_threads[NR_THREADS]; 
    uint16_t zero_position=-1; // default is no zero element in secret key
    unsigned char exists_zero_position=0; // default is no zero element in secret key
    pthread_t threads[NR_THREADS];
    struct Arguments argumentsForThread[NR_THREADS];
    unsigned int error;
    int i;
    char c;
    char out_eq_type=0;
    
    if(argc<4)
        return FILE_OPEN_ERROR;
    
    char filename_output[strlen(argv[3])+3]; //+3 for storing the filename of the NR_THREADS
    
    if(!strcmp("lin",argv[1])){
        out_eq_type=1;
    }else if(!strcmp("quad",argv[1])){
        out_eq_type=2;
    }else
        return INPUT_ERROR;

    
    // Input file has the format ct faulty_hash index1 index2
    // Open File to split equally on the number of threads
    fp_input = fopen(argv[2], "r");
    if (!fp_input)
        return FILE_OPEN_ERROR;
    
    
    
    if(argc==5){
        zero_position=(uint16_t)(atoi(argv[4]));
        exists_zero_position=1;
//         printf("zero_position %d\n", zero_position);
    }
    
    
    
    /** Split the file readings equally to all threads **/
     
    // Extract characters from file and store in character c
    uint32_t total_lines=0, total_lines2=0;
    //for (char c = getc(fp_input); c != EOF; c = getc(fp_input))
    while((c=fgetc(fp_input))!=EOF) {
        if (c == '\n') // Increment if this character is newline
            total_lines++;
    }
    printf("The file has %d lines\n", total_lines);
    
    fseek(fp_input,0,SEEK_END); // Find out how many lines/elements has the file
    uint32_t len = ftell(fp_input);
    printf("Total size of file.txt = %d bytes\n", len);
    total_lines2=len/(crypto_kem_CIPHERTEXTBYTES*2+32*2+4*2+4); //format ct_bytes*2+" "+hash_Bytes*2+" "+index1_4bytes" "+index2_4bytes+"\n" characters in each line
    printf("Total number of lines = %d\n\n",total_lines2);
    
    if(total_lines != total_lines2)
        return COMPUTATION_ERROR; // Content of input file does not fit the required format
    
    
    fclose(fp_input);
    
    /****** Start the Threads ******/
    // Every thread needs to read at least the total number of lines read_nr_lines= total_lines/NR_THREADS
    // The remaining total_lines%NR_THREADS need to be taken from the first few threads
    // First thread starts reading from the bottom of the file
    
    clock_t tic = clock();
    
    uint32_t read_nr_lines=0;
    uint32_t old_position=len;
    // Creates the Threads with different arguments
    for (i=0; i<NR_THREADS; ++i)
    {
        read_nr_lines=total_lines/NR_THREADS;
        if(i<(total_lines%NR_THREADS)) // First total_lines MOD NR_THREADS threads need to take over one more line
            read_nr_lines++;
        
        argumentsForThread[i].out_eq_type=out_eq_type;
        argumentsForThread[i].filename_input=argv[2];
        argumentsForThread[i].filename_output_thread=malloc(sizeof((argv[3]))+3*sizeof(char));
        argumentsForThread[i].zero_position=zero_position;
        argumentsForThread[i].exists_zero_position=exists_zero_position;
        
        argumentsForThread[i].thread=i;
        printf("Thread: %d\nread_nr_lines: %d\n",argumentsForThread[i].thread,read_nr_lines);
        sprintf(argumentsForThread[i].filename_output_thread,"%s%.3d", argv[3],i);
        argumentsForThread[i].beginReadingAtPosition=old_position-(read_nr_lines*(crypto_kem_CIPHERTEXTBYTES*2+32*2+4*2+4)); // see format above
        old_position=argumentsForThread[i].beginReadingAtPosition;
        argumentsForThread[i].read_nr_lines=read_nr_lines;
        printf("filename_output_thread: %s\nBeginReadingAtPosition: %u\n\n", argumentsForThread[i].filename_output_thread,argumentsForThread[i].beginReadingAtPosition );
        error = pthread_create(&threads[i], NULL, findErrorVector, (void *)&argumentsForThread[i]);
        
        
        if(error) // Check if needed!
            return THREAD_ERROR;
        read_nr_lines=0;
    }
    
    printf("All Threads started!\n\n");

    fflush(NULL);
    // Waits for threads to finish
    int return_value=0;
    for (i=0; i<NR_THREADS; ++i) {
        error = pthread_join(threads[i], (void*)&return_value);
        if(error) // Check if needed!
            return THREAD_ERROR;
        if(return_value!=0)
            return return_value;
    }
    
    clock_t toc = clock();
    printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
  
    /***** FILE CONCATENATION *************/
    // Each thread creates its own output file
    // Combine the NR_THREADS output files into one file
    
    fp_output = fopen(argv[3], "w");
    if (!fp_output)
        return FILE_OPEN_ERROR;
    if(out_eq_type==1 && exists_zero_position)
        fprintf(fp_output,"%d\n", zero_position);
    
    for(i=NR_THREADS-1;i>=0;i--){
        sprintf(filename_output,"%s%.3d",argv[3],i);
        printf("%s\n",filename_output);
        fp_threads[i] = fopen(filename_output, "r");
        if (!fp_threads[i])
            return FILE_OPEN_ERROR;
        // Copy content of file character by character
        while ((c = fgetc(fp_threads[i])) != EOF)
            fputc(c, fp_output);
        
        fclose(fp_threads[i]);
    }
    
    fclose(fp_output);
    
    return 0;
}

