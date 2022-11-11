# Key-Recovery Fault Attack Simulation on Classic McEliece KEM

This repository contains all the code to simulate our fault attack on Classic McEliece as described in our [preprint](https://eprint.iacr.org/2022/1529).


It contains C-code to simulate the fault injections on the [Classic McEliece](https://classic.mceliece.org/nist.html) Round 3 C-code implementations and to brute-force the input of the hashed session keys. A SageMath script can use these inputs to build a system of polynomial equations, find suitable solutions, and calculate an alternative secret key, using a C-program for validation.

    .
    ├── 01_AttackSimulation                 # Classic McEliece Round 3 code with attack simulation
    │   └── ...                             # Simulations for ELPb and ELPz for different parameter sets
    ├── 02_DehashingBruteForce              # BruteForce to find input of Hash output
    │   └── ...                             # Code for different parameter sets
    ├── 03_ComputeAlternativeSecretKey
    │   ├── comparePK                       # C program to compare public keys
    │   └── comp_alt_key.sage               # SageMath program for calculating an alternative secret key
    └── README.md

### Prerequisite
The implementations run on Linux systems with x86-64 processors integrating the Advanced Vector Extension (AVX2).  
It needs `gcc`, `make` and `SageMath`. (We used `gcc version 12.2.0-1` and `SageMath version 9.7`.)  



### 1. Attack Simulation
The programs in `01_AttackSimulation` contain the main simulation code for different parameters and fault attack szenarios ELPb and ELPz (see the paper for more details). The programs simulate the fault injections, output the faulty session keys, and the linear and quadratic equations and call the SageMath script that calculates an alternative secret key.  

The simulation code in `01_AttackSimulation` uses the AVX implementation of Classic McEliece for faster running time. Therefore, the code runs only on x86-64 processors with integrated AVX2 Extension.


The simulation is started by executing the following commands (where XXXX stands for the different parameter sets):  

    make                        # Execute in folder attack_mcelieceXXXXXX_avx, generates executable in parent folder
    cd ..
    ./runSimulationXXXX.sh      # Execute in parent folder parametersXXXX  


The generated output files are `input_hashes_lin.txt` and `input_hashes_quad.txt` containing the hashes of the faulty outputs. These files are needed for the brute-force program in `02_DehashingBruteForce`. Additionally, there are the files `input_sage_lin.sage` and `input_sage_quad.sage` generated from the main simulation that contain the linear and quadratic equations such that the solving of the system of polynomial equations and the computation of an alternative secret key can directly be called from the simulation program. All output files will be stored in the folder `parametersXXXX`.  

The simulation will automatically run the `runKeyRecovery.sh` script that calls `comp_alt_key.sage`. (Make sure that the symbolic links of `comp_alt_key.sage` and `comparePK_XXXX` in folders `parametersXXXX ` point to the right files.)

The alternative secret key can then be found in the file `output_sage_a.txt` (support) and `output_sage_g.txt` (Goppa polynomial). 


### 2. Brute-force on Hashed Output
The programms in `02_DehashingBruteForce` prove the feasibility of our attack and show that the output session keys can be de-hashed in acceptable time. (For the ELPb case the dehashing of a full attack takes up to a few days on an AMD EPYC 7543P 3.3 GHz with 64 threads. For the ELPz case it takes less than a second.)  
Sample input files are provided, but if you want to run the code for a specific instance then copy the files `input_hashes_lin.txt` and `input_hashes_quad.txt` from the attack simulation output to the corresponding dehashing folder and adjust or remove the position of the zero element in `run.sh`.  
In the case of ELPb manually delete the first line (that contains the position of the zero element) in the output file `input_sage_quad.sage` such that the quadratic equations can be properly read by the SageMath script.

### 3. Key-Recovery Calculation
The SageMath script `comp_alt_key.sage` is called by the main simulation code (in `01_AttackSimulation`) and does not need to be run manually. It reads linear and quadratic equations from files `input_sage_lin.sage` and `input_sage_quad.sage`, and outputs the support of the Goppa code in file `output_sage_a.txt` and the Goppa polynomial in file `output_sage_g.txt`.  
The `comparePK` program compares two public keys and is automatically called from the SageMath script.  
Binarys for `comparePK` can be generated using `make comparePK_XXXX` where XXXX stands for the different parameter sets.
