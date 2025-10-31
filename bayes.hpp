#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string>


void Bayes(int wsqA,int wsqT,int wsqC, int wsqG, int crqA, int crqT, int crqC, int crqG, char refbase, unsigned int pos, char* chrom,
    int w_A, int w_T, int w_C, int w_G, int c_A, int c_T, int c_C, int c_G, std::string& genotypemaybe, double &qual);
