#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bayes.hpp"
#include "fisher.hpp"

extern int nLayerMax;
extern int vQualMin;
extern int nLayerMin;
extern int nLayerMax;
extern float vSnpRate;
extern float minhomfreq;
extern float pvalue_cutoff;
int minhetfreq = vSnpRate;
extern int minread2;

//void printSNP(char* chrom, int pos, char refbase, char* var, double pvalue, int w_A, int w_T, int w_C, int w_G, int c_A, int c_T, int c_C, int c_G, 
//    int genoqual, char* filter){
void printSNP(char* chrom, int pos, char refbase, double pvalue, int w_A, int w_T, int w_C, int w_G, int c_A, int c_T, int c_C, int c_G, 
    char* filter){
    printf("%s\t%d\t%c\t%s\t%f\t%d,%d,%d,%d\t%d,%d,%d,%d\n", chrom, pos, refbase, filter, pvalue,w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G);
}

int genotype(int wsqA,int wsqT,int wsqC, int wsqG, int crqA, int crqT, int crqC, int crqG, char refbase, unsigned int pos, char* chrom,
    int w_A, int w_T, int w_C, int w_G, int c_A, int c_T, int c_C, int c_G, std::string genotypemaybe, double genoqual)
{
    //print "lines[0]\tlines[1]\trefbase\tgenotypemaybe\n";

    unsigned int totaldepth=w_A+w_T+w_C+w_G+c_A+c_T+c_C+c_G;
    if(totaldepth > nLayerMax){
        return 0;
    }
    int Ctotaldepth=c_C+c_A+c_T+c_G;
    int Wtotaldepth=w_A+w_T+w_C+w_G;
    if(genotypemaybe == "AA"){//genotypeis AA
        if(refbase == 'A'){
            return 0;   
        }
        if(refbase == 'T'){//T>AA
            double qvalue=(wsqA>crqA)?wsqA:crqA;
            int expectvar = (int) pow(0.1, (qvalue/10))/3;
            int var=w_A+c_A;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,var, totaldepth);
            
            if(totaldepth >= nLayerMin  && qvalue >= vQualMin && var >=minread2 ){
                //sprintf("%.2f", f)
                double T2A= (double) var/totaldepth;

                if(T2A>=minhomfreq){
                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS");
                }else{
                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                }
                
            }else{
                double T2A=(double) var/totaldepth;
                printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");              
            }                                   
        }
        if(refbase == 'C'){//C>AA
            double qvalue=(wsqA>crqA)?wsqA:crqA;
            //int depth=w_C+c_C+w_A+c_A;
            int var=w_A+c_A;   
            int expectvar = (double) pow(0.1, (qvalue/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,var, totaldepth);
            if(totaldepth >= nLayerMin  && qvalue >= vQualMin && var >=minread2 ){
                                double C2A=(double)var/totaldepth;
                                if(C2A>=minhomfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS");
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                                }   
    
                        }else{
                double C2A= (double) var/totaldepth;
                printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                        }                   
        }

        if(refbase == 'G'){//G>AA
            double qvalue=wsqA;
                        int depth=w_G+w_A;
                        int var=w_A;
            double G2A;
            if(depth>0){
                G2A=(double) var/depth;
            }else{
                G2A=0;
            }
            int expectvar = (double) pow(0.1, (qvalue/10))/3;
            double pvalue = fishers_exact((int)expectvar*depth,depth,var, depth);

                        if(depth >= nLayerMin  && qvalue >= vQualMin && var >=minread2 ){
                                if(G2A>=minhomfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS"); 
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                                }

                        }else{
                            printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                        }
        }
    }

    if(genotypemaybe == "AT"){ //genotype is AT
        if(refbase == 'A'){//A>AT
                double qvalue=(wsqT>crqT)?wsqT:crqT;;
            //int depth=w_A+w_T+c_A+c_T;
            int var=w_T+c_T;   
            int expectvar = (double) pow(0.1, (qvalue/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,var, totaldepth);
            if(totaldepth >= nLayerMin  && qvalue >= vQualMin && var >=minread2 ){
                                double A2T=(double) var/totaldepth;
                                if(A2T>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS"); 
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                                }
                        }else{
                double A2T=(double) var/totaldepth;
                printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                        } 
        }   

                if(refbase == 'T'){//T>AT          
                    double qvalue=(wsqA>crqA)?wsqA:crqA;
                        //int depth=w_A+w_T+c_A+c_T;
                        int var=c_A+w_A;
                        int expectvar = (double) pow(0.1, (qvalue/10))/3;
                    double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,var, totaldepth);

                        if(totaldepth >= nLayerMin  && qvalue >= vQualMin && var >=minread2 ){
                                double T2A=(double) var/totaldepth;
                                if(T2A>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS"); 
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                                }

                        }else{
                double T2A=(double) var/totaldepth;
                                printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                        }
                }   
    if(refbase == 'C'){//C>AT
            double qvalueA=(wsqA>crqA)?wsqA:crqA;
            double qvalueT=crqT;
            int varA=w_A+c_A;
            int varT=c_T;
            int depth=w_A+c_A+c_T;
            int expectvar = (double) pow(0.1, ((qvalueA+qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varA+varT, totaldepth);

            if(totaldepth >= nLayerMin  && qvalueA >= vQualMin && qvalueT>=vQualMin && varA >=minread2 && varT>=minread2 ){
                double C2A=(double) varA/totaldepth;
                double C2T=(double) varT/Ctotaldepth;
                                if(C2A>=minhetfreq && C2T>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS"); 
                                }else{
                                   printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                                }   

                        }else{
                                double C2A=(double) varT/totaldepth;
                                double C2T=(double) varT/totaldepth;
                                printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                        }       
            
                }   
                if(refbase == 'G'){//G>AT  
            double qvalueA=wsqA;
                        double qvalueT=(wsqT>crqT)?wsqT:crqT;
                        int varA=w_A;
                        int varT=c_T+w_T;
                        int depth=w_A+c_T+w_T;
            int expectvar = (double) pow(0.1, ((qvalueA+qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varA+varT, totaldepth);

                        if(totaldepth >= nLayerMin  && qvalueA >= vQualMin && qvalueT>=vQualMin && varA >=minread2 && varT>=minread2 ){
                                double G2A=(double) varA/Wtotaldepth;
                                double G2T=(double) varT/totaldepth;
                                if(G2A>=minhetfreq && G2T>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS"); 
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                                }

                        }else{
                                double G2A=(double) varA/totaldepth;
                                double G2T=(double) varT/totaldepth;
                                printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                        }
                
                }                   
    }
//AC 
    if(genotypemaybe == "AC"){ //genotype is AC
                if(refbase == 'A'){ //A>AC
                        double qvalueC=(wsqC>crqC)?wsqC:crqC;
                        int varC=w_C+c_C+w_T;
                        //int depth=w_A+c_C+w_C;
            double A2C=(double) varC/totaldepth;
            int expectvar = (double) pow(0.1, (qvalueC/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varC, totaldepth);

                        if(totaldepth >= nLayerMin  && qvalueC >= vQualMin  && varC >=minread2 ){
                                if(A2C>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS"); 
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                                }   

                        }else{
                            printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                        }               
                }
                if(refbase == 'T'){//T>AC ////////////////////////modify totaldepth
            double qvalueA=(wsqA>crqA)?wsqA:crqA;
                        double qvalueC=(wsqC>crqC)?wsqC:crqC;
                        int varA=w_A+c_A;
                        int varC=c_C+w_C;
                        //int depth=w_A+c_C+w_C+c_A+c_T;
            double T2A=(double) varA/totaldepth;
                        double T2C=(double) varC/totaldepth;
            int expectvar = (double) pow(0.1, ((qvalueC+qvalueA)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varC+varA, totaldepth);

                        if(totaldepth >= nLayerMin  && qvalueA >= vQualMin && qvalueC>=vQualMin && varA >=minread2 && varC>=minread2 ){
                                if(T2A>=minhetfreq && T2C>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS"); 
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                                }   

                        }else{
                            printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                        }       
        
                }
                if(refbase == 'C'){//C>AC
            double qvalueA=(wsqA>crqA)?wsqA:crqA;
                        int varA=w_A+c_A;
                        int depth=w_A+c_C+w_C;
                        double C2A=(double) varA/totaldepth;
                        int expectvar = (double) pow(0.1, ((qvalueA)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varA, totaldepth);

                        if(totaldepth >= nLayerMin  && qvalueA >= vQualMin  && varA >=minread2 ){
                                if(C2A>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS"); 
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                                }

                        }else{
                            printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                        }           
                }
                if(refbase == 'G'){//G>AC
            double qvalueA=wsqA;
                        double qvalueC=(wsqC>crqC)?wsqC:crqC;
                        int varA=w_A;
                        int varC=c_C+w_C;
                        int depth=w_A+c_C+w_C+c_A;
                        double G2A=(double) varA/Wtotaldepth;
                        double G2C=(double) varC/totaldepth;
                        int expectvar = (double) pow(0.1, ((qvalueC+qvalueA)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varC+varA, totaldepth);

                        if(totaldepth >= nLayerMin  && qvalueA >= vQualMin && qvalueC>=vQualMin && varA >=minread2 && varC>=minread2 ){
                                if(G2A>=minhetfreq && G2C>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS"); 
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                                }

                        }else{
                            printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                        }
                }
        }   
    
//AG 
    if(genotypemaybe == "AG"){
                if(refbase == 'A'){//A>AG
            double qvalueG=(wsqG>crqG)?wsqG:crqG;
                        int varG=w_G+c_G;
                        int depth=w_A+c_G+w_G+w_T+c_T+w_C+c_C;
                        double A2G=(double) varG/depth;
                        int expectvar = (double) pow(0.1, ((qvalueG)/10))/3;
            double pvalue = fishers_exact((int)expectvar*depth,depth,varG, depth);
                        if(depth >= nLayerMin  && qvalueG >= vQualMin  && varG >=minread2 ){
                                if(A2G>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS"); 
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                                }   

                        }else{
                            printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                        }               
                }
                if(refbase == 'T'){//T>AG
            double qvalueA=wsqA;
                        double qvalueG=(wsqG>crqG)?wsqG:crqG;
                        int varA=w_A;
                        int varG=c_G+w_G;
                        int depth=w_A+c_G+w_G+c_A+c_T+w_T;
                        double T2A=(double) varA/Wtotaldepth;
                        double T2G=(double) varG/totaldepth;
                        int expectvar = (double) pow(0.1, ((qvalueG+qvalueA)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varG+varA, totaldepth);
                        if(totaldepth >= nLayerMin  && qvalueA >= vQualMin && qvalueG>=vQualMin && varA >=minread2 && varG>=minread2 ){
                                if(T2A>=minhetfreq && T2G>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS"); 
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                                }

                        }else{
                            printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                        }
    
                }
                if(refbase == 'C'){//C>AG
            double qvalueA=wsqA;
                        double qvalueG=(wsqG>crqG)?wsqG:crqG;
                        int varA=w_A;
                        int varG=c_G+w_G;
                        int depth=w_A+c_C+w_C+c_A+c_T+w_T;
                        double C2A=(double) varA/Wtotaldepth;
                        double C2G=(double) varG/totaldepth;
                        int expectvar = (double) pow(0.1, ((qvalueG+qvalueA)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varG+varA, totaldepth);
                        if(depth >= nLayerMin  && qvalueA >= vQualMin && qvalueG>=vQualMin && varA >=minread2 && varG>=minread2 ){
                                if(C2A>=minhetfreq && C2G>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS"); 
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                                }

                        }else{
                            printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                        }
    
                }
                if(refbase == 'G'){//G>AG
            double qvalueA=wsqA;
                        int varA=w_A;
                        int depth=w_G+w_A+w_T+w_C;
                        //int depth=c_G+w_G+w_A+w_T+c_T+w_C+c_C;
                        double G2A=(double) varA/depth;
                        int expectvar = (double) pow(0.1, ((qvalueA)/10))/3;
            double pvalue = fishers_exact((int)expectvar*depth,depth,varA, depth);
                        if(depth >= nLayerMin  && qvalueA >= vQualMin  && varA >=minread2 ){
                                if(G2A>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS"); 
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                                }

                        }else{
                            printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                        }
    
                }
        }   
//TT
    if(genotypemaybe == "TT"){
                if(refbase == 'A'){//A>TT
            double qvalueT=crqT;
                        int depth=w_A+c_A+c_T;
                        int varT=c_T+w_T;
        int expectvar = (double) pow(0.1, ((qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varT, totaldepth);

                        if(totaldepth >= nLayerMin  && qvalueT >= vQualMin && varT >=minread2 ){
                                double A2T=(double) varT/totaldepth;
                                //print "\nA2T\t", A2T, "\t", varT, "\t", totaldepth;
                                if(A2T>=minhomfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS"); 
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                                }

                        }else{
                                double A2T=(double)varT/totaldepth;
                                printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                        }
                }
                if(refbase == 'T'){
            return 0; 
                }
                if(refbase == 'C'){//C>TT
            double qvalueT=crqT;
                        int depth=c_T+c_G;
                        int varT=c_T;
            double C2T;
            if(depth>0){
                C2T=(double) varT/Ctotaldepth;
            }else{
                C2T=0;
                
            }
            int expectvar = (double) pow(0.1, ((qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*Ctotaldepth,Ctotaldepth,varT, Ctotaldepth);
            
                        if(depth >= nLayerMin  && qvalueT >= vQualMin && varT >=minread2 ){
                                if(C2T>=minhomfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS"); 
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                                }

                        }else{
                            printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                        }
    
                }
                if(refbase == 'G'){//G>TT
            double qvalueT=crqT;
                        int depth=w_G+c_G+c_T;
                        int varT=c_T+w_T;
            int expectvar = (double) pow(0.1, ((qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varT, totaldepth);
                        if(depth >= nLayerMin  && qvalueT >= vQualMin && varT >=minread2 ){
                                double G2T=(double) varT/totaldepth;
                                if(G2T>=minhomfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS"); 
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                                }

                        }else{
                                double G2T=(double) varT/totaldepth;
                                printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low"); 
                        }
    
                }
        }  

//CC     
        if(genotypemaybe == "CC"){
                if(refbase == 'A'){//A>CC
            double qvalueC=(crqC>wsqC)?crqC:wsqC;
                        int depth=w_A+c_A+c_C+w_C+w_G+w_G+c_T;
                        int varC=w_C+c_C;
            int expectvar = (double) pow(0.1, ((qvalueC)/10))/3;
            double pvalue = fishers_exact((int)expectvar*depth,depth,varC, depth);

                        if(depth >= nLayerMin  && qvalueC >= vQualMin && varC >=minread2 ){
                                double A2C=(double) varC/depth;
                                if(A2C>=minhomfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS");
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                                }

                        }else{
                                double A2C=(double) varC/depth;
                                printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                        }   
                }   
                if(refbase == 'T'){//T>CC
            double qvalueC=(crqC>wsqC)?crqC:wsqC;
                        int depth=c_T+c_C+w_C+w_A+c_A+w_G+c_G;
                        int varC=w_C+c_C;
            int expectvar = (double) pow(0.1, ((qvalueC)/10))/3;
            double pvalue = fishers_exact((int)expectvar*depth,depth,varC, depth);

                        if(depth >= nLayerMin  && qvalueC >= vQualMin && varC >=minread2 ){
                                double T2C=(double) varC/depth;
                                double T2Cc=(double) c_C/Ctotaldepth;
                                if(T2C>=minhomfreq || T2Cc>=minhomfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS");
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                                }

                        }else{
                                double T2C=(double) varC/depth;
                                printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                        }

                }   
                if(refbase == 'C'){//C>CC  
            return 0;
                }   
                if(refbase == 'G'){//G>CC
            double qvalueC=(crqC>wsqC)?crqC:wsqC;
                        int depth=w_G+c_G+c_C+w_C+w_A+c_A+c_T+w_T;
                        int varC=w_C+c_C+w_T;
            int expectvar = (double) pow(0.1, ((qvalueC)/10))/3;
            double pvalue = fishers_exact((int)expectvar*depth,depth,varC, depth);

                        if(depth >= nLayerMin  && qvalueC >= vQualMin && varC >=minread2 ){
                                double G2C=(double) varC/depth;
                                if(G2C>=minhomfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS");
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                                }

                        }else{
                                double G2C=(double) varC/depth;
                                printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                        }
                }   
        }   
//GG    
        if(genotypemaybe == "GG"){
                if(refbase == 'A'){//A>GG
            double qvalueG=(crqG>wsqG)?crqG:wsqG;
                        int depth=w_A+c_G+w_G+w_T+c_T+w_C+c_C;
                        int varG=w_G+c_G;
            int expectvar = (double) pow(0.1, ((qvalueG)/10))/3;
            double pvalue = fishers_exact((int)expectvar*depth,depth,varG, depth);

                        if(depth >= nLayerMin  && qvalueG >= vQualMin && varG >=minread2 ){
                                double A2G=(double) varG/depth;
                                double A2Gw=(double) w_G/Wtotaldepth;
                                if(A2G>=minhomfreq || A2Gw>=minhomfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS");
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                                }

                        }else{
                                double A2G=(double) varG/depth;
                                printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                        }
            
                }   
                if(refbase == 'T'){//T>GG
            double qvalueG=(crqG>wsqG)?crqG:wsqG;
                        int depth=c_T+w_T+c_G+w_G+w_A+w_C+c_C;
                        int varG=w_G+c_G+c_A;
            int expectvar = (double) pow(0.1, ((qvalueG)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varG, totaldepth);

                        if(totaldepth >= nLayerMin  && qvalueG >= vQualMin && varG >=minread2 ){
                                double T2G=(double) varG/totaldepth;
                                if(T2G>=minhomfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS");
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                                }

                        }else{
                                double T2G=(double) varG/totaldepth;
                                printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                        }
    
                } 

  
                if(refbase == 'C'){//C>GG
            double qvalueG=(crqG>wsqG)?crqG:wsqG;
                        int depth=c_T+w_T+c_G+w_G+w_A+w_C+c_C;
                        int varG=w_G+c_G+c_A;
                        int expectvar = (double) pow(0.1, ((qvalueG)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varG, totaldepth);

                        if(totaldepth >= nLayerMin  && qvalueG >= vQualMin && varG >=minread2 ){
                                double C2G=(double) varG/totaldepth;
                                if(C2G>=minhomfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS");
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                                }

                        }else{
                                double C2G=(double) varG/totaldepth;
                                printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                        }
    
                }   
                if(refbase == 'G'){
            return 0;
                }   
        }    
//CT 
    if(genotypemaybe == "CT"){
                if(refbase == 'A'){//A>CT
            double qvalueC=(wsqC>crqC)?wsqC:crqC;
                        double qvalueT=crqT;
                        int varC=w_C+c_C;
                        int varT=c_T;
                        //int depth=w_A+c_A+c_T;
                         int depth=totaldepth-w_T;
                         int expectvar = (double) pow(0.1, ((qvalueC+qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*depth,depth,varC+varT, depth);

                        if(depth >= nLayerMin  && qvalueC >= vQualMin && qvalueT>=vQualMin && varC >=minread2 && varT>=minread2 ){
                                double A2C=(double) varC/depth;
                                double A2T=(double) varT/Ctotaldepth;
                                if(A2C>=minhetfreq && A2T>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS");
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                                }

                        }else{
                                double A2C= (double) varC/depth;
                                double A2T= (double) varT/depth;
                                printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                        }

                }
                if(refbase == 'T'){ //T>CT
            double qvalueC=(wsqC>crqC)?wsqC:crqC;
                        int varC=w_C+c_C;
                        //int depth=c_G+w_G+w_A+w_T+c_T+w_C+c_C;
            int depth=totaldepth-w_T;
                        double T2C= (double) varC/depth;
                        int expectvar = (double) pow(0.1, ((qvalueC)/10))/3;
            double pvalue = fishers_exact((int)expectvar*depth,depth,varC, depth);

                        if(depth >= nLayerMin  && qvalueC >= vQualMin  && varC >=minread2 ){
                                if(T2C>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS");
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                                }

                        }else{
                            printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                        }
                }
                if(refbase == 'C'){ //C>CT
            double qvalueT=crqT;
                        int varT=c_T;
                        //int depth=c_G+w_G+w_A+w_T+c_T+w_C+c_C;
                        int depth=totaldepth-w_T;
                        double C2T=(double) varT/Ctotaldepth;
                        int expectvar = (double) pow(0.1, ((qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*Ctotaldepth,Ctotaldepth,varT, Ctotaldepth);

                        if(depth >= nLayerMin  && qvalueT >= vQualMin  && varT >=minread2 ){
                                if(C2T>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS");
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                                }

                        }else{
                            printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                        }

                }
                if(refbase == 'G'){//G>CT
            double qvalueC=(wsqC>crqC)?wsqC:crqC;
                        double qvalueT=crqT;
                        int varC=w_C+c_C;
                        int varT=c_T;
                        //int depth=w_A+c_A+c_T;
                        int depth=totaldepth-w_T;
                        int expectvar = (double) pow(0.1, ((qvalueC+qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*depth,depth,varC+varT, depth);

                        if(depth >= nLayerMin  && qvalueC >= vQualMin && qvalueT>=vQualMin && varC >=minread2 && varT>=minread2 ){
                                double G2C=(double) varC/depth;
                                double G2T=(double) varT/Ctotaldepth;
                                if(G2C>=minhetfreq && G2T>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS");
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                                }

                        }else{
                                double G2C= (double) varC/depth;
                                double G2T= (double) varT/depth;
                                printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                        }           
                }
        }   
//GT
    if(genotypemaybe == "GT"){
                if(refbase == 'A'){ //A>GT
            double qvalueG=(wsqG>crqG)?wsqG:crqG;
                        double qvalueT=crqT;
                        int varG=w_G+c_G;
                        int varT=c_T+w_T;
                        //int depth=w_A+c_A+c_T;
                        int depth=totaldepth-c_A;
                        int expectvar = (double) pow(0.1, ((qvalueG+qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varG+varT, totaldepth);

                        if(depth >= nLayerMin  && qvalueG >= vQualMin && qvalueT>=vQualMin && varG >=minread2 && varT>=minread2 ){
                                double A2G= (double) varG/depth;
                                double A2T= (double) varT/totaldepth;
                                if(A2G>=minhetfreq && A2T>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS");
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                                }

                        }else{
                                double A2G= (double) varG/depth;
                                double A2T= (double) varT/depth;
                                printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                        }

                }
                if(refbase == 'T'){//T>G
            double qvalueG=(wsqG>crqG)?wsqG:crqG;
                        int varG=w_G+c_G+c_A;
                        //int depth=w_A+c_A+c_T;
                        int depth=w_A+w_T+w_C+w_G+c_T+c_C+c_G;
                        int expectvar = (double) pow(0.1, ((qvalueG)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varG, totaldepth);

                        if(depth >= nLayerMin  && qvalueG >= vQualMin  && varG >=minread2 ){
                                double T2G= (double) varG/totaldepth;
                                if(T2G>=minhetfreq ){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS");
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                                }

                        }else{
                                double T2G= (double) varG/totaldepth;
                                
                        }   
                }
                if(refbase == 'C'){//C>GT
            double qvalueG=(wsqG>crqG)?wsqG:crqG;
                        double qvalueT=crqT;
                        int varG=w_G+c_G+c_A;
                        int varT=c_T;
                        //int depth=w_A+c_A+c_T;
                        int depth=totaldepth-w_T-c_A;
                        int expectvar = (double) pow(0.1, ((qvalueG+qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varG+varT, totaldepth);

                        if(depth >= nLayerMin  && qvalueG >= vQualMin && qvalueT>=vQualMin && varG >=minread2 && varT>=minread2 ){
                                double C2G= (double)varG/totaldepth;
                                double C2T= (double)varT/Ctotaldepth;
                                if(C2G>=minhetfreq && C2T>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS");
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                                }

                        }else{
                                double C2G= (double)varG/totaldepth;
                                double C2T= (double)varT/Ctotaldepth;
                                printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                        }

                }
                if(refbase == 'G'){//G>GT
            double qvalueT=crqT;
                        int varT=c_T+w_T;
                        //int depth=w_A+c_A+c_T;
                        int depth=totaldepth;
                        int expectvar = (double)pow(0.1, ((qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varT, totaldepth);

                        if(depth >= nLayerMin  && qvalueT >= vQualMin  && varT >=minread2 ){
                                double G2T= (double)varT/totaldepth;
                                if(G2T>=minhetfreq ){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS");
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                                }

                        }else{
                                double G2T= (double)varT/totaldepth;
                                printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                        }
                }
        }   
//CG
    if(genotypemaybe == "CG"){
                if(refbase == 'A'){//A>CG
            double qvalueG=(wsqG>crqG)?wsqG:crqG;
                        double qvalueC=(wsqT>crqT)?wsqT:crqT;
                        int varG=w_G+c_G;
                        int varC=c_C+w_C+w_T;
                        //int depth=w_A+c_A+c_T;
                        int depth=totaldepth-c_A;
                        int expectvar = (double)pow(0.1, ((qvalueG+qvalueC)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varG+varC, totaldepth);

                        if(totaldepth >= nLayerMin  && qvalueG >= vQualMin && qvalueC>=vQualMin && varG >=minread2 && varC>=minread2 ){
                                double A2G= (double)varG/depth;
                                double A2C= (double)varC/totaldepth;
                                if(A2G>=minhetfreq && A2C>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS");
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                                }
                        }else{
                                double A2G= (double)varG/depth;
                                double A2C= (double)varC/totaldepth;
                                printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                        }
    
                }
                if(refbase == 'T'){//T>CG
            double qvalueG=(wsqG>crqG)?wsqG:crqG;
                        double qvalueC=(wsqT>crqT)?wsqT:crqT;
                        int varG=w_G+c_G;
                        int varC=c_C+w_C;
                        //int depth=w_A+c_A+c_T;
                        int depth=totaldepth-w_T-c_A;
                        int expectvar = (double)pow(0.1, ((qvalueG+qvalueC)/10))/3;
            double pvalue = fishers_exact((int)expectvar*depth,depth,varG+varC, depth);

                        if(depth >= nLayerMin  && qvalueG >= vQualMin && qvalueC>=vQualMin && varG >=minread2 && varC>=minread2 ){
                                double T2G= (double)varG/depth;
                                double T2C= (double)varC/depth;
                                if(T2G>=minhetfreq && T2C>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS");
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                                }

                        }else{
                                double T2G= (double)varG/depth;
                                double T2C= (double)varC/depth;
                                printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                        }
                }
                if(refbase == 'C'){//C>CG
            double qvalueG=wsqG;
                        int varG=w_G+c_G+c_A;
                        //int depth=c_G+w_G+w_A+w_T+c_T+w_C+c_C;
                        int depth=totaldepth-w_T-c_A;
                        double C2G= (double)varG/totaldepth;
                        int expectvar = (double)pow(0.1, ((qvalueG)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varG, totaldepth);

                        if(depth >= nLayerMin  && qvalueG >= vQualMin  && varG >=minread2 ){
                                if(C2G>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS");
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                                }

                        }else{
                            printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                        }           
                }
                if(refbase == 'G'){//G>CG
            double qvalueC=crqC;
                        int varC=c_C+w_C+w_T;
                        //int depth=c_G+w_G+w_A+w_T+c_T+w_C+c_C;
                        int depth=totaldepth-w_T-c_A;
                        double G2C= (double)varC/totaldepth;
                        int expectvar = (double)pow(0.1, ((qvalueC)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varC, totaldepth);

                        if(depth >= nLayerMin  && qvalueC >= vQualMin  && varC >=minread2 ){
                                if(G2C>=minhetfreq){
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "PASS");
                                }else{
                                    printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                                }
                        }else{
                            printSNP(chrom, pos, refbase, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, "Low");
                        }
    
                }
        }   
    //unless(refbase=~/[ACGT]/i){
        return 0;
        //print "lines[0]\tlines[1]\t\.\trefbase\t0\tSuper\tNN\t\.\t".join("\t",@lines[3..6])."\n";
    //}

}