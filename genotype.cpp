#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bayes.hpp"
#include "fisher.hpp"
#include "tstudenttest.hpp"

extern int maxcover;
extern int minquali;
extern int mincover;
extern float minhetfreq;
extern float minhomfreq;
extern float pvalue_cutoff;
//int minhetfreq = minhetfreq;
extern int minread2;
extern int minvarread;
extern float minvarrate;
extern int othercover;
extern int printLowQ;
extern float tvalcutoff;

float ttest(int totaldepth, int ad){
    if(totaldepth<2){
        //fprintf(stderr, "T test, tc error.\n");
        return 1;
    }

    float array1[totaldepth];
    for(int i=0; i<totaldepth; i++) array1[i]=0;
    for(int i=0; i<ad; i++){
        array1[i]=1;
    }
    float array2[totaldepth];
    for(int i=0; i<totaldepth; i++) array2[i]=0;
    return tstudtest(array1, array2, totaldepth, totaldepth);
}

int printSNP(FILE* posFptr, char* chrom, int pos, char refbase, char* altbase, int genoqual, int w_A, int w_T, int w_C, int w_G, int c_A, int c_T, int c_C, int c_G, 
    int filter, int wsqA,int wsqT,int wsqC, int wsqG, int crqA, int crqT, int crqC, int crqG, int totaldepth, int adf, int adr, int ad, 
    int gt, double pvalue, double varfreq, double varfreq2, int var1, int var2, int efcoverage1, int efcoverage2){
//gt 1,2,3 0/1 1/2 1/1
    float ptstuval1=1;
    float ptstuval2=1;
    if(var1>=minread2) ptstuval1 = ttest(efcoverage1, var1);
    if(var2>=minread2) ptstuval2 = ttest(efcoverage2, var2);
    float ptstuval = ttest(totaldepth, ad);
    if( (filter==1 && pvalue < pvalue_cutoff && (ptstuval1<tvalcutoff || ptstuval2<tvalcutoff)) || 
     (ad>minvarread && (ptstuval1<tvalcutoff || ptstuval2<tvalcutoff)) ){//|| (filter==0 && ad>minvarread && (varfreq>minvarrate || varfreq2>minvarrate) && pvalue < pvalue_cutoff+0.01) ){
        //fprintf(stderr, "\n---highQ %d----\n", gt);
        if(gt==1) {
            fprintf(posFptr, "%s\t%d\t.\t%c\t%s\t%d\tPASS\tNS=1;DP=%d;ADF=%d;ADR=%d;AD=%d;\tGT:PVAL:BSD:BSQ:ALFR\t0/1:%.3f:%d,%d,%d,%d,%d,%d,%d,%d:%d,%d,%d,%d,%d,%d,%d,%d:%.2f\t%.3f\t%.3f\t%.3f\n", chrom, pos, refbase, altbase, genoqual, totaldepth, adf, adr, ad, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G,wsqA, wsqT, wsqC, wsqG, crqA, crqT, crqC, crqG, varfreq, ptstuval, ptstuval1, ptstuval2);
        }else if(gt==2){ // && ad>=8){
            if(varfreq<0.1 && varfreq2>0.5){
                fprintf(posFptr, "%s\t%d\t.\t%c\t%c\t%d\tPASS\tNS=1;DP=%d;ADF=%d;ADR=%d;AD=%d;\tGT:PVAL:BSD:BSQ:ALFR\t1/1:%.3f:%d,%d,%d,%d,%d,%d,%d,%d:%d,%d,%d,%d,%d,%d,%d,%d:%.2f\t%.3f\t%.3f\t%.3f\n", chrom, pos, refbase, altbase[2], genoqual, totaldepth, adf, adr, ad, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G,wsqA, wsqT, wsqC, wsqG, crqA, crqT, crqC, crqG, varfreq2, ptstuval, ptstuval1, ptstuval2);
            }else if(varfreq2<0.1 && varfreq>0.5){
                fprintf(posFptr, "%s\t%d\t.\t%c\t%c\t%d\tPASS\tNS=1;DP=%d;ADF=%d;ADR=%d;AD=%d;\tGT:PVAL:BSD:BSQ:ALFR\t1/1:%.3f:%d,%d,%d,%d,%d,%d,%d,%d:%d,%d,%d,%d,%d,%d,%d,%d:%.2f\t%.3f\t%.3f\t%.3f\n", chrom, pos, refbase, altbase[0], genoqual, totaldepth, adf, adr, ad, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G,wsqA, wsqT, wsqC, wsqG, crqA, crqT, crqC, crqG, varfreq, ptstuval, ptstuval1, ptstuval2);

            }
            else fprintf(posFptr, "%s\t%d\t.\t%c\t%s\t%d\tPASS\tNS=1;DP=%d;ADF=%d;ADR=%d;AD=%d;\tGT:PVAL:BSD:BSQ:ALFR\t1/2:%.3f:%d,%d,%d,%d,%d,%d,%d,%d:%d,%d,%d,%d,%d,%d,%d,%d:%.2f,%.2f\t%.3f\t%.3f\t%.3f\n", chrom, pos, refbase, altbase, genoqual, totaldepth, adf, adr, ad, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G,wsqA, wsqT, wsqC, wsqG, crqA, crqT, crqC, crqG, varfreq, varfreq2, ptstuval, ptstuval1, ptstuval2);
        //}else if(printLowQ && gt==2 && ad<8){
        //    fprintf(posFptr, "%s\t%d\t.\t%c\t%s\t%d\tLow\tNS=1;DP=%d;ADF=%d;ADR=%d;AD=%d;\tGT:PVAL:BSD:BSQ:ALFR\t1/2:%.3f:%d,%d,%d,%d,%d,%d,%d,%d:%d,%d,%d,%d,%d,%d,%d,%d:%.2f,%.2f\t%.3f\t%.3f\t%.3f\n", chrom, pos, refbase, altbase, genoqual, totaldepth, adf, adr, ad, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G,wsqA, wsqT, wsqC, wsqG, crqA, crqT, crqC, crqG, varfreq, varfreq2, ptstuval, ptstuval1, ptstuval2);
        }else if(gt==3) 
            fprintf(posFptr, "%s\t%d\t.\t%c\t%s\t%d\tPASS\tNS=1;DP=%d;ADF=%d;ADR=%d;AD=%d;\tGT:PVAL:BSD:BSQ:ALFR\t1/1:%.3f:%d,%d,%d,%d,%d,%d,%d,%d:%d,%d,%d,%d,%d,%d,%d,%d:%.2f\t%.3f\t%.3f\t%.3f\n", chrom, pos, refbase, altbase, genoqual, totaldepth, adf, adr, ad, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G,wsqA, wsqT, wsqC, wsqG, crqA, crqT, crqC, crqG, varfreq, ptstuval, ptstuval1, ptstuval2);
        else if(printLowQ) fprintf(stderr, "%s\t%d\t%c\t%s can not define the genotype type\n", chrom, pos, refbase, altbase);
        return 1;
    }else if(printLowQ){ //Low Qual
       //fprintf(stderr, "\n---lowQ %d----\n", gt);
	   //if((filter==1 || pvalue < pvalue_cutoff)){
        if(gt==1)
            fprintf(posFptr, "%s\t%d\t.\t%c\t%s\t%d\tLow\tNS=1;DP=%d;ADF=%d;ADR=%d;AD=%d;\tGT:PVAL:BSD:BSQ:ALFR\t0/1:%.3f:%d,%d,%d,%d,%d,%d,%d,%d:%d,%d,%d,%d,%d,%d,%d,%d:%.2f\t%.3f\t%.3f\t%.3f\n", chrom, pos, refbase, altbase, genoqual, totaldepth, adf, adr, ad, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G,wsqA, wsqT, wsqC, wsqG, crqA, crqT, crqC, crqG, varfreq, ptstuval, ptstuval1, ptstuval2);
        else if(gt==2)
            fprintf(posFptr, "%s\t%d\t.\t%c\t%s\t%d\tLow\tNS=1;DP=%d;ADF=%d;ADR=%d;AD=%d;\tGT:PVAL:BSD:BSQ:ALFR\t1/2:%.3f:%d,%d,%d,%d,%d,%d,%d,%d:%d,%d,%d,%d,%d,%d,%d,%d:%.2f,%.2f\t%.3f\t%.3f\t%.3f\n", chrom, pos, refbase, altbase, genoqual, totaldepth, adf, adr, ad, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G,wsqA, wsqT, wsqC, wsqG, crqA, crqT, crqC, crqG, varfreq, varfreq2, ptstuval, ptstuval1, ptstuval2);
        else if(gt==3)
            fprintf(posFptr, "%s\t%d\t.\t%c\t%s\t%d\tLow\tNS=1;DP=%d;ADF=%d;ADR=%d;AD=%d;\tGT:PVAL:BSD:BSQ:ALFR\t1/1:%.3f:%d,%d,%d,%d,%d,%d,%d,%d:%d,%d,%d,%d,%d,%d,%d,%d:%.2f\t%.3f\t%.3f\t%.3f\n", chrom, pos, refbase, altbase, genoqual, totaldepth, adf, adr, ad, pvalue, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G,wsqA, wsqT, wsqC, wsqG, crqA, crqT, crqC, crqG, varfreq, ptstuval, ptstuval1, ptstuval2);
        else fprintf(stderr, "%s\t%d\t%c\t%s can not define the genotype type\n", chrom, pos, refbase, altbase);
	//}
    }
    return 0;

}

int genotype(FILE* posFptr, int wsqA,int wsqT,int wsqC, int wsqG, int crqA, int crqT, int crqC, int crqG, char refbase, unsigned int pos, char* chrom,
    int w_A, int w_T, int w_C, int w_G, int c_A, int c_T, int c_C, int c_G, std::string genotypemaybe, double genoqual)
{
    //print "lines[0]\tlines[1]\trefbase\tgenotypemaybe\n";
    int filpass=0;
    unsigned int totaldepth=w_A+w_T+w_C+w_G+c_A+c_T+c_C+c_G;
    if(totaldepth > maxcover){
        return 0;
    }
    char altbase[4];
    int adf=0; int adr=0; int ad=0;

    int Ctotaldepth=c_C+c_A+c_T+c_G;
    int Wtotaldepth=w_A+w_T+w_C+w_G;
    if(genotypemaybe == "AA"){//genotypeis AA
        if(refbase == 'A'){
            return 0;   
        }
        altbase[0]='A';altbase[1]='\0';

        if(refbase == 'T'){//T>AA
            double qvalue=(wsqA>crqA)?wsqA:crqA;
            double expectvar = (int) pow(0.1, (qvalue/10))/3;
            int var=w_A+c_A;
            adf=w_A;
            adr=c_A;
            ad=var;
            if(w_G+c_G+w_C+c_C > othercover) return 0;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,var, totaldepth);
            
            if(totaldepth >= mincover  && qvalue >= minquali && var >=minread2 ){
                //sprintf("%.2f", f)
                double T2A= (double) var/totaldepth;

                if(T2A>=minhomfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, T2A, 0, var, 0, totaldepth, totaldepth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, T2A, 0, var, 0, totaldepth, totaldepth);
                }
                
            }else{
                double T2A=(double) var/totaldepth;
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, T2A, 0, var, 0, totaldepth, totaldepth);
            }                                   
        }
        if(refbase == 'C'){//C>AA
            double qvalue=(wsqA>crqA)?wsqA:crqA;
            //int depth=w_C+c_C+w_A+c_A;
            int var=w_A+c_A;   
            adf=w_A;
            adr=c_A;
            ad=var;
            if(w_G+c_G > othercover) return 0;
            double expectvar = (double) pow(0.1, (qvalue/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,var, totaldepth);
            if(totaldepth >= mincover  && qvalue >= minquali && var >=minread2 ){
                double C2A=(double)var/totaldepth;
                if(C2A>=minhomfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, C2A, 0, var, 0, totaldepth, totaldepth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, C2A, 0, var, 0, totaldepth, totaldepth);
                }   
    
            }else{
                double C2A= (double) var/totaldepth;
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, C2A, 0, var, 0, totaldepth, totaldepth);
            }                   
        }

        if(refbase == 'G'){//G>AA
            double qvalue=wsqA;
            int depth=w_G+w_A;
            int var=w_A;
            adf=w_A;
            adr=0;
            ad=var;
            if(w_T+c_T+w_C+c_C > othercover) return 0;
            double G2A;
            if(depth>0){
                G2A=(double) var/depth;
            }else{
                G2A=0;
            }
            double expectvar = (double) pow(0.1, (qvalue/10))/3;
            double pvalue = fishers_exact((int)expectvar*depth,depth,var, depth);

            if(depth >= mincover  && qvalue >= minquali && var >=minread2 ){
                if(G2A>=minhomfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, G2A, 0, var, 0, depth, depth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, G2A, 0, var, 0, depth, depth); 
                }

            }else{ // HHHHH 01
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, G2A, 0, var, 0, depth, depth);
            }
        }
    }

    if(genotypemaybe == "AT"){ //genotype is AT
        if(w_G+c_G+w_C+c_C > othercover) return 0;
        if(refbase == 'A'){//A>AT
            altbase[0]='T';altbase[1]='\0';

            double qvalue=(wsqT>crqT)?wsqT:crqT;
            //int depth=w_A+w_T+c_A+c_T;
            int var=w_T+c_T;
            adf=w_T;
            adr=c_T;
            ad=var;
            double expectvar = (double) pow(0.1, (qvalue/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,var, totaldepth);
            if(totaldepth >= mincover  && qvalue >= minquali && var >=minread2 ){
                double A2T=(double) var/totaldepth;
                if(A2T>=minhetfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, A2T, 0, var, 0, totaldepth, totaldepth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, A2T, 0, var, 0, totaldepth, totaldepth); 
                }
            }else{
                double A2T=(double) var/totaldepth;
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, A2T, 0, var, 0, totaldepth, totaldepth);
            } 
        }   

        if(refbase == 'T'){//T>AT 
            altbase[0]='A';altbase[1]='\0';        
            double qvalue=(wsqA>crqA)?wsqA:crqA;
            //int depth=w_A+w_T+c_A+c_T;
            int var=c_A+w_A;
            adf=w_A;
            adr=c_A;
            ad=var;
            double expectvar = (double) pow(0.1, (qvalue/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,var, totaldepth);

            if(totaldepth >= mincover  && qvalue >= minquali && var >=minread2 ){
                double T2A=(double) var/totaldepth;
                if(T2A>=minhetfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, T2A, 0, var, 0, totaldepth, totaldepth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, T2A, 0, var, 0, totaldepth, totaldepth);
                }

                }else{
                    double T2A=(double) var/totaldepth;
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, T2A, 0, var, 0, totaldepth, totaldepth); 
                }
        }   
        if(refbase == 'C'){//C>AT
            altbase[0]='A';altbase[1]=',';altbase[2]='T';altbase[3]='\0';
            double qvalueA=(wsqA>crqA)?wsqA:crqA;
            double qvalueT=crqT;
            int varA=w_A+c_A;
            int varT=c_T;
            adf=w_A;
            adr=c_A+c_T;
            ad=varA+varT;
            int depth=w_A+c_A+c_T;
            double expectvar = (double) pow(0.1, ((qvalueA+qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varA+varT, totaldepth);

            if(totaldepth >= mincover  && qvalueA >= minquali && qvalueT>=minquali && varA >=minread2 && varT>=minread2 ){
                double C2A=(double) varA/totaldepth;
                double C2T=(double) varT/Ctotaldepth;
                if(C2A>=minhetfreq && C2T>=minhetfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, C2A, C2T, varA, varT, totaldepth, Ctotaldepth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, C2A, C2T, varA, varT, totaldepth, Ctotaldepth); 
                }   
            }else{
                double C2A=(double) varT/totaldepth;
                double C2T=(double) varT/totaldepth;
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, C2A, C2T, varA, varT, totaldepth, Ctotaldepth);
            }       
            
        }   
        if(refbase == 'G'){//G>AT 
            altbase[0]='A';altbase[1]=',';altbase[2]='T';altbase[3]='\0'; 
            double qvalueA=wsqA;
            double qvalueT=(wsqT>crqT)?wsqT:crqT;
            int varA=w_A;
            int varT=c_T+w_T;
            adf=w_A+w_T;
            adr=c_T;
            ad=varA+varT;
            int depth=w_A+c_T+w_T;
            double expectvar = (double) pow(0.1, ((qvalueA+qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varA+varT, totaldepth);
            //qw
            if(totaldepth >= mincover  && qvalueA >= minquali && qvalueT>=minquali && varA >=minread2 && varT>=minread2 ){
                double G2A=(double) varA/Wtotaldepth;
                double G2T=(double) varT/totaldepth;
                if(G2A>=minhetfreq && G2T>=minhetfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, G2A, G2T, varA, varT, Wtotaldepth, totaldepth); 
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, G2A, G2T, varA, varT, Wtotaldepth, totaldepth); 
                }

            }else{
                double G2A=(double) varA/totaldepth;
                double G2T=(double) varT/totaldepth;
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, G2A, G2T, varA, varT, Wtotaldepth, totaldepth);
            }
                
        }                   
    }
//AC 
    if(genotypemaybe == "AC"){ //genotype is AC
        if(c_T+w_G+c_G > othercover) return 0;
        if(refbase == 'A'){ //A>AC
            altbase[0]='C';altbase[1]='\0'; 
            double qvalueC=(wsqC>crqC)?wsqC:crqC;
            if(qvalueC < wsqT) qvalueC=wsqT;
            int varC=w_C+c_C+w_T;
            adf=w_C+w_T;
            adr=c_C;
            ad=varC;
            //int depth=w_A+c_C+w_C;
            double A2C=(double) varC/totaldepth;
            double expectvar = (double) pow(0.1, (qvalueC/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varC, totaldepth);

            if(totaldepth >= mincover  && qvalueC >= minquali  && varC >=minread2 ){
                if(A2C>=minhetfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, A2C, 0, varC, 0, totaldepth, totaldepth); 
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, A2C, 0, varC, 0, totaldepth, totaldepth); 
                }   

            }else{
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, A2C, 0, varC, 0, totaldepth, totaldepth); 
            }               
        }
        if(refbase == 'T'){//T>AC ////////////////////////modify totaldepth
            altbase[0]='A';altbase[1]=',';altbase[2]='C';altbase[3]='\0'; 
            double qvalueA=(wsqA>crqA)?wsqA:crqA;
            double qvalueC=(wsqC>crqC)?wsqC:crqC;
            int varA=w_A+c_A;
            int varC=c_C+w_C;
            int depth=c_T+c_C+w_C+w_A+c_A+w_G+c_G;
            adf=w_A+w_C;
            adr=c_A+c_C;
            ad=varA+varC;
            //int depth=w_A+c_C+w_C+c_A+c_T;
            double T2A=(double) varA/totaldepth;
            double T2C=(double) varC/depth;
            double expectvar = (double) pow(0.1, ((qvalueC+qvalueA)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varC+varA, totaldepth);

            if(totaldepth >= mincover  && qvalueA >= minquali && qvalueC>=minquali && varA >=minread2 && varC>=minread2 ){
                if(T2A>=minhetfreq && T2C>=minhetfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, T2A, T2C, varA, varC, totaldepth, depth);  
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, T2A, T2C, varA, varC, totaldepth, depth);
                }   
            }else{
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, T2A, T2C, varA, varC, totaldepth, depth);
            }       
        }
        if(refbase == 'C'){//C>AC
            altbase[0]='A';altbase[1]='\0';
            double qvalueA=(wsqA>crqA)?wsqA:crqA;
            int varA=w_A+c_A;
            int depth=w_A+c_C+w_C;
            adf=w_A;
            adr=c_A;
            ad=varA;
            double C2A=(double) varA/totaldepth;
            double expectvar = (double) pow(0.1, ((qvalueA)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varA, totaldepth);

            if(totaldepth >= mincover  && qvalueA >= minquali  && varA >=minread2 ){
                if(C2A>=minhetfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, C2A, 0, varA, 0, totaldepth, totaldepth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, C2A, 0, varA, 0, totaldepth, totaldepth);
                }
            }else{
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, C2A, 0, varA, 0, totaldepth, totaldepth);
            }           
        }
        if(refbase == 'G'){//G>AC
            altbase[0]='A';altbase[1]=',';altbase[2]='C';altbase[3]='\0'; 
            double qvalueA=wsqA;
            double qvalueC=(wsqC>crqC)?wsqC:crqC;
            int varA=w_A;
            int varC=c_C+w_C;
            adf=w_A+w_C;
            adr=c_C;
            ad=varA+varC;
            int depth=w_A+c_C+w_C+c_A;
            double G2A=(double) varA/Wtotaldepth;
            double G2C=(double) varC/totaldepth;
            double expectvar = (double) pow(0.1, ((qvalueC+qvalueA)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varC+varA, totaldepth);

            if(totaldepth >= mincover  && qvalueA >= minquali && qvalueC>=minquali && varA >=minread2 && varC>=minread2 ){
                if(G2A>=minhetfreq && G2C>=minhetfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, G2A, G2C, varA, varC, Wtotaldepth, totaldepth);
                }else{
                   filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, G2A, G2C, varA, varC, Wtotaldepth, totaldepth);
                }

            }else{
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, G2A, G2C, varA, varC, Wtotaldepth, totaldepth);
            }
        }
    }   
    
//AG 
    if(genotypemaybe == "AG"){
        if(w_T+c_T+w_C+c_C > othercover) return 0;
        if(refbase == 'A'){//A>AG
            altbase[0]='G';altbase[1]='\0'; 
            double qvalueG=(wsqG>crqG)?wsqG:crqG;
            int varG=w_G+c_G;
            adf=w_G;
            adr=c_G;
            ad=varG;
            int depth=w_A+c_G+w_G+w_T+c_T+w_C+c_C;
            double A2G=(double) varG/depth;
            double expectvar = (double) pow(0.1, ((qvalueG)/10))/3;
            double pvalue = fishers_exact((int)expectvar*depth,depth,varG, depth);
            if(depth >= mincover  && qvalueG >= minquali  && varG >=minread2 ){
                if(A2G>=minhetfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, A2G, 0, varG, 0, depth, totaldepth); 
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, A2G, 0, varG, 0, depth, totaldepth); 
                }   

            }else{
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, A2G, 0, varG, 0, depth, totaldepth);
            }               
        }
        if(refbase == 'T'){//T>AG
            altbase[0]='A';altbase[1]=',';altbase[2]='G';altbase[3]='\0'; 
            double qvalueA=wsqA;
            double qvalueG=(wsqG>crqG)?wsqG:crqG;
            int varA=w_A;
            int varG=c_G+w_G;
            adf=w_A+w_G;
            adr=c_G;
            ad=varA+varG;
            int depth=w_A+c_G+w_G+c_A+c_T+w_T;
            double T2A=(double) varA/Wtotaldepth;
            double T2G=(double) varG/totaldepth;
            double expectvar = (double) pow(0.1, ((qvalueG+qvalueA)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varG+varA, totaldepth);
            if(totaldepth >= mincover  && qvalueA >= minquali && qvalueG>=minquali && varA >=minread2 && varG>=minread2 ){
              if(T2A>=minhetfreq && T2G>=minhetfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, T2A, T2G, varA, varG, Wtotaldepth, totaldepth); 
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, T2A, T2G, varA, varG, Wtotaldepth, totaldepth);
                }

                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, T2A, T2G, varA, varG, Wtotaldepth, totaldepth);
                }
    
            }
        if(refbase == 'C'){//C>AG
            altbase[0]='A';altbase[1]=',';altbase[2]='G';altbase[3]='\0'; //2
            double qvalueA=wsqA;
            double qvalueG=(wsqG>crqG)?wsqG:crqG;
            int varA=w_A;
            int varG=c_G+w_G;
            adf=w_A+w_G;
            adr=c_G;
            ad=varA+varG;
            int depth=w_A+c_C+w_C+c_A+c_T+w_T;
            double C2A=(double) varA/Wtotaldepth;
            double C2G=(double) varG/totaldepth;
            double expectvar = (double) pow(0.1, ((qvalueG+qvalueA)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varG+varA, totaldepth);
            if(totaldepth >= mincover  && qvalueA >= minquali && qvalueG>=minquali && varA >=minread2 && varG>=minread2 ){
                if(C2A>=minhetfreq && C2G>=minhetfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, C2A, C2G, varA, varG, Wtotaldepth, totaldepth); 
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, C2A, C2G, varA, varG, Wtotaldepth, totaldepth); 
                }

            }else{
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, C2A, C2G, varA, varG, Wtotaldepth, totaldepth); 
            }
    
        }
        if(refbase == 'G'){//G>AG
            altbase[0]='A';altbase[1]='\0'; //1
            double qvalueA=wsqA;
            int varA=w_A;
            adf=w_A;
            adr=0;
            ad=varA;
            //int depth=w_G+w_A+w_T+w_C;
            //int depth=c_G+w_G+w_A+w_T+c_T+w_C+c_C;
            double G2A=(double) varA/Wtotaldepth;
            double expectvar = (double) pow(0.1, ((qvalueA)/10))/3;
            double pvalue = fishers_exact((int)expectvar*Wtotaldepth,Wtotaldepth,varA, Wtotaldepth);
            if(Wtotaldepth >= mincover  && qvalueA >= minquali  && varA >=minread2 ){
                if(G2A>=minhetfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, G2A, 0, varA, 0, Wtotaldepth, totaldepth); 
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, G2A, 0, varA, 0, Wtotaldepth, totaldepth); 
                }

            }else{
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, G2A, 0, varA, 0, Wtotaldepth, totaldepth); 
            }
    
        }
    }   
//TT
    if(genotypemaybe == "TT"){
        altbase[0]='T';altbase[1]='\0'; //3
        if(refbase == 'A'){//A>TT
            if(w_G+c_G+w_C+c_C > othercover) return 0;
            double qvalueT=(crqT>wsqT)?crqT:wsqT;
            //int depth=w_A+c_A+c_T;
            int varT=c_T+w_T;
            adf=w_T;
            adr=c_T;
            ad=varT;
            double expectvar = (double) pow(0.1, ((qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varT, totaldepth);

            if(totaldepth >= mincover  && qvalueT >= minquali && varT >=minread2 ){
                    double A2T=(double) varT/totaldepth;
                    //print "\nA2T\t", A2T, "\t", varT, "\t", totaldepth;
                    if(A2T>=minhomfreq){
                        filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, A2T, 0, varT, 0, totaldepth, totaldepth); 
                    }else{
                        filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, A2T, 0, varT, 0, totaldepth, totaldepth); 
                    }

            }else{
                    double A2T=(double)varT/totaldepth;
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, A2T, 0, varT, 0, totaldepth, totaldepth);
            }
        }
        if(refbase == 'T'){
            return 0; 
        }
        if(refbase == 'C'){//C>TT
            if(w_G+c_G+w_A+c_A > othercover) return 0;
            double qvalueT=crqT;
            int depth=c_T+c_C;
            int varT=c_T;
            adf=0;
            adr=c_T;
            ad=varT;
            double C2T;
            if(depth>0){
                C2T=(double) varT/Ctotaldepth;
            }else{
                C2T=0;
                
            }
            double expectvar = (double) pow(0.1, ((qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*Ctotaldepth,Ctotaldepth,varT, Ctotaldepth);
            
            if(totaldepth >= mincover  && qvalueT >= minquali && varT >=minread2 ){
                if(C2T>=minhomfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, C2T, 0, varT, 0, Ctotaldepth, totaldepth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, C2T, 0, varT, 0, Ctotaldepth, totaldepth);
                }

            }else{
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, C2T, 0, varT, 0, Ctotaldepth, totaldepth); 
            }
    
        }
        if(refbase == 'G'){//G>TT
            if(w_A+w_C+c_C > othercover) return 0;
            double qvalueT=(crqT>wsqT)?crqT:wsqT;
            int depth=w_G+c_G+c_T+w_T+c_A;
            int varT=c_T+w_T;
            adf=w_T;
            adr=c_T;
            ad=varT;
            double expectvar = (double) pow(0.1, ((qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varT, totaldepth);
            if(totaldepth >= mincover  && qvalueT >= minquali && varT >=minread2 ){
                double G2T=(double) varT/totaldepth;
                if(G2T>=minhomfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, G2T, 0, varT, 0, totaldepth, totaldepth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, G2T, 0, varT, 0, totaldepth, totaldepth); 
                }

                }else{
                    double G2T=(double) varT/totaldepth;
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, G2T, 0, varT, 0, totaldepth, totaldepth);
                }
    
            }
        }  

//CC     
    if(genotypemaybe == "CC"){
        altbase[0]='C';altbase[1]='\0'; //3
        if(refbase == 'A'){//A>CC
            if(c_T+w_G+c_G > othercover) return 0;
            double qvalueC=(crqC>wsqC)?crqC:wsqC;
            int depth=w_A+c_A+c_C+w_C+w_G+w_G+c_T;
            int varC=w_C+c_C;
            adf=w_C;
            adr=c_C;
            ad=varC;
            double expectvar = (double) pow(0.1, ((qvalueC)/10))/3;
            double pvalue = fishers_exact((int)expectvar*depth,depth,varC, depth);

            if(depth >= mincover  && qvalueC >= minquali && varC >=minread2 ){
                double A2C=(double) varC/depth;
                if(A2C>=minhomfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, A2C, 0, varC, 0, depth, totaldepth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, A2C, 0, varC, 0, depth, totaldepth);
                }

            }else{
                double A2C=(double) varC/depth;
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, A2C, 0, varC, 0, depth, totaldepth);
            }   
        }   
        if(refbase == 'T'){//T>CC
            if(w_A+c_A+w_G+c_G > othercover) return 0;
            double qvalueC=(crqC>wsqC)?crqC:wsqC;
            int depth=c_T+c_C+w_C+w_A+c_A+w_G+c_G;
            int varC=w_C+c_C;
            adf=w_C;
            adr=c_C;
            ad=varC;
            double expectvar = (double) pow(0.1, ((qvalueC)/10))/3;
            double pvalue = fishers_exact((int)expectvar*depth,depth,varC, depth);

            if(totaldepth >= mincover  && qvalueC >= minquali && varC >=minread2 ){
                double T2C=(double) varC/depth;
                double T2Cc=(double) c_C/Ctotaldepth;
                if(T2C>=minhomfreq || T2Cc>=minhomfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, T2C, 0, c_C, varC, Ctotaldepth, depth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, T2C, 0, c_C, varC, Ctotaldepth, depth);
                }

            }else{
                double T2C=(double) varC/depth;
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, T2C, 0, c_C, varC, Ctotaldepth, depth);
            }

        }   
        if(refbase == 'C'){//C>CC  
            return 0;
        }   
        if(refbase == 'G'){//G>CC
            if(c_T+w_A > othercover) return 0;
            double qvalueC=(crqC>wsqC)?crqC:wsqC;
            int depth=w_G+c_G+c_C+w_C+w_A+c_A+c_T+w_T;
            int varC=w_C+c_C+w_T;
            adf=w_C+w_T;
            adr=c_C;
            ad=varC;
            double expectvar = (double) pow(0.1, ((qvalueC)/10))/3;
            double pvalue = fishers_exact((int)expectvar*depth,depth,varC, depth);

            if(depth >= mincover  && qvalueC >= minquali && varC >=minread2 ){
                double G2C=(double) varC/depth;
                if(G2C>=minhomfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, G2C, 0, varC, 0, depth, totaldepth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, G2C, 0, varC, 0, depth, totaldepth);
                }

            }else{
                double G2C=(double) varC/depth;
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, G2C, 0, varC, 0, depth, totaldepth);
            }
        }   
    }   
//GG    
    if(genotypemaybe == "GG"){
        altbase[0]='G';altbase[1]='\0'; //3
        if(refbase == 'A'){//A>GG
            if(w_T+c_T+w_C+c_C > othercover) return 0;
            double qvalueG=(crqG>wsqG)?crqG:wsqG;
            int depth=w_A+c_G+w_G+w_T+c_T+w_C+c_C;
            int varG=w_G+c_G;
            adf=w_G;
            adr=c_G;
            ad=varG;
            double expectvar = (double) pow(0.1, ((qvalueG)/10))/3;
            double pvalue = fishers_exact((int)expectvar*depth,depth,varG, depth);

            if(depth >= mincover  && qvalueG >= minquali && varG >=minread2 ){
                double A2G=(double) varG/depth;
                double A2Gw=(double) w_G/Wtotaldepth;
                if(A2G>=minhomfreq || A2Gw>=minhomfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, A2G, 0, w_G, varG, Wtotaldepth, depth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, A2G, 0, w_G, varG, Wtotaldepth, depth);
                }

            }else{
                double A2G=(double) varG/depth;
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, A2G, 0, w_G, varG, Wtotaldepth, depth);
            }
            
        }   
        if(refbase == 'T'){//T>GG
            if(w_A+w_C+c_C > othercover) return 0;
            double qvalueG=(crqG>wsqG)?crqG:wsqG;
            if(qvalueG<crqA) qvalueG=crqA;
            int depth=c_T+w_T+c_G+w_G+w_A+w_C+c_C;
            int varG=w_G+c_G+c_A;
            adf=w_G;
            adr=c_G+c_A;
            ad=varG;
            double expectvar = (double) pow(0.1, ((qvalueG)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varG, totaldepth);

            if(totaldepth >= mincover  && qvalueG >= minquali && varG >=minread2 ){
                double T2G=(double) varG/totaldepth;
                if(T2G>=minhomfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, T2G, 0, varG, 0, totaldepth, totaldepth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, T2G, 0, varG, 0, totaldepth, totaldepth);
                }
            }else{
                double T2G=(double) varG/totaldepth;
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, T2G, 0, varG, 0, totaldepth, totaldepth);
            }
    
        } 
        if(refbase == 'C'){//C>GG
            if(c_T+w_A > othercover) return 0;
            double qvalueG=(crqG>wsqG)?crqG:wsqG;
            if(qvalueG<crqA) qvalueG=crqA;
            int depth=c_T+w_T+c_G+w_G+w_A+w_C+c_C;
            int varG=w_G+c_G+c_A;
            adf=w_G;
            adr=c_G+c_A;
            ad=varG;
            double expectvar = (double) pow(0.1, ((qvalueG)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varG, totaldepth);

            if(totaldepth >= mincover  && qvalueG >= minquali && varG >=minread2 ){
                double C2G=(double) varG/totaldepth;
                if(C2G>=minhomfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, C2G, 0, varG, 0, totaldepth, totaldepth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, C2G, 0, varG, 0, totaldepth, totaldepth);
                }

            }else{
                double C2G=(double) varG/totaldepth;
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 3, pvalue, C2G, 0, varG, 0, totaldepth, totaldepth);
            }
    
        }   
        if(refbase == 'G'){
            return 0;
        }   
    }    
//CT 
    if(genotypemaybe == "CT"){
        if(w_A+c_A+w_G+c_G > othercover) return 0;
        if(refbase == 'A'){//A>CT
            altbase[0]='C';altbase[1]=',';altbase[2]='T';altbase[3]='\0'; //2
            double qvalueC=(wsqC>crqC)?wsqC:crqC;
            double qvalueT=crqT;
            int varC=w_C+c_C;
            int varT=c_T;
            adf=w_C;
            adr=c_C+c_T;
            ad=varC+varT;
            //int depth=w_A+c_A+c_T;
            int depth=totaldepth-w_T;
            double expectvar = (double) pow(0.1, ((qvalueC+qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*depth,depth,varC+varT, depth);

            if(totaldepth >= mincover  && qvalueC >= minquali && qvalueT>=minquali && varC >=minread2 && varT>=minread2 ){
                double A2C=(double) varC/depth;
                double A2T=(double) varT/Ctotaldepth;
                if(A2C>=minhetfreq && A2T>=minhetfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, A2C, A2T, varC, varT, depth, Ctotaldepth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, A2C, A2T, varC, varT, depth, Ctotaldepth);
                }

            }else{
                double A2C= (double) varC/depth;
                double A2T= (double) varT/depth;
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, A2C, A2T, varC, varT, depth, Ctotaldepth);
            }

        }
        if(refbase == 'T'){ //T>CT
            altbase[0]='C';altbase[1]='\0'; //1
            double qvalueC=(wsqC>crqC)?wsqC:crqC;
            int varC=w_C+c_C;
            adf=w_C;
            adr=c_C;
            ad=varC;
            //int depth=c_G+w_G+w_A+w_T+c_T+w_C+c_C;
            int depth=totaldepth-w_T;
            double T2C= (double) varC/depth;
            double expectvar = (double) pow(0.1, ((qvalueC)/10))/3;
            double pvalue = fishers_exact((int)expectvar*depth,depth,varC, depth);

            if(totaldepth >= mincover  && qvalueC >= minquali  && varC >=minread2 ){
                if(T2C>=minhetfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, T2C, 0, varC, 0, depth, totaldepth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, T2C, 0, varC, 0, depth, totaldepth);
                }

            }else{
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, T2C, 0, varC, 0, depth, totaldepth);
            }
        }
        if(refbase == 'C'){ //C>CT
            altbase[0]='T';altbase[1]='\0'; //1
            double qvalueT=crqT;
            int varT=c_T;
            adf=0;
            adr=c_T;
            ad=varT;
            //int depth=c_G+w_G+w_A+w_T+c_T+w_C+c_C;
            int depth=totaldepth-w_T;
            double C2T=(double) varT/Ctotaldepth;
            double expectvar = (double) pow(0.1, ((qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*Ctotaldepth,Ctotaldepth,varT, Ctotaldepth);

            if(totaldepth >= mincover  && qvalueT >= minquali  && varT >=minread2 ){
                if(C2T>=minhetfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, C2T, 0, varT, 0, Ctotaldepth, totaldepth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, C2T, 0, varT, 0, Ctotaldepth, totaldepth);
                }

            }else{
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, C2T, 0, varT, 0, Ctotaldepth, totaldepth);
            }

        }
        if(refbase == 'G'){//G>CT
            altbase[0]='C';altbase[1]=',';altbase[2]='T';altbase[3]='\0'; //2
            double qvalueC=(wsqC>crqC)?wsqC:crqC;
            double qvalueT=crqT;
            int varC=w_C+c_C;
            int varT=c_T;
            adf=w_C;
            adr=c_C+c_T;
            ad=varC+varT;
            //int depth=w_A+c_A+c_T;
            int depth=totaldepth-w_T;
            double expectvar = (double) pow(0.1, ((qvalueC+qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*depth,depth,varC+varT, depth);

            if(totaldepth >= mincover  && qvalueC >= minquali && qvalueT>=minquali && varC >=minread2 && varT>=minread2 ){
                    double G2C=(double) varC/depth;
                    double G2T=(double) varT/Ctotaldepth;
                    if(G2C>=minhetfreq && G2T>=minhetfreq){
                        filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, G2C, G2T, varC, varT, depth, Ctotaldepth);
                    }else{
                        filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, G2C, G2T, varC, varT, depth, Ctotaldepth);
                    }

            }else{
                double G2C= (double) varC/depth;
                double G2T= (double) varT/depth;
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, G2C, G2T, varC, varT, depth, Ctotaldepth);
            }           
        }
    }   
//GT
    if(genotypemaybe == "GT"){
        if(w_A+c_C+w_C > othercover) return 0;
        if(refbase == 'A'){ //A>GT
            altbase[0]='G';altbase[1]=',';altbase[2]='T';altbase[3]='\0'; //2
            double qvalueG=(wsqG>crqG)?wsqG:crqG;
            double qvalueT=(wsqT>crqT)?wsqT:crqT;
            int varG=w_G+c_G;
            int varT=c_T+w_T;
            adf=w_G+w_T;
            adr=c_G+c_T;
            ad=varG+varT;
            //int depth=w_A+c_A+c_T;
            int depth=totaldepth-c_A;
            double expectvar = (double) pow(0.1, ((qvalueG+qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varG+varT, totaldepth);

            if(totaldepth >= mincover  && qvalueG >= minquali && qvalueT>=minquali && varG >=minread2 && varT>=minread2 ){
                double A2G= (double) varG/depth;
                double A2T= (double) varT/totaldepth;
                if(A2G>=minhetfreq && A2T>=minhetfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, A2G, A2T, varG, varT, depth, totaldepth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, A2G, A2T, varG, varT, depth, totaldepth);
                }

            }else{
                double A2G= (double) varG/depth;
                double A2T= (double) varT/depth;
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, A2G, A2T, varG, varT, depth, totaldepth);
            }

        }
        if(refbase == 'T'){//T>G
            altbase[0]='G';altbase[1]='\0'; //1
            double qvalueG=(wsqG>crqG)?wsqG:crqG;
            if(qvalueG<crqA) qvalueG=crqA;
            int varG=w_G+c_G+c_A;
            adf=w_G;
            adr=c_G+c_A;
            ad=varG;
            //int depth=w_A+c_A+c_T;
            int depth=w_A+w_T+w_C+w_G+c_T+c_C+c_G;
            double expectvar = (double) pow(0.1, ((qvalueG)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varG, totaldepth);

            if(totaldepth >= mincover  && qvalueG >= minquali  && varG >=minread2 ){
                double T2G= (double) varG/totaldepth;
                if(T2G>=minhetfreq ){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, T2G, 0, varG, 0, totaldepth, totaldepth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, T2G, 0, varG, 0, totaldepth, totaldepth);
                }

            }else{
                    double T2G= (double) varG/totaldepth;
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, T2G, 0, varG, 0, totaldepth, totaldepth); 
            }   
        }
        if(refbase == 'C'){//C>GT
            altbase[0]='G';altbase[1]=',';altbase[2]='T';altbase[3]='\0'; //2
            double qvalueG=(wsqG>crqG)?wsqG:crqG;
            if(qvalueG<crqA) qvalueG=crqA;
            double qvalueT=crqT;
            int varG=w_G+c_G+c_A;
            int varT=c_T;
            adf=w_G;
            adr=c_G+c_A+c_T;
            ad=varG+varT;
            //int depth=w_A+c_A+c_T;
            int depth=totaldepth-w_T-c_A;
            double expectvar = (double) pow(0.1, ((qvalueG+qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varG+varT, totaldepth);

            if(totaldepth >= mincover  && qvalueG >= minquali && qvalueT>=minquali && varG >=minread2 && varT>=minread2 ){
                double C2G= (double)varG/totaldepth;
                double C2T= (double)varT/Ctotaldepth;
                if(C2G>=minhetfreq && C2T>=minhetfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, C2G, C2T, varG, varT, totaldepth, Ctotaldepth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, C2G, C2T, varG, varT, totaldepth, Ctotaldepth);
                }

            }else{
                double C2G= (double)varG/totaldepth;
                double C2T= (double)varT/Ctotaldepth;
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, C2G, C2T, varG, varT, totaldepth, Ctotaldepth);
            }

        }
        if(refbase == 'G'){//G>GT
            altbase[0]='T';altbase[1]='\0'; //1
            double qvalueT=(wsqT>crqT)?wsqT:crqT;
            int varT=c_T+w_T;
            adf=w_T;
            adr=c_T;
            ad=varT;
            //int depth=w_A+c_A+c_T;
            int depth=totaldepth;
            double expectvar = (double)pow(0.1, ((qvalueT)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varT, totaldepth);

            if(depth >= mincover  && qvalueT >= minquali  && varT >=minread2 ){
                double G2T= (double)varT/totaldepth;
                if(G2T>=minhetfreq ){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, G2T, 0, varT, 0, totaldepth, totaldepth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, G2T, 0, varT, 0, totaldepth, totaldepth);
                }

            }else{
                double G2T= (double)varT/totaldepth;
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, G2T, 0, varT, 0, totaldepth, totaldepth);
            }
        }
    }   
//CG
    if(genotypemaybe == "CG"){
        if(w_A+c_T > othercover) return 0;
        if(refbase == 'A'){//A>CG
            altbase[0]='C';altbase[1]=',';altbase[2]='G';altbase[3]='\0'; //2
            double qvalueG=(wsqG>crqG)?wsqG:crqG;
            double qvalueC=(wsqC>crqC)?wsqC:crqC;
            if(qvalueC < wsqT) qvalueC = wsqT;
            int varG=w_G+c_G;
            int varC=c_C+w_C+w_T;
            adf=w_G+w_C+w_T;
            adr=c_G+c_C;
            ad=varG+varC;
            //int depth=w_A+c_A+c_T;
            int depth=totaldepth-c_A;
            double expectvar = (double)pow(0.1, ((qvalueG+qvalueC)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varG+varC, totaldepth);

            if(totaldepth >= mincover  && qvalueG >= minquali && qvalueC>=minquali && varG >=minread2 && varC>=minread2 ){
                double A2G= (double)varG/depth;
                double A2C= (double)varC/totaldepth;
                if(A2G>=minhetfreq && A2C>=minhetfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, A2G, A2C, varG, varC, depth, totaldepth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, A2G, A2C, varG, varC, depth, totaldepth);
                }
            }else{
                double A2G= (double)varG/depth;
                double A2C= (double)varC/totaldepth;
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, A2G, A2C, varG, varC, depth, totaldepth);
            }

        }
        if(refbase == 'T'){//T>CG
            altbase[0]='C';altbase[1]=',';altbase[2]='G';altbase[3]='\0'; //2
            double qvalueG=(wsqG>crqG)?wsqG:crqG;
            if(qvalueG<crqA) qvalueG=crqA;
            double qvalueC=(wsqC>crqC)?wsqC:crqC;
            int varG=w_G+c_G;
            int varC=c_C+w_C;
            adf=w_G+w_C;
            adr=c_G+c_C;
            ad=varG+varC;
            //int depth=w_A+c_A+c_T;
            int depth=totaldepth-w_T-c_A;
            double expectvar = (double)pow(0.1, ((qvalueG+qvalueC)/10))/3;
            double pvalue = fishers_exact((int)expectvar*depth,depth,varG+varC, depth);

            if(totaldepth >= mincover  && qvalueG >= minquali && qvalueC>=minquali && varG >=minread2 && varC>=minread2 ){
                double T2G= (double)varG/depth;
                double T2C= (double)varC/depth;
                if(T2G>=minhetfreq && T2C>=minhetfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, T2C, T2G, varG, varC, depth, depth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, T2C, T2G, varG, varC, depth, depth);
                }

            }else{
                double T2G= (double)varG/depth;
                double T2C= (double)varC/depth;
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 2, pvalue, T2C, T2G, varG, varC, depth, depth);
            }
        }
        if(refbase == 'C'){//C>CG
            altbase[0]='G';altbase[1]='\0'; //1
            double qvalueG=(wsqG>crqG)?wsqG:crqG;
            if(qvalueG<crqA) qvalueG=crqA;
            int varG=w_G+c_G+c_A;
            adf=w_G;
            adr=c_G+c_A;
            ad=varG;
            //int depth=c_G+w_G+w_A+w_T+c_T+w_C+c_C;
            int depth=totaldepth-w_T-c_A;
            double C2G= (double)varG/totaldepth;
            double expectvar = (double)pow(0.1, ((qvalueG)/10))/3;
            double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varG, totaldepth);

            if(totaldepth >= mincover  && qvalueG >= minquali  && varG >=minread2 ){
                if(C2G>=minhetfreq){
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, C2G, 0, varG, 0, totaldepth, totaldepth);
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, C2G, 0, varG, 0, totaldepth, totaldepth);
                }

            }else{
                filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, C2G, 0, varG, 0, totaldepth, totaldepth);
            }           
        }
        if(refbase == 'G'){//G>CG
            altbase[0]='C';altbase[1]='\0'; //1
            double qvalueC=(wsqC>crqC)?wsqC:crqC;
            if(qvalueC<wsqT) qvalueC = wsqT;
                int varC=c_C+w_C+w_T;
                adf=w_C+w_T;
                adr=c_C;
                ad=varC;
                //int depth=c_G+w_G+w_A+w_T+c_T+w_C+c_C;
                int depth=totaldepth-w_T-c_A;
                double G2C= (double)varC/totaldepth;
                double expectvar = (double)pow(0.1, ((qvalueC)/10))/3;
                double pvalue = fishers_exact((int)expectvar*totaldepth,totaldepth,varC, totaldepth);

                if(totaldepth >= mincover  && qvalueC >= minquali  && varC >=minread2 ){
                    if(G2C>=minhetfreq){
                        filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 1 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, G2C, 0, varC, 0, totaldepth, totaldepth);
                    }else{
                        filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, G2C, 0, varC, 0, totaldepth, totaldepth);
                    }
                }else{
                    filpass = printSNP(posFptr, chrom, pos, refbase, altbase, genoqual, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, 0 , wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG,totaldepth, adf, adr, ad, 1, pvalue, G2C, 0, varC, 0, totaldepth, totaldepth);
                }
    
        }
    }   
    //unless(refbase=~/[ACGT]/i){
    return filpass;
        //print "lines[0]\tlines[1]\t\.\trefbase\t0\tSuper\tNN\t\.\t".join("\t",@lines[3..6])."\n";
    //}

}
