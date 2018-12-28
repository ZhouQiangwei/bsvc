/*
* Created on Fri Sep 21 10:29:50 2018
* @author: qwzhou
* @email: qiangwei.zhou2013@gmail.com
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <math.h>
#include "bayes.hpp"

double getmin(double a1,double a2,double a3,double a4, double c1,double c2,double c3,double c4){
    double min=1;
    if(a1<min) min=a1;
    if(a2<min) min=a2;
    if(a3<min) min=a3;
    if(a4<min) min=a4;
    if(c1<min) min=c1;
    if(c2<min) min=c2;
    if(c3<min) min=c3;
    if(c4<min) min=c4;
    return min;
}

double Factorial(int aa, int tt, int cc, int gg)
{
    int total=aa+tt+cc+gg;
    int naa=0;
    int ntt=0;int ncc=0;int ngg=0;int ntotal=0;
    if(aa<=1){
        naa=0; 
    }else{
        for(int i=1;i<=aa;i++){
            naa+=log(i);
        }
    }
    if(tt<=1){
                ntt=0; 
        }else{
            for(int i=1;i<=tt;i++){
                ntt+=log(i);
            }   
        }   
    if(cc<=1){
                ncc=0; 
        }else{
            for(int i=1;i<=cc;i++){
                ncc+=log(i);
            }   
        }   
    if(gg<=1){
        ngg=0; 
        }else{
            for(int i=1;i<=gg;i++){
                ngg+=log(i);
            }   
        }       
    if(total<=1){
        ntotal=0; 
        }else{
            for(int i=1;i<=total;i++){
                ntotal+=log(i);
            }    
        }   
    //print "\nntotal\t", ntotal, "\t", naa+ncc+ntt+ngg;
    double nn=ntotal-naa-ncc-ntt-ngg;
    return nn; 
}

void processgenotypesingle(double a, int& validN, double& first, double& second,double& fenmu, std::string& bestgeno, std::string& secondgeno, std::string genoNow){
    if(a!=0) validN++;
    else return;
    if(a>first) {
        second=first;
        secondgeno=bestgeno;
        first=a;
        bestgeno = genoNow;
    }else if(a>second) {
        second=a;
        secondgeno=genoNow;
    }
    fenmu+=pow(2.7, a);
}

void processGenotype(double AA,double AC,double AT,double AG,double CC,double CG,double CT,double GG,double GT,double TT,double aa,double ac,double at,double ag,double cc,double cg,double ct,double gg,double gt,double tt, int& validN,
    double& first, double& second, double & fenmu, std::string& bestgeno, std::string& secondgeno){
    if(AA!=0)
        processgenotypesingle(aa, validN, first, second, fenmu, bestgeno, secondgeno, "AA");
    if(AC!=0)
        processgenotypesingle(ac, validN, first, second, fenmu, bestgeno, secondgeno, "AC");
    if(AT!=0)
        processgenotypesingle(at, validN, first, second, fenmu, bestgeno, secondgeno, "AT");
    if(AG!=0)
        processgenotypesingle(ag, validN, first, second, fenmu, bestgeno, secondgeno, "AG");
    if(CC!=0)
        processgenotypesingle(cc, validN, first, second, fenmu, bestgeno, secondgeno, "CC");
    if(CG!=0)
        processgenotypesingle(cg, validN, first, second, fenmu, bestgeno, secondgeno, "CG");
    if(CT!=0)
        processgenotypesingle(ct, validN, first, second, fenmu, bestgeno, secondgeno, "CT");
    if(GG!=0)
        processgenotypesingle(gg, validN, first, second, fenmu, bestgeno, secondgeno, "GG");
    if(GT!=0)
        processgenotypesingle(gt, validN, first, second, fenmu, bestgeno, secondgeno, "GT");
    if(TT!=0)
        processgenotypesingle(tt, validN, first, second, fenmu, bestgeno, secondgeno, "TT");
}


void Bayes(int wsqA,int wsqT,int wsqC, int wsqG, int crqA, int crqT, int crqC, int crqG, char refbase, unsigned int pos, char* chrom,
    int w_A, int w_T, int w_C, int w_G, int c_A, int c_T, int c_C, int c_G, std::string& genotypemaybe, double &qual)
{
    int mincovercut=0;
    //print "line[0]\n";
    //double ptransition=0.00066;
    //double ptransversion=0.00033;
    //error of base
   //print "\n", pos, "\t" , refbase, "\n";
    //floor(a * 10000.000f + 0.5) / 10000.000f; /*保留小数点后四位*/
    double baseqWA=pow(0.1, (wsqA/10));
    double baseqCA=pow(0.1, (crqA/10));
    double baseqWT=pow(0.1, (wsqT/10));
    double baseqCT=pow(0.1, (crqT/10));
    double baseqWC=pow(0.1, (wsqC/10));
    double baseqCC=pow(0.1, (crqC/10));
    double baseqWG=pow(0.1, (wsqG/10));
    double baseqCG=pow(0.1, (crqG/10));
    double minQ=getmin(baseqWA,baseqCA,baseqWT, baseqCT, baseqWC, baseqCC, baseqWG, baseqCG);
    if(minQ==1 ){
        genotypemaybe="NN";
        qual=0;
        return;// "gntpmaybe\tqualerr"; 
    }


    //P(Di|g=Gj)
    double nn=Factorial(w_A,c_T,c_C,w_G);
    double aa=0;double ac=0;double at=0;double ag=0;double cc=0;double cg=0;double ct=0;double gg=0;double gt=0;double tt=0;
    int genomaycutoff=1; // 0

    //1、AA
    char ref=refbase;
    if(ref != 'A'){
        if(wsqA>0 && w_A>mincovercut && ref ==  'G'){
            double a_aa=1-baseqWA;
            double other=baseqWA/3;
            double nn=Factorial(w_A,c_T,c_C,w_G);

            aa = nn+ w_A*log(a_aa) + (c_T+w_T+w_C+c_C+w_G+c_G)*log(other) ;
        }
        if((wsqA>0 && w_A>mincovercut || crqA>0 && c_A>mincovercut) && ref != 'G'){
            double a_aa=1-baseqWA;
            double other=baseqWA/3;

            aa = nn+ (w_A + c_A)*log(a_aa) + (c_T+w_T+w_C+c_C+w_G+c_G)*log(other);
        }
    }
    //2、GG
    if(ref != 'G'){
        if((wsqG>0 && w_G>mincovercut|| crqG>0 && c_G>mincovercut || crqA>0 && c_A>mincovercut) && ref != 'A'){
            double gg_q;
            double other;
            if(wsqG>0 && w_G>mincovercut){
                gg_q = 1-baseqWG;
                other=baseqWG/3;
            }else if(crqG>0 && c_G>mincovercut){
                gg_q = 1-baseqCG;
                other=baseqCG/3;
            }else if(crqA>0 && c_A>mincovercut){
                gg_q = 1-baseqCA;
                other=baseqCA/3;
            }
            gg =  nn+ (w_G+c_G+c_A)*log(gg_q) + (w_A+c_T+w_T+c_C+w_C)*log(other);
        }
        if((wsqG>0 && w_G>mincovercut || crqG>0 && c_G>mincovercut) && ref ==  'A'){
            double gg_q;
            double other;
            if(wsqG>0 && w_G>mincovercut){
                gg_q = 1-baseqWG;
                other=baseqWG/3;
            }else if(crqG>0 && c_G>mincovercut){
                gg_q = 1-baseqCG;
                other=baseqCG/3;
            }
            gg = nn+ (w_G+c_G)*log(gg_q) + (w_A+c_T+w_T+c_C+w_C)*log(other);
        }
    }
    //3、TT
    if(ref != 'T'){
        if(crqT>0 && c_T>mincovercut && ref ==  'C'){
            double tt_q=1-baseqCT;
            double other=baseqCT/3;
            tt= nn+ c_T*log(tt_q)+ (w_A+c_C+w_G+c_G+w_C+c_A)*log(other);
        }
        if((wsqT>0 && w_T>mincovercut || crqT>0 && c_T>mincovercut) && ref != 'C'){
            double tt_q;
            double other;
            if(wsqT>0 && w_T>mincovercut){
                tt_q = 1-baseqWT;
                other=baseqWT/3;
            }else if(crqT>0 && c_T>mincovercut){
                tt_q = 1-baseqCT;
                other=baseqCT/3;
            }
            tt=  nn+ (w_T+c_T)*log(tt_q)+ (w_A+c_C+w_G+c_G+w_C+c_A)*log(other);
        }
    }
    //4、CC
    if(ref != 'C'){
        if((crqC>0 && c_C>mincovercut || wsqC>0 && w_C>mincovercut) && ref ==  'T'){
            double cc_q;
            double other;
            if(wsqC>0 && w_C>mincovercut){
                cc_q = 1-baseqWC;
                other=baseqWC/3;
            }else if(crqC>0 && c_C>mincovercut){
                cc_q = 1-baseqCC;
                other=baseqCC/3;
            }
            cc = nn+ (c_C+w_C)*log(cc_q) + (w_A+c_A+c_T+w_G+c_G)*log(other);
        }
        if((crqC>0 && c_C>mincovercut || wsqC>0 && w_C>mincovercut || wsqT>0 && w_T>mincovercut) && ref != 'T'){
            double cc_q;
            double other;
            if(wsqC>0 && w_C>mincovercut){
                cc_q = 1-baseqWC;
                other=baseqWC/3;
            }else if(crqC>0 && c_C>mincovercut){
                cc_q = 1-baseqCC;
                other=baseqCC/3;
            }else if(wsqT>0 && w_T>mincovercut){
                cc_q=1-baseqWT;
                other=baseqWT/3;
            }
            cc = nn+ (c_C+w_C+w_T)*log(cc_q) + (w_A+c_A+c_T+w_G+c_G)*log(other);
        }
    }
    //if(CC==TT) {//看C的数目。

    //5、AC
    if((wsqA>0 && w_A>mincovercut|| crqA>0 && c_A>mincovercut) &&ref != 'G'){
        double ac_q;
        double other;
        if(wsqA>0 && w_A>mincovercut){
            ac_q = (1-baseqWA)/2;
            other=baseqWA/3;
        }else if(crqA>0 && c_A>mincovercut){
            ac_q = (1-baseqCA)/2;
            other=baseqCA/3;
        }
        if((crqC>0 && c_C>mincovercut || wsqC>0 && w_C>mincovercut) && ref ==  'T'){
            ac=nn+ (c_C+w_C+w_A+c_A)*log(ac_q)+(c_T+w_G+c_G)*log(other);
        }
        if((crqC>0 && c_C>mincovercut || wsqC>0 && w_C>mincovercut || wsqT>0 && w_T>mincovercut) && ref != 'T'){
            ac=nn+ (c_C+w_C+w_T+w_A+c_A)*log(ac_q)+(c_T+w_G+c_G)*log(other);
        }
    }
    if(wsqA>0 && w_A>mincovercut && ref ==  'G'){
        double ac_q = (1-baseqWA)/2;
        double other=baseqWA/3;
        if((crqC>0 && c_C>mincovercut || wsqC>0 && w_C>mincovercut || wsqT>0 && w_T>mincovercut)){
            ac=nn+ (c_C+w_C+w_T+w_A)*log(ac_q)+(c_T+w_G+c_G)*log(other);
        }
    }
    //// 6、AG
    if( (wsqA>0 && w_A>mincovercut|| crqA>0 && c_A>mincovercut) && ref != 'G'){
        double ag_q;
        double other;
        if(wsqA>0 && w_A>mincovercut){
            ag_q = (1-baseqWA)/2;
            other=baseqWA/3;
        }else if(crqA>0 && c_A>mincovercut){
            ag_q = (1-baseqCA)/2;
            other=baseqCA/3;
        }
        if(crqG>0 && c_G>mincovercut || wsqG>0 && w_G>mincovercut){
            ag = nn+ (w_A +c_A+ w_G + c_G)*log(ag_q) + (c_C+c_T+w_C+ w_T)*log(other);
        }
    }
    if(wsqA>0 && w_A>mincovercut && ref ==  'G'){
        double ag_q = (1-baseqWA)/2;
        double other=baseqWA/3;
        if(crqG>0 && c_G>mincovercut || wsqG>0 && w_G>mincovercut){
            ag = nn+ (w_A + w_G + c_G)*log(ag_q) + (c_C+c_T+w_C+ w_T)*log(other);
        }
    }
     
    //7、AT
    if((wsqA>0 && w_A>mincovercut|| crqA>0 && c_A>mincovercut)&& ref != 'G'){
        double at_q;
        double other;
        if(wsqA>0 && w_A>mincovercut){
            at_q = (1-baseqWA)/2;
            other=baseqWA/3;
        }else if(crqA>0 && c_A>mincovercut){
            at_q = (1-baseqCA)/2;
            other=baseqCA/3;
        }
        if((wsqT >0 || crqT>0 && c_T>mincovercut)&& ref != 'C'){
            at = nn+ (w_A+c_A+w_T+c_T)*log(at_q) + (c_C+w_G+w_C+c_G)*log(other);
        }
        if(crqT>0 && c_T>mincovercut && ref ==  'C'){
            at = nn+ (w_A+c_A+c_T)*log(at_q) + (c_C+w_G+w_C+c_G)*log(other);
        }
    }
    if(wsqA>0 && w_A>mincovercut&& ref ==  'G'){
        double at_q = (1-baseqWA)/2;
        double other=baseqWA/3;
        if(wsqT>0 && w_T>mincovercut || crqT>0 && c_T>mincovercut){
            at = nn+ (w_A+w_T+c_T)*log(at_q) + (c_C+w_G+w_C+c_G)*log(other);
        }
    }
     
    //8、TC
    if(crqT>0 && c_T>mincovercut && ref ==  'C'){
        double tc_q = (1-baseqCT)/2;
        double other=baseqCT/3;
        if(wsqC>0 && w_C>mincovercut || crqC>0 && c_C>mincovercut){
            ct=  nn+ (w_C+c_T+c_C)*log(tc_q) + (w_A+w_G+c_G+c_A)*log(other);
        }
    }
    if((crqT>0 && c_T>mincovercut||wsqT>0 && w_T>mincovercut) && ref != 'C'){
        double tc_q;
        double other;
        if(crqT>0 && c_T>mincovercut){
            tc_q = (1-baseqCT)/2;
            other=baseqCT/3;
        }else if(wsqT>0 && w_T>mincovercut){
            tc_q = (1-baseqWT)/2;
            other=baseqWT/3;
        }
        if(wsqC>0 && w_C>mincovercut || crqC>0 && c_C>mincovercut){
            ct=  nn+ (w_C+w_T+c_T+c_C)*log(tc_q) + (w_A+w_G+c_G+c_A)*log(other);
        }
    }
     
    //9、TG
    if((crqT>0 && c_T>mincovercut||wsqT>0 && w_T>mincovercut) && ref != 'C' ){
        double tg_q;
        double other;
        if(crqT>0 && c_T>mincovercut){
            tg_q = (1-baseqCT)/2;
            other=baseqCT/3;
        }else if(wsqT>0 && w_T>mincovercut){
            tg_q = (1-baseqWT)/2;
            other=baseqWT/3;
        }
        if((wsqG>0 && w_G>mincovercut || crqG>0 && c_G>mincovercut|| crqA>0 && c_A>mincovercut) && ref != 'A'){
            gt= nn+ (c_T+w_T+w_G+c_G+c_A)*log(tg_q) + (w_A+c_C+w_C)*log(other);        
        }
        if((wsqG>0 && w_G>mincovercut || crqG>0 && c_G>mincovercut) && ref ==  'A'){
            gt= nn+ (c_T+w_T+w_G+c_G)*log(tg_q) + (w_A+c_C+w_C)*log(other);        
        }
    }
    if(crqT>0 && c_T>mincovercut && ref ==  'C'){
        double tg_q = (1-baseqCT)/2;
        double other=baseqCT/3;
        if(wsqG>0 && w_G>mincovercut || crqG>0 && c_G>mincovercut||crqA>0 && c_A>mincovercut){
            gt=  nn+ (c_T+w_G+c_G+c_A)*log(tg_q) + (w_A+c_C+w_C)*log(other);       
        }
    }
     
    //10、CG
    if((crqC>0 && c_C>mincovercut || wsqC>0 && w_C>mincovercut) && ref ==  'T' ){
        double cg_q;
        double other;
        if(crqC>0 && c_C>mincovercut){
            cg_q = (1-baseqCC)/2;
            other=baseqCC/3;
        }else if(wsqC>0 && w_C>mincovercut){
            cg_q = (1-baseqWC)/2;
            other=baseqWC/3;
        }
        cg= nn+ (c_C+w_G+w_C+c_G+c_A)*log(cg_q) + (w_A+c_T)*log(other);  
    }
    if((crqC>0 && c_C>mincovercut || wsqC>0 && w_C>mincovercut || wsqT>0 && w_T>mincovercut) && ref != 'T'){
        double cg_q;
        double other;
        if(crqC>0 && c_C>mincovercut){
            cg_q = (1-baseqCC)/2;
            other=baseqCC/3;
        }else if(wsqC>0 && w_C>mincovercut){
            cg_q = (1-baseqWC)/2;
            other=baseqWC/3;
        }else if(wsqT>0 && w_T>mincovercut){
            cg_q = (1-baseqWT)/2;
            other=baseqWT/3;
        }
        if((w_G>0 || c_G>0|| w_T>0) && ref ==  'A'){
            cg= nn+ (c_C+w_G+w_C+c_G+w_T)*log(cg_q) + (w_A+c_T)*log(other);  
        }
        if((wsqG>0 && w_G>mincovercut || crqG>0 && c_G>mincovercut||crqA>0 && c_A>mincovercut) && ref != 'A'){
            cg= nn+ (c_C+w_G+w_C+c_G+ w_T+c_A)*log(cg_q) + (w_A+c_T)*log(other);  
        }
    }

    
    

    //P(D|g=Gj)                  
    //sum(P(Gj)*P(D|g=Gj))
    double fenmu=0;
    double AA=aa;double AC=ac;double AT=at;double AG=ag;double CC=cc;double CG=cg;double CT=ct;double GG=gg;double GT=gt;double TT=tt;
    
    if(refbase ==  'A'){
        aa+=log(0.985);
        tt+=log(0.000083);
        cc+=log(0.000083);
        gg+=log(0.00033);
        at+=log(0.00017);
        ac+=log(0.00017);
        ag+=log(0.000667);
        ct+=(log(2.78) - 8*log(10));
        gt+=(log(1.1) - 7*log(10));
        cg+=(log(1.1) - 7*log(10));    
    }    
    if(refbase ==  'T'){
        aa+=log(0.000083);
        tt+=log(0.985);
        cc+=log(0.00033);
        gg+=log(0.000083);
        at+=log(0.00017);
        ac+=(log(1.1) - 7*log(10));
        ag+=(log(2.78) - 8*log(10));
        ct+=log(0.000667);
        gt+=log(0.00017);
        cg+=(log(1.1) - 7*log(10));
    }    
    if(refbase ==  'C'){
        aa+=log(0.000083);
        tt+=log(0.00033);
        cc+=log(0.985);
        gg+=log(0.000083);
        at+=(log(1.1) - 7*log(10));
        ac+=log(0.00017);
        ag+=(log(2.78) - 8*log(10));
        ct+=log(0.000667);
        gt+=(log(1.1) - 7*log(10));
        cg+=log(0.00017);
    }       
    if(refbase ==  'G'){
        aa+=log(0.00033);
        tt+=log(0.000083);
        cc+=log(0.000083);
        gg+=log(0.9985);
        at+=(log(1.1) - 7*log(10));
        ac+=(log(1.1) - 7*log(10));
        ag+=(log(6.67) - 4*log(10));
        ct+=(log(2.78) - 8*log(10));
        gt+=(log(1.67) - 4*log(10));
        cg+=(log(1.67) - 4*log(10));
    }

    int validN=0;
    double bestvalue=-10000;double secondvalue=-10000;
    std::string bestgeno="NN"; std::string secondgeno="NN";
    processGenotype(AA, AC,AT,AG,CC,CG,CT,GG,GT,TT, aa,ac,at,ag,cc,cg,ct,gg,gt,tt, validN, bestvalue, secondvalue, fenmu, bestgeno, secondgeno);

    //version 0.42
    if( (w_C+c_C<1 || w_G+c_G<1) && bestvalue==secondvalue){
        if(bestgeno=="CC" && secondgeno=="TT"){
            if(w_C+c_C<1) bestgeno=secondgeno;
        }else if(bestgeno=="GG" && secondgeno=="AA"){
            if(w_G+c_G<1) bestgeno=secondgeno;
        }else if(bestgeno=="CG" && secondgeno=="GT"){
            if(w_C+c_C<1) bestgeno=secondgeno;
        }else if(bestgeno=="CG" && secondgeno=="AC"){
            if(w_G+c_G<1) bestgeno=secondgeno;
        }else if(bestgeno=="GT" && secondgeno=="AT"){
            if(w_G+c_G<1) bestgeno=secondgeno;
        }else if(bestgeno=="AC" && secondgeno=="AT"){
            if(w_C+c_C<1) bestgeno=secondgeno;
        }
    }
    //printf("\n%s %s %f %f %f %d\n", bestgeno.c_str(), secondgeno.c_str(), bestvalue, secondvalue, fenmu, validN);
    //if(validN>1) {
    //    print "\nsecond\t", sort[1], "\t", hash2{sort[1]}. "\n";
    //}
    double prob=0;
    if(validN==0){
        genotypemaybe="NN";
        qual=0;
    }else{
        genotypemaybe=bestgeno;
        double first=pow(2.7, bestvalue);
        if(validN>1){
            if(fenmu==0){
                qual=1000;
            }else{
                prob=1-first/fenmu;
                if(prob==0){

                    qual=1000;
                }else{
                    qual=-10*log(prob)/log(10);
                }
            }
        }else if(validN==1){
            if(bestgeno ==  "AA"){
                double hom=w_A; 
                prob=1-1/(1+pow(0.5,hom));
            }
            else if(bestgeno ==  "TT"){
                double hom=c_T;     
                prob=1-1/(1+pow(0.5,hom));
            }   
            else if(bestgeno ==  "CC"){
                double hom=w_C+c_C;
                prob=1-1/(1+pow(0.5,hom));
            }
            else if(bestgeno ==  "TT"){
                double hom=w_G+c_G;
                prob=1-1/(1+pow(0.5,hom));
            }else{
                prob=1;
            }

            if(prob==0){
                qual=1000; 
            }else{
                qual=-10*log(prob)/log(10);
            }
        }
    }
    
    qual=int(qual);

    return;
    //return genotypemaybe;
    //print "bestgeno\thash{bestgeno}\tsort[1]\thash{sort[1]}\n";
}
