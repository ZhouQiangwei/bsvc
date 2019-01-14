/*
* Created on Fri Sep 21 10:29:50 2018
* @author: qwzhou
* @email: qiangwei.zhou2013@gmail.com
* partail refered form bssnper code
*/
#include "bamprocess.h"
#include <string.h>
#include <string>
#include "bayes.hpp"
#include "genotype.hpp"

extern unsigned int longestchr;
extern int minquali;
extern int mincover;
extern int maxcover;
extern float minhetfreq;
extern float errorrate;
extern unsigned int mapqThr;
extern int meth;
extern int snp;
extern int NTHREAD;
extern char bamFileName[1024];
extern char** chrName;
extern int minvarread;
extern int methmincover;
extern int multiout;

int ithreadschr = 0;
pthread_mutex_t meth_counter_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t output_mutex = PTHREAD_MUTEX_INITIALIZER;

void *npsnpAnalysis(void *arg){
    ARGS *args = (struct ARGS*)arg;
    //if(meth)
    //    methProcess(bamFileName, args->methCgFileName, args->methChgFileName, args->methChhFileName, args->hashTable, args->chrSeqArray, args->chrLen, args->chrCnt, minquali, mincover, mapqThr);

    //Process SNP
    if(snp){
        // SNP outfile
        if(!multiout && args->snpFptr == NULL) {
            fprintf(stderr, "Could not open SNP output file!\n");
            exit(1);
        }
        char tempoutfile[1000];
        if(NTHREAD>1){
            int processchrom=0;
            for(; ithreadschr < args->chrCnt;){
                pthread_mutex_lock(&meth_counter_mutex);
                processchrom=ithreadschr;
                ithreadschr++;
                pthread_mutex_unlock(&meth_counter_mutex);
                // SNP outfile
                if(multiout){
                    sprintf(tempoutfile, "%s.%s.vcf", args->snpFileName, chrName[processchrom]);
                    FILE* snptempfp = fopen(tempoutfile, "w");
                    sprintf(tempoutfile, "%s.%s.vcf", args->methFileName, chrName[processchrom]);
                    FILE* methtempfp;
		    if(meth==1) methtempfp = fopen(tempoutfile, "w");
                    snpProcess_multiop(methtempfp, snptempfp, bamFileName, args->hashTable, args->chrSeqArray, args->chrLen, args->chrCnt, minquali, maxcover, minhetfreq, errorrate, mapqThr, chrName[processchrom]);
                    fclose(snptempfp);
                    fclose(methtempfp);
                }else
                    snpProcess(args->methFptr, args->snpFptr, bamFileName, args->hashTable, args->chrSeqArray, args->chrLen, args->chrCnt, minquali, maxcover, minhetfreq, errorrate, mapqThr, chrName[processchrom]);
            }
        }
        else{
            snpProcess_singlet(args->methFptr, args->snpFptr, bamFileName, args->hashTable, args->chrSeqArray, args->chrLen, args->chrCnt, minquali, maxcover, minhetfreq, errorrate, mapqThr);
        }
        
    }
}

void printMeth(FILE* methFptr, int len, char* curChr, unsigned int* w_Mm_CG, unsigned int* w_Mc_CG, unsigned int* w_Mq_CG, unsigned int* c_Mm_CG, unsigned int* c_Mc_CG, unsigned int* c_Mq_CG, char* tag, int mincover)
{
    int i;

    // Record meth sites
    if(strcmp(tag, "CG") == 0) {
        for(i = 0; i < len - 1; i++) {
            if(w_Mc_CG[i] + c_Mc_CG[i+1] >= mincover)
                if(w_Mc_CG[i] != 0 && c_Mc_CG[i+1] != 0)
                    fprintf(methFptr, "%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", curChr, i + 1, tag, w_Mm_CG[i], w_Mc_CG[i], (unsigned int)(w_Mq_CG[i]/w_Mc_CG[i]), c_Mm_CG[i+1], c_Mc_CG[i+1], (unsigned int)(c_Mq_CG[i+1]/c_Mc_CG[i+1]));
                else if(w_Mc_CG[i] == 0)
                    fprintf(methFptr, "%s\t%d\t%s\t.\t.\t.\t%d\t%d\t%d\n", curChr, i + 1, tag, c_Mm_CG[i+1], c_Mc_CG[i+1], (unsigned int)(c_Mq_CG[i+1]/c_Mc_CG[i+1]));
                else
                    fprintf(methFptr, "%s\t%d\t%s\t%d\t%d\t%d\t.\t.\t.\n", curChr, i + 1, tag, w_Mm_CG[i], w_Mc_CG[i], (unsigned int)(w_Mq_CG[i]/w_Mc_CG[i]));
        }
    }
    else if(strcmp(tag, "CHG")){
        for(i = 0; i < len - 2; i++) {
            if(w_Mc_CG[i] + c_Mc_CG[i+2] >= mincover)
                if(w_Mc_CG[i] != 0 && c_Mc_CG[i+2] != 0)
                    fprintf(methFptr, "%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", curChr, i + 1, tag, w_Mm_CG[i], w_Mc_CG[i], (unsigned int)(w_Mq_CG[i]/w_Mc_CG[i]), c_Mm_CG[i+2], c_Mc_CG[i+2], (unsigned int)(c_Mq_CG[i+2]/c_Mc_CG[i+2]));
                else if(w_Mc_CG[i] == 0)
                    fprintf(methFptr, "%s\t%d\t%s\t.\t.\t.\t%d\t%d\t%d\n", curChr, i + 3, tag, c_Mm_CG[i+2], c_Mc_CG[i+2], (unsigned int)(c_Mq_CG[i+2]/c_Mc_CG[i+2]));
                else
                    fprintf(methFptr, "%s\t%d\t%s\t%d\t%d\t%d\t.\t.\t.\n", curChr, i + 1, tag, w_Mm_CG[i], w_Mc_CG[i], (unsigned int)(w_Mq_CG[i]/w_Mc_CG[i]));
        }
    }
    else {
        for(i = 0; i < len - 2; i++) {
            if(w_Mc_CG[i] + c_Mc_CG[i+2] >= mincover) {
                if(w_Mc_CG[i] > 0)
                    fprintf(methFptr, "%s\t%d\t%s\t%d\t%d\t%d\t.\t.\t.\n", curChr, i + 1, tag, w_Mm_CG[i], w_Mc_CG[i], (unsigned int)(w_Mq_CG[i]/w_Mc_CG[i]));
                if(c_Mc_CG[i+2] > 0)
                    fprintf(methFptr, "%s\t%d\t%s\t.\t.\t.\t%d\t%d\t%d\n", curChr, i + 3, tag, c_Mm_CG[i+2], c_Mc_CG[i+2], (unsigned int)(c_Mq_CG[i+2]/c_Mc_CG[i+2]));
            }
        }

    }
}

void printSnp(FILE* methFptr, FILE* posFptr, char** chrSeqArray, int idx, int len, float minhetfreq, char* curChr, unsigned short *w_A, unsigned short *w_T, unsigned short *w_C, unsigned short *w_G, unsigned short *c_A, unsigned short *c_T, unsigned short *c_C, unsigned short *c_G, unsigned short *w_Aq, unsigned short *w_Tq, unsigned short *w_Cq, unsigned short *w_Gq, unsigned short *c_Aq, unsigned short *c_Tq, unsigned short *c_Cq, unsigned short *c_Gq, 
    unsigned short *w_Q, unsigned short *c_Q)
{
    int i, j, m,ccover,wcover;
    int v, n, v2, n2;
    float vRate[3];

    //meth v
    char context[5];
    float mratio=0;
    int ct,ga;
    // Record snp sites
    for(i = 0; i < len; i++)
    {

        m = w_A[i]+w_T[i]+w_C[i]+w_G[i]+c_A[i]+c_T[i]+c_C[i]+c_G[i];
        ccover=c_A[i]+c_T[i]+c_C[i]+c_G[i];
        wcover=w_A[i]+w_T[i]+w_C[i]+w_G[i];
        if(m<1) continue;
        // No snp pt number
        switch(chrSeqArray[idx][i])
        {
        case 'A':
            // A>C=(C+WsT)/total reads
            v = w_C[i] + c_C[i] + w_T[i];
            n = m;
            vRate[0] = (float)v/n;
            // A>T=T/total reads
            v = w_T[i] + c_T[i];
            n = m;
            vRate[1] = (float)v/n;
            // A>G=WsG/(WsG+WsA) ==> G/m
            v = w_G[i];
            n = w_G[i] + w_A[i];
            v2 = w_G[i] + c_G[i];
            n2 = m;
            vRate[2] = (float)v/n>(float)v2/n2?(float)v/n:(float)v2/n2;
            break;
        case 'T':
            // T>A=A/total reads
            v = w_A[i] + c_A[i];
            n = m;
            vRate[0] = (float)v/n;
            // T>C=CrC/(CrT+CrC) ==> C/(C+CrT)
            v = c_C[i];
            n = c_T[i] + c_C[i];
            v2 = c_C[i] + w_C[i];
            n2 = m;
            vRate[1] = (float)v/n>(float)v2/n2?(float)v/n:(float)v2/n2;
            // T>G=(G+CrA)/total reads
            v = w_G[i] + c_G[i] + c_A[i];
            n = m;
            vRate[2] = (float)v/n;
            break;
        case 'C':
            // C>A=A/total reads
            v = w_A[i] + c_A[i];
            n = m;
            vRate[0] = (float)v/n;
            // C>T=CrT/(CrC+CrT)
            v = c_T[i];
            n = c_C[i] + c_T[i];
            vRate[1] = (float)v/n;
            // C>G=(G+CrA)/total reads
            v = w_G[i] + c_G[i] + c_A[i];
            n = m;
            vRate[2] = (float)v/n;
            break;
        case 'G':
            // G>A=WsA/(WsA+WsG)
            v = w_A[i];
            n = w_A[i] + w_G[i];
            vRate[0] = (float)v/n;
            // G>T=T/total reads
            v = w_T[i] + c_T[i];
            n = m;
            vRate[1] = (float)v/n;
            // G>C=(WsT+C)/total reads
            v = w_T[i] + w_C[i] + c_C[i];
            n = m;
            vRate[2] = (float)v/n;
            break;
        }
        int filterpass=0;
        // Filtering
        char refbase=chrSeqArray[idx][i];
        unsigned int wsqA= (unsigned int)((float)w_Aq[i]/w_A[i]+0.5);
        unsigned int wsqT= (unsigned int)((float)w_Tq[i]/w_T[i]+0.5);
        unsigned int wsqC= (unsigned int)((float)w_Cq[i]/w_C[i]+0.5);
        unsigned int wsqG= (unsigned int)((float)w_Gq[i]/w_G[i]+0.5);
        unsigned int crqA= (unsigned int)((float)c_Aq[i]/c_A[i]+0.5);
        unsigned int crqT=  (unsigned int)((float)c_Tq[i]/c_T[i]+0.5);
        unsigned int crqC= (unsigned int)((float)c_Cq[i]/c_C[i]+0.5);
        unsigned int crqG=  (unsigned int)((float)c_Gq[i]/c_G[i]+0.5);
        std::string genotypemaybe="NN";
        if(vRate[0] > minhetfreq || vRate[1] > minhetfreq || vRate[2] > minhetfreq )
        {
            double qual=1;
            Bayes(wsqA,wsqT,wsqC, wsqG, crqA, crqT, crqC, crqG, refbase, i+1, curChr, w_A[i], w_T[i], w_C[i], w_G[i], c_A[i], c_T[i], c_C[i], c_G[i], genotypemaybe, qual);
            filterpass = genotype(posFptr, wsqA, wsqT, wsqC,  wsqG,  crqA,  crqT,  crqC,  crqG,  refbase, i+1, curChr,
     w_A[i],  w_T[i],  w_C[i],  w_G[i],  c_A[i], c_T[i], c_C[i], c_G[i], genotypemaybe, qual);
        }

        //Meth
        if(meth && m>=methmincover && m<maxcover){
            context[0]=refbase;context[1]='\0';
            //chrM  1   -   CHH 0   434 0.000000    Cq,Tq refbase,genotype
            if(filterpass==0){
                
                if(refbase=='C'){
                    ct=w_C[i]+w_T[i];
                    if(ct>methmincover){
                        if(i+1 < len && chrSeqArray[idx][i+1]=='G'){
                            context[0]='C';
                            context[1]='G';
                            context[2]='\0';
                        }else if(i+2 < len && chrSeqArray[idx][i+2]=='G'){
                            context[0]='C';
                            context[1]='H';
                            context[2]='G';
                            context[3]='\0';
                        }else if(i+2 < len){
                            context[0]='C';
                            context[1]='H';
                            context[2]='H';
                            context[3]='\0';
                        }
                        fprintf(methFptr, "%s\t%d\t+\t%s\t%d\t%d\t%f\t%d,%d\t%d,%d\t%c\n",curChr, i+1, context, w_C[i], ct, (float)w_C[i]/(ct), wsqC,wsqT,wcover, ccover, refbase);
                    }
                }
                else if(refbase=='G'){
                    ga=c_G[i]+c_A[i];
                    if(ga>methmincover){
                        if(i > 0 && chrSeqArray[idx][i-1]=='C'){
                            context[0]='C';
                            context[1]='G';
                            context[2]='\0';
                        }else if(i > 1 && chrSeqArray[idx][i-2]=='C'){
                            context[0]='C';
                            context[1]='H';
                            context[2]='G';
                            context[3]='\0';
                        }else if(i > 1){
                            context[0]='C';
                            context[1]='H';
                            context[2]='H';
                            context[3]='\0';
                        }
                        fprintf(methFptr, "%s\t%d\t-\t%s\t%d\t%d\t%f\t%d,%d\t%d,%d\t%c\n",curChr, i+1, context, c_G[i], ga, (float)c_G[i]/ga, crqG,crqA,wcover, ccover,refbase);
                    }
                }
            }else if(filterpass==1 && genotypemaybe!="CT" && genotypemaybe!="AG"){
                
                if(genotypemaybe=="CG"){
                    ct=w_C[i]+w_T[i];
                    if(ct>methmincover){
                        if(i+1 < len && chrSeqArray[idx][i+1]=='G'){
                            context[0]='C';
                            context[1]='G';
                            context[2]='\0';
                        }else if(i+2 < len && chrSeqArray[idx][i+2]=='G'){
                            context[0]='C';
                            context[1]='H';
                            context[2]='G';
                            context[3]='\0';
                        }else if(i+2 < len){
                            context[0]='C';
                            context[1]='H';
                            context[2]='H';
                            context[3]='\0';
                        }
                        fprintf(methFptr, "%s\t%d\t+\t%s\t%d\t%d\t%f\t%d,%d\t%d,%d\t%c,%s\n",curChr, i+1, context, w_C[i], ct, (float)w_C[i]/ct, wsqC,wsqT,wcover, ccover,refbase,genotypemaybe.c_str());
                    }
                    ga=c_G[i]+c_A[i];
                    if(ga>methmincover){
                        if(i > 0 && chrSeqArray[idx][i-1]=='C'){
                            context[0]='C';
                            context[1]='G';
                            context[2]='\0';
                        }else if(i > 1 && chrSeqArray[idx][i-2]=='C'){
                            context[0]='C';
                            context[1]='H';
                            context[2]='G';
                            context[3]='\0';
                        }else if(i > 1){
                            context[0]='C';
                            context[1]='H';
                            context[2]='H';
                            context[3]='\0';
                        }
                        fprintf(methFptr, "%s\t%d\t-\t%s\t%d\t%d\t%f\t%d,%d\t%d,%d\t%c,%s\n",curChr, i+1, context, c_G[i], ga, (float)c_G[i]/ga, crqG,crqA,wcover, ccover,refbase,genotypemaybe.c_str());
                    }
                }else if(genotypemaybe[0]=='C' || genotypemaybe[1] == 'C'){
                    ct=w_C[i]+w_T[i];
                    if(ct>methmincover){
                        if(i+1 < len && chrSeqArray[idx][i+1]=='G'){
                            context[0]='C';
                            context[1]='G';
                            context[2]='\0';
                        }else if(i+2 < len && chrSeqArray[idx][i+2]=='G'){
                            context[0]='C';
                            context[1]='H';
                            context[2]='G';
                            context[3]='\0';
                        }else if(i+2 < len){
                            context[0]='C';
                            context[1]='H';
                            context[2]='H';
                            context[3]='\0';
                        }
                        fprintf(methFptr, "%s\t%d\t+\t%s\t%d\t%d\t%f\t%d,%d\t%d,%d\t%c,%s\n",curChr, i+1, context, w_C[i], ct, (float)w_C[i]/ct, wsqC,wsqT,wcover, ccover,refbase,genotypemaybe.c_str());
                    }
                }else if(genotypemaybe[0]=='G' || genotypemaybe[1] == 'G'){
                    ga=c_G[i]+c_A[i];
                    if(ga>methmincover){
                        if(i > 0 && chrSeqArray[idx][i-1]=='C'){
                            context[0]='C';
                            context[1]='G';
                            context[2]='\0';
                        }else if(i > 1 && chrSeqArray[idx][i-2]=='C'){
                            context[0]='C';
                            context[1]='H';
                            context[2]='G';
                            context[3]='\0';
                        }else if(i > 1){
                            context[0]='C';
                            context[1]='H';
                            context[2]='H';
                            context[3]='\0';
                        }
                        fprintf(methFptr, "%s\t%d\t-\t%s\t%d\t%d\t%f\t%d,%d\t%d,%d\t%c,%s\n",curChr, i+1, context, c_G[i], ga, (float)c_G[i]/ga, crqG,crqA,wcover, ccover,refbase,genotypemaybe.c_str());
                    }
                }

            }
        }

    }
}

int parseBuffer(bam_header_t *header, bam1_t *b, MapRecord* record, unsigned int mapqThr)
{
    int i, num, plen, flag;
    int seqPtr, seqBufPtr;
    char nchar;
    char *pdest;

    const bam1_core_t *c = &b->core;
    uint8_t *s = bam1_seq(b);
    uint8_t *t = bam1_qual(b);

    //kstring_t str;
    //str.l = str.m = 0; str.s = 0;

    /////////////////////////////////
    // Load record
    /////////////////////////////////
    // {{
    // 1 - QNAME
    strcpy(record->qname, bam1_qname(b));
    // 2 - FLAG
    record->flag = c->flag;
    if((record->flag & 0x04) != 0) {
        // Unmap
        return 1;
    }
    if((record->flag & 0x10) == 0) {
        record->strand = '+';
    }
    else {
        record->strand = '-';
    }
    if((record->flag & 0x40) != 0) {
        record->r12 = 1;
    }
    else {
        record->r12 = 2;
    }
    // 3 - RNAME
    strcpy(record->chrome, header->target_name[c->tid]);
    // 4 - POS
    record->offset = c->pos + 1;
    // 5 - MAPQ
    record->mapq = c->qual;
    if(c->qual < mapqThr)
        return 1;
    // 6 - CIGAR
    //char* tmp=(char*)calloc(1000, sizeof(char));
    char tmp[256];
    int j=0;

    //int cl=0;
    for (i = 0; i < c->n_cigar; ++i){
        sprintf(tmp, "%d%c", bam1_cigar(b)[i]>>BAM_CIGAR_SHIFT, "MIDNSHP"[bam1_cigar(b)[i]&BAM_CIGAR_MASK]);
        //cl=strlen(tmp)-1;
        //if(strstr(tmp, "M")!=NULL || strstr(tmp, "I")!=NULL || strstr(tmp, "D")!=NULL || strstr(tmp, "S")!=NULL) {
        //if(cl>0 && (tmp[cl]=='M' || tmp[cl]=='I' || tmp[cl]=='D' || tmp[cl]=='S')){
        for(int l= 0; l<strlen(tmp); ++l,j++) record->cigar[j] = tmp[l];
        //}else return 1;
        //if(i==0)
        //    strcpy(record->cigar, tmp);
        //else
        //    strcat(record->cigar, tmp);
    }

    record->cigar[j]= '\0';
    //record->cigar = str.s;
    
    //free(tmp);
    //tmp=NULL;
    // 10 - SEQ
    for(i = 0; i < c->l_qseq; ++i)
        record->seqBuf[i] = bam_nt16_rev_table[bam1_seqi(s, i)];
    record->seqBuf[i] = '\0';
    // 11 - QUAL
    if(t[0] == 0xff) {
        strcpy(record->qualBuf, "*");
    }
    else {
        for(i = 0; i < c->l_qseq; ++i)
            record->qualBuf[i] = (char)(t[i] + 33);
    }
    // }}

    /////////////////////////////////
    // To upper
    /////////////////////////////////
    // {{
    for(i = 0; i < strlen(record->seqBuf); i++)
        record->seqBuf[i] = toupper(record->seqBuf[i]);
    // }}

    /////////////////////////////////
    // Cigar analysis
    /////////////////////////////////
    // {{
    pdest = record->cigar;
    nchar = nxtChar(pdest, &plen);
    seqPtr = 0;
    seqBufPtr = 0;
    while(nchar != 'Z') {
        // Digit parsing
        for(i = 0; i < plen; i++)
            record->comBuf[i] = pdest[i];
        record->comBuf[plen] = '\0';
        num = atoi(record->comBuf);
        // Sequence handling
        switch(nchar) {
            case 'M':
                for(i = 0; i < num; i++) {
                    record->seq[seqPtr] = record->seqBuf[seqBufPtr];
                    record->qual[seqPtr] = record->qualBuf[seqBufPtr];
                    seqPtr++;
                    seqBufPtr++;
                }
                break;
            case 'D':
                for(i = 0; i < num; i++) {
                    record->seq[seqPtr] = 'N';
                    record->qual[seqPtr] = '5';
                    seqPtr++;
                }
                break;
            case 'I':
                seqBufPtr = seqBufPtr + num;
                break;
            case 'S':
                for(i = 0; i < num; i++) {
                    //record->seq[seqPtr] = 'N';
                    //record->qual[seqPtr] = record->qualBuf[seqBufPtr];
                    //seqPtr++;
                    seqBufPtr++;
                }
                break;
            case 'H':
                for(i = 0; i < num; i++) {
                    //record->seq[seqPtr] = 'N';
                    //record->qual[seqPtr] = record->qualBuf[seqBufPtr];
                    //seqPtr++;
                    seqBufPtr++;
                }
                break;
            default:
                fprintf(stderr, "Error in cigar parsing!\n");
		break;
                exit(1);
        }

        pdest = pdest + plen + 1;
        nchar = nxtChar(pdest, &plen);
    }
    record->len = seqPtr;
    record->seq[seqPtr] = '\0';
    record->qual[seqPtr] = '\0';
    // }}
    return 0;
}

char nxtChar(char* cigar, int* len)
{
    // "len" is the length of digit part
    int tLen, mLen;
    char mChar;
    char* pdest;

    mLen = 1000;
    mChar = 'Z';

    pdest = strstr(cigar, "M");
    if(pdest != NULL) {
        tLen = (int)(pdest - cigar);
        if(tLen < mLen) {
            mLen = tLen;
            mChar = 'M';
        }
    }

    pdest = strstr(cigar, "D");
    if(pdest != NULL) {
        tLen = (int)(pdest - cigar);
        if(tLen < mLen) {
            mLen = tLen;
            mChar = 'D';
        }
    }

    pdest = strstr(cigar, "I");
    if(pdest != NULL) {
        tLen = (int)(pdest - cigar);
        if(tLen < mLen) {
            mLen = tLen;
            mChar = 'I';
        }
    }

    pdest = strstr(cigar, "S");
    if(pdest != NULL) {
        tLen = (int)(pdest - cigar);
        if(tLen < mLen) {
            mLen = tLen;
            mChar = 'S';
        }
    }

    pdest = strstr(cigar, "H");
    if(pdest != NULL) {
        tLen = (int)(pdest - cigar);
        if(tLen < mLen) {
            mLen = tLen;
            mChar = 'H';
        }
    }

    if(mLen < 1000)
        *len = mLen;

    return mChar;
}

void dispRecord(MapRecord* record)
{
    int i;

    for(i = 0; i < 40; i++)
        fprintf(stderr, "*");
    fprintf(stderr, "\n");

    fprintf(stderr, "qname: %s\n", record->qname);
    fprintf(stderr, "strand: %c\n", record->strand);
    fprintf(stderr, "chrome: %s\n", record->chrome);
    fprintf(stderr, "flag: %d\n", record->flag);
    fprintf(stderr, "offset: %d\n", record->offset);
    fprintf(stderr, "mapq: %d\n", record->mapq);
    fprintf(stderr, "cigar: %s\n", record->cigar);
    fprintf(stderr, "seq: %s\n", record->seq);
    fprintf(stderr, "qual: %s\n", record->qual);
    fprintf(stderr, "len: %d\n", record->len);
    fprintf(stderr, "r12: %d\n", record->r12);

    for(i = 0; i < 40; i++)
        fprintf(stderr, "*");
    fprintf(stderr, "\n");
}

int seqReverse(char* seq, int seqLen)
{
    int i, j;
    char tmp;

    // Step 1, Reverse
    for(i = 0; i < seqLen / 2; i++) {
        tmp = seq[i];
        seq[i] = seq[seqLen - 1 - i];
        seq[seqLen - 1 - i] = tmp;
    }
    // Step 2, Change
    for(i = 0; i < seqLen; i++) {
        switch(seq[i]) {
            case 'A':
                seq[i] = 'T';
                break;
            case 'T':
                seq[i] = 'A';
                break;
            case 'C':
                seq[i] = 'G';
                break;
            case 'G':
                seq[i] = 'C';
                break;
            case 'N':
                seq[i] = 'N';
                break;
            default:
                fprintf(stderr, "Sequence content error!\n");
                for(j = 0; j < seqLen; j++)
                    fprintf(stderr, "%c", seq[j]);
                fprintf(stderr, "\n");
                return 1;
        }
    }

    return 0;
}

void methProcess(char* bamFileName, char* methCgFileName, char* methChgFileName, char* methChhFileName, HashNode** hashTable, char** chrSeqArray, int* chrLen, int chrCnt, int minquali, int mincover, unsigned int mapqThr)
{
    int i, j, off, idx, len;
    int m, n, cnt, iread;
    int* chrDone = (int*)calloc(chrCnt, sizeof(int));

    MapRecord* record = (MapRecord*)malloc(sizeof(MapRecord));
    record->qname = (char*)malloc(sizeof(char) * 1000);
    record->chrome = (char*)malloc(sizeof(char) * 1000);
    record->cigar = (char*)malloc(sizeof(char) * 1000);
    //record->tmp = (char*)malloc(sizeof(char) * 1000);
    record->seq = (char*)malloc(sizeof(char) * 1000);
    record->seqBuf = (char*)malloc(sizeof(char) * 1000);
    record->qual = (char*)malloc(sizeof(char) * 1000);
    record->qualBuf = (char*)malloc(sizeof(char) * 1000);
    record->comBuf = (char*)malloc(sizeof(char) * 1000);

    unsigned int *w_Mm_CG = NULL;
    unsigned int *w_Mc_CG = NULL;
    unsigned int *w_Mq_CG = NULL;
    unsigned int *c_Mm_CG = NULL;
    unsigned int *c_Mc_CG = NULL;
    unsigned int *c_Mq_CG = NULL;

    unsigned int *w_Mm_CHG = NULL;
    unsigned int *w_Mc_CHG = NULL;
    unsigned int *w_Mq_CHG = NULL;
    unsigned int *c_Mm_CHG = NULL;
    unsigned int *c_Mc_CHG = NULL;
    unsigned int *c_Mq_CHG = NULL;

    unsigned int *w_Mm_CHH = NULL;
    unsigned int *w_Mc_CHH = NULL;
    unsigned int *w_Mq_CHH = NULL;
    unsigned int *c_Mm_CHH = NULL;
    unsigned int *c_Mc_CHH = NULL;
    unsigned int *c_Mq_CHH = NULL;

    int rowCnt = 0;
    int finCnt = 0;
    char curChr[100] = "chr1234567890";
    uint8_t *s, *t;
    bam_header_t *header;
    fprintf(stderr, "Meth process begins...\n");

    bamFile in = bam_open(bamFileName, "r");
    if(in == NULL) {
        fprintf(stderr, "Cannot open bam file!\n");
        exit(1);
    }

    FILE* methCgFptr = fopen(methCgFileName, "w");
    if(methCgFptr == NULL) {
        fprintf(stderr, "Could not open file %s!\n", methCgFileName);
        exit(1);
    }
    fprintf(methCgFptr, "#CHROM\tPOS\tCONTEXT\tWatson-METH\tWatson-COVERAGE\tWatson-QUAL\tCrick-METH\tCrick-COVERAGE\tCrick-QUAL\n");

    FILE* methChgFptr = fopen(methChgFileName, "w");
    if(methChgFptr == NULL) {
        fprintf(stderr, "Could not open file %s!\n", methChgFileName);
        exit(1);
    }
    fprintf(methChgFptr, "#CHROM\tPOS\tCONTEXT\tWatson-METH\tWatson-COVERAGE\tWatson-QUAL\tCrick-METH\tCrick-COVERAGE\tCrick-QUAL\n");

    FILE* methChhFptr = fopen(methChhFileName, "w");
    if(methChhFptr == NULL) {
        fprintf(stderr, "Could not open file %s!\n", methChhFileName);
        exit(1);
    }
    fprintf(methChhFptr, "#CHROM\tPOS\tCONTEXT\tWatson-METH\tWatson-COVERAGE\tWatson-QUAL\tCrick-METH\tCrick-COVERAGE\tCrick-QUAL\n");

    bam1_t* b = bam_init1();
    if(b == NULL) {
        fprintf(stderr, "Cannot init bam structure!\n");
        exit(1);
    }

    header = bam_header_read(in);
    while(bam_read1(in, b) >= 0) {
        // Parse record
        if(parseBuffer(header, b, record, mapqThr) == 1)
            continue;

        // Check chrome
        if(strcmp(curChr, record->chrome) != 0) {
            fprintf(stderr, "Chrome updated to %s!\n", record->chrome);
            // Save old chrome statics results
            if(rowCnt > 0) {
                // Update
                finCnt = rowCnt;
                // Print
                printMeth(methCgFptr, len, curChr, w_Mm_CG, w_Mc_CG, w_Mq_CG, c_Mm_CG, c_Mc_CG, c_Mq_CG, "CG", mincover);
                printMeth(methChgFptr, len, curChr, w_Mm_CHG, w_Mc_CHG, w_Mq_CHG, c_Mm_CHG, c_Mc_CHG, c_Mq_CHG, "CHG", mincover);
                printMeth(methChhFptr, len, curChr, w_Mm_CHH, w_Mc_CHH, w_Mq_CHH, c_Mm_CHH, c_Mc_CHH, c_Mq_CHH, "CHH", mincover);
                // Memory gathering for x_Mx_CG
                free(w_Mm_CG);
                free(w_Mc_CG);
                free(w_Mq_CG);
                free(c_Mm_CG);
                free(c_Mc_CG);
                free(c_Mq_CG);
                // Memory gathering for x_Mx_CHG
                free(w_Mm_CHG);
                free(w_Mc_CHG);
                free(w_Mq_CHG);
                free(c_Mm_CHG);
                free(c_Mc_CHG);
                free(c_Mq_CHG);
                // Memory gathering for x_Mx_CHH
                free(w_Mm_CHH);
                free(w_Mc_CHH);
                free(w_Mq_CHH);
                free(c_Mm_CHH);
                free(c_Mc_CHH);
                free(c_Mq_CHH);
            }

            // Update current chrome
            strcpy(curChr, record->chrome);

            // Get chrome length
            idx = hash_table_lookup(hashTable, curChr);
            if(idx == -1) {
                fprintf(stderr, "%s not found in chrome name array.\n", curChr);
                exit(1);
            }
            len = chrLen[idx];

            // Check if chrome has been processed
            if(chrDone[idx] != 0) {
                fprintf(stderr, "%s has already been processed. The bam file is not sorted.\n", curChr);
                exit(1);
            }
            chrDone[idx] = 1;

            // Memory alloction for x_Mx_CG
            if(!(w_Mm_CG = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
                fprintf(stderr, "Not enough memory for w_Mm_CG of %s.\n", curChr);
                exit(1);
            }
            if(!(w_Mc_CG = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
                fprintf(stderr, "Not enough memory for w_Mc_CG of %s.\n", curChr);
                exit(1);
            }
            if(!(w_Mq_CG = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
                fprintf(stderr, "Not enough memory for w_Mq_CG of %s.\n", curChr);
                exit(1);
            }
            if(!(c_Mm_CG = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
                fprintf(stderr, "Not enough memory for c_Mm_CG of %s.\n", curChr);
                exit(1);
            }
            if(!(c_Mc_CG = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
                fprintf(stderr, "Not enough memory for c_Mc_CG of %s.\n", curChr);
                exit(1);
            }
            if(!(c_Mq_CG = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
                fprintf(stderr, "Not enough memory for c_Mq_CG of %s.\n", curChr);
                exit(1);
            }

            // Memory alloction for x_Mx_CHG
            if(!(w_Mm_CHG = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
                fprintf(stderr, "Not enough memory for w_Mm_CG of %s.\n", curChr);
                exit(1);
            }
            if(!(w_Mc_CHG = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
                fprintf(stderr, "Not enough memory for w_Mc_CG of %s.\n", curChr);
                exit(1);
            }
            if(!(w_Mq_CHG = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
                fprintf(stderr, "Not enough memory for w_Mq_CG of %s.\n", curChr);
                exit(1);
            }
            if(!(c_Mm_CHG = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
                fprintf(stderr, "Not enough memory for c_Mm_CG of %s.\n", curChr);
                exit(1);
            }
            if(!(c_Mc_CHG = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
                fprintf(stderr, "Not enough memory for c_Mc_CG of %s.\n", curChr);
                exit(1);
            }
            if(!(c_Mq_CHG = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
                fprintf(stderr, "Not enough memory for c_Mq_CG of %s.\n", curChr);
                exit(1);
            }

            // Memory alloction for x_Mx_CG
            if(!(w_Mm_CHH = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
                fprintf(stderr, "Not enough memory for w_Mm_CG of %s.\n", curChr);
                exit(1);
            }
            if(!(w_Mc_CHH = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
                fprintf(stderr, "Not enough memory for w_Mc_CG of %s.\n", curChr);
                exit(1);
            }
            if(!(w_Mq_CHH = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
                fprintf(stderr, "Not enough memory for w_Mq_CG of %s.\n", curChr);
                exit(1);
            }
            if(!(c_Mm_CHH = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
                fprintf(stderr, "Not enough memory for c_Mm_CG of %s.\n", curChr);
                exit(1);
            }
            if(!(c_Mc_CHH = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
                fprintf(stderr, "Not enough memory for c_Mc_CG of %s.\n", curChr);
                exit(1);
            }
            if(!(c_Mq_CHH = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
                fprintf(stderr, "Not enough memory for c_Mq_CG of %s.\n", curChr);
                exit(1);
            }
        }

        /////////////////////////////////////////////////
        // Filter Step 1: Base quality
        /////////////////////////////////////////////////
        for(i = 0; i < record->len; i++) {
            if(record->seq[i] != 'N' && (unsigned short)(record->qual[i] - 33) < minquali)
                record->seq[i] = 'N';
        }
        /////////////////////////////////////////////////
        // Additional Step 1: Check and record CG methy
        /////////////////////////////////////////////////
        off = record->offset - 1;
        if((record->r12 == 1 && record->strand == '+') || (record->r12 == 2 && record->strand == '-')) {
            // Watson chain
            for(i = 0; i < record->len; i++) {
                if(off + 1 < chrLen[idx]) {
                    if(record->seq[i] == 'C') {
                        if(chrSeqArray[idx][off] == 'C' && chrSeqArray[idx][off+1] == 'G') {
                            w_Mm_CG[off]++;
                            w_Mc_CG[off]++;
                            w_Mq_CG[off] += record->qual[i] - 33;
                        }
                    }
                    else if(record->seq[i] == 'T') {
                        if(chrSeqArray[idx][off] == 'C' && chrSeqArray[idx][off+1] == 'G') {
                            w_Mc_CG[off]++;
                            w_Mq_CG[off] += record->qual[i] - 33;
                        }
                    }
                }
                off++;
            }
        }
        else {
            // Crick chain
            for(i = 0; i < record->len; i++) {
                if(off - 1 >= 0) {
                    if(record->seq[i] == 'G') {
                        if(chrSeqArray[idx][off] == 'G' && chrSeqArray[idx][off-1] == 'C') {
                            c_Mm_CG[off]++;
                            c_Mc_CG[off]++;
                            c_Mq_CG[off] += record->qual[i] - 33;
                        }
                    }
                    else if(record->seq[i] == 'A') {
                        if(chrSeqArray[idx][off] == 'G' && chrSeqArray[idx][off-1] == 'C') {
                            c_Mc_CG[off]++;
                            c_Mq_CG[off] += record->qual[i] - 33;
                        }
                    }
                }
                off++;
            }
        }
        /////////////////////////////////////////////////
        // Additional Step 2: Check and record CHG methy
        /////////////////////////////////////////////////
        off = record->offset - 1;
        if((record->r12 == 1 && record->strand == '+') || (record->r12 == 2 && record->strand == '-')) {
            // Watson chain
            for(i = 0; i < record->len; i++) {
                if(off + 2 < chrLen[idx]) {
                    if(record->seq[i] == 'C') {
                        if(chrSeqArray[idx][off] == 'C' && chrSeqArray[idx][off+1] != 'G' && chrSeqArray[idx][off+2] == 'G') {
                            w_Mm_CHG[off]++;
                            w_Mc_CHG[off]++;
                            w_Mq_CHG[off] += record->qual[i] - 33;
                        }
                    }
                    else if(record->seq[i] == 'T') {
                        if(chrSeqArray[idx][off] == 'C' && chrSeqArray[idx][off+1] != 'G' && chrSeqArray[idx][off+2] == 'G') {
                            w_Mc_CHG[off]++;
                            w_Mq_CHG[off] += record->qual[i] - 33;
                        }
                    }
                }
                off++;
            }
        }
        else {
            // Crick chain
            for(i = 0; i < record->len; i++) {
                if(off - 2 >= 0) {
                    if(record->seq[i] == 'G') {
                        if(chrSeqArray[idx][off] == 'G' && chrSeqArray[idx][off-1] != 'C' && chrSeqArray[idx][off-2] == 'C') {
                            c_Mm_CHG[off]++;
                            c_Mc_CHG[off]++;
                            c_Mq_CHG[off] += record->qual[i] - 33;
                        }
                    }
                    else if(record->seq[i] == 'A') {
                        if(chrSeqArray[idx][off] == 'G' && chrSeqArray[idx][off-1] != 'C' && chrSeqArray[idx][off-2] == 'C') {
                            c_Mc_CHG[off]++;
                            c_Mq_CHG[off] += record->qual[i] - 33;
                        }
                    }
                }
                off++;
            }
        }
        /////////////////////////////////////////////////
        // Additional Step 3: Check and record CHH methy
        /////////////////////////////////////////////////
        off = record->offset - 1;
        if((record->r12 == 1 && record->strand == '+') || (record->r12 == 2 && record->strand == '-')) {
            // Watson chain
            for(i = 0; i < record->len; i++) {
                /*
                if(off == 257561 && strcmp(curChr, "chr1") == 0) {
                        fprintf(stderr, "\n\n");
                        for(j = 0; j < 10; j++)
                            fprintf(stderr, "%c", chrSeqArray[idx][off+j]);
                        fprintf(stderr, "\n\n");
                        fprintf(stderr, "%d, %c, %c, %d\n", record->r12, record->strand, record->seq[i],  chrLen[idx]);
                        exit(1);
                }
                */
                if(off + 2 < chrLen[idx]) {
                    if(record->seq[i] == 'C') {
                        if(chrSeqArray[idx][off] == 'C' && chrSeqArray[idx][off+1] != 'G' && chrSeqArray[idx][off+2] != 'G') {
                            w_Mm_CHH[off]++;
                            w_Mc_CHH[off]++;
                            w_Mq_CHH[off] += record->qual[i] - 33;
                        }
                    }
                    else if(record->seq[i] == 'T') {
                        if(chrSeqArray[idx][off] == 'C' && chrSeqArray[idx][off+1] != 'G' && chrSeqArray[idx][off+2] != 'G') {
                            w_Mc_CHH[off]++;
                            w_Mq_CHH[off] += record->qual[i] - 33;
                        }
                    }
                }
                off++;
            }
        }
        else {
            // Crick chain
            for(i = 0; i < record->len; i++) {
                if(off - 2 >= 0) {
                    if(record->seq[i] == 'G') {
                        if(chrSeqArray[idx][off] == 'G' && chrSeqArray[idx][off-1] != 'C' && chrSeqArray[idx][off-2] != 'C') {
                            c_Mm_CHH[off]++;
                            c_Mc_CHH[off]++;
                            c_Mq_CHH[off] += record->qual[i] - 33;
                        }
                    }
                    else if(record->seq[i] == 'A') {
                        if(chrSeqArray[idx][off] == 'G' && chrSeqArray[idx][off-1] != 'C' && chrSeqArray[idx][off-2] != 'C') {
                            c_Mc_CHH[off]++;
                            c_Mq_CHH[off] += record->qual[i] - 33;
                        }
                    }
                }
                off++;
            }
        }

        // Update row counter
        rowCnt++;
        if(rowCnt % 500000 == 0)
            fprintf(stderr, "%d records have been dealed.\n", rowCnt);
    }

    // Last batch
    if(finCnt < rowCnt) {
        // Print
        printMeth(methCgFptr, len, curChr, w_Mm_CG, w_Mc_CG, w_Mq_CG, c_Mm_CG, c_Mc_CG, c_Mq_CG, "CG", mincover);
        printMeth(methChgFptr, len, curChr, w_Mm_CHG, w_Mc_CHG, w_Mq_CHG, c_Mm_CHG, c_Mc_CHG, c_Mq_CHG, "CHG", mincover);
        printMeth(methChhFptr, len, curChr, w_Mm_CHH, w_Mc_CHH, w_Mq_CHH, c_Mm_CHH, c_Mc_CHH, c_Mq_CHH, "CHH", mincover);
        // Memory gathering for x_Mx_CG
        free(w_Mm_CG);
        free(w_Mc_CG);
        free(w_Mq_CG);
        free(c_Mm_CG);
        free(c_Mc_CG);
        free(c_Mq_CG);
        // Memory gathering for x_Mx_CHG
        free(w_Mm_CHG);
        free(w_Mc_CHG);
        free(w_Mq_CHG);
        free(c_Mm_CHG);
        free(c_Mc_CHG);
        free(c_Mq_CHG);
        // Memory gathering for x_Mx_CHH
        free(w_Mm_CHH);
        free(w_Mc_CHH);
        free(w_Mq_CHH);
        free(c_Mm_CHH);
        free(c_Mc_CHH);
        free(c_Mq_CHH);
    }

    bam_header_destroy(header);
    bam_close(in);
    bam_destroy1(b);
    fclose(methCgFptr);
    fclose(methChgFptr);
    fclose(methChhFptr);

    fprintf(stderr, "Meth process ends...\n");
}


void memerror(int nline){
    fprintf(stderr, "Not enough memory for SNP of %d.\n", nline);
}

//void init_mem_snp(unsigned short *w_A, unsigned short *w_T ,unsigned short *w_C ,unsigned short *w_G ,unsigned short *c_A ,
//    unsigned short *c_T ,unsigned short *c_C ,unsigned short *c_G, unsigned short *w_Aq ,unsigned short *w_Tq ,unsigned short *w_Cq ,
//    unsigned short *w_Gq , unsigned short *c_Aq ,unsigned short *c_Tq ,unsigned short *c_Cq ,unsigned short *c_Gq, unsigned short *w_Q, 
//    unsigned short *c_Q, unsigned int lenchr)
void init_mem_snp(unsigned short *&w_A, unsigned short *&w_T ,unsigned short *&w_C ,unsigned short *&w_G ,unsigned short *&c_A ,
    unsigned short *&c_T ,unsigned short *&c_C ,unsigned short *&c_G, unsigned short *&w_Aq ,unsigned short *&w_Tq ,unsigned short *&w_Cq ,
    unsigned short *&w_Gq , unsigned short *&c_Aq ,unsigned short *&c_Tq ,unsigned short *&c_Cq ,unsigned short *&c_Gq, unsigned short *&w_Q, 
    unsigned short *&c_Q, unsigned int lenchr){
    fprintf(stderr, "Initial SNP array, len: %d\n", lenchr);
    unsigned int i=0;
    for(i = 0; i < lenchr; i++){
        w_A[i]=0;w_T[i]=0;w_C[i]=0;w_G[i]=0;c_A[i]=0;c_T[i]=0;c_C[i]=0;c_G[i]=0;
        w_Aq[i]=0;w_Tq[i]=0;w_Cq[i]=0;w_Gq[i]=0;c_Aq[i]=0;c_Tq[i]=0;c_Cq[i]=0;c_Gq[i]=0;
        w_Q[i]=0;c_Q[i]=0;
    }
}

void creat_mem_snp(unsigned short *&w_A, unsigned short *&w_T ,unsigned short *&w_C ,unsigned short *&w_G ,unsigned short *&c_A ,
    unsigned short *&c_T ,unsigned short *&c_C ,unsigned short *&c_G, unsigned short *&w_Aq ,unsigned short *&w_Tq ,unsigned short *&w_Cq ,
    unsigned short *&w_Gq , unsigned short *&c_Aq ,unsigned short *&c_Tq ,unsigned short *&c_Cq ,unsigned short *&c_Gq, unsigned short *&w_Q, 
    unsigned short *&c_Q, unsigned int lenchr){
    fprintf(stderr, "Memory alloction for SNP/meth, array len: %u\n", lenchr);
    // Memory alloction for x_X
    if(!(w_A = (unsigned short*)calloc(lenchr, sizeof(unsigned short)))) {
        memerror(1);
        exit(1);
    }
    if(!(w_T = (unsigned short*)calloc(lenchr, sizeof(unsigned short)))) {
        memerror(2);
        exit(1);
    }
    if(!(w_C = (unsigned short*)calloc(lenchr, sizeof(unsigned short)))) {
        memerror(3);
        exit(1);
    }
    if(!(w_G = (unsigned short*)calloc(lenchr, sizeof(unsigned short)))) {
        memerror(4);
        exit(1);
    }
    if(!(c_A = (unsigned short*)calloc(lenchr, sizeof(unsigned short)))) {
        memerror(5);
        exit(1);
    }
    if(!(c_T = (unsigned short*)calloc(lenchr, sizeof(unsigned short)))) {
        memerror(6);
        exit(1);
    }
    if(!(c_C = (unsigned short*)calloc(lenchr, sizeof(unsigned short)))) {
        memerror(7);
        exit(1);
    }
    if(!(c_G = (unsigned short*)calloc(lenchr, sizeof(unsigned short)))) {
        memerror(8);
        exit(1);
    }

    // Memory alloction for x_Xq
    if(!(w_Aq = (unsigned short*)calloc(lenchr, sizeof(unsigned short)))) {
        memerror(9);
        exit(1);
    }
    if(!(w_Tq = (unsigned short*)calloc(lenchr, sizeof(unsigned short)))) {
        memerror(10);
        exit(1);
    }
    if(!(w_Cq = (unsigned short*)calloc(lenchr, sizeof(unsigned short)))) {
        memerror(11);
        exit(1);
    }
    if(!(w_Gq = (unsigned short*)calloc(lenchr, sizeof(unsigned short)))) {
        memerror(12);
        exit(1);
    }
    if(!(c_Aq = (unsigned short*)calloc(lenchr, sizeof(unsigned short)))) {
        memerror(13);
        exit(1);
    }
    if(!(c_Tq = (unsigned short*)calloc(lenchr, sizeof(unsigned short)))) {
        memerror(14);
        exit(1);
    }
    if(!(c_Cq = (unsigned short*)calloc(lenchr, sizeof(unsigned short)))) {
        memerror(15);
        exit(1);
    }
    if(!(c_Gq = (unsigned short*)calloc(lenchr, sizeof(unsigned short)))) {
        memerror(16);
        exit(1);
    }

    // Memory alloction for x_Q
    if(!(w_Q = (unsigned short*)calloc(lenchr, sizeof(unsigned short)))) {
        memerror(25);
        exit(1);
    }
    if(!(c_Q = (unsigned short*)calloc(lenchr, sizeof(unsigned short)))) {
        memerror(26);
        exit(1);
    }
}

void snpProcess(FILE* methFptr, FILE* snpFptr, char* bamFileName, HashNode** hashTable, char** chrSeqArray, int* chrLen, int chrCnt, int minquali, 
    int maxcover, float minhetfreq, float errorrate, unsigned int mapqThr, char* processChrom)
{
    int i, j, off, idx, len;
    int m, n, cnt, iread;
    int* chrDone = (int*)calloc(chrCnt, sizeof(int));

    MapRecord* record = (MapRecord*)malloc(sizeof(MapRecord));
    record->qname = (char*)malloc(sizeof(char) * 1000);
    record->chrome = (char*)malloc(sizeof(char) * 1000);
    record->cigar = (char*)malloc(sizeof(char) * 1000);
    //record->tmp = (char*)malloc(sizeof(char) * 1000);
    record->seq = (char*)malloc(sizeof(char) * 1000);
    record->seqBuf = (char*)malloc(sizeof(char) * 1000);
    record->qual = (char*)malloc(sizeof(char) * 1000);
    record->qualBuf = (char*)malloc(sizeof(char) * 1000);
    record->comBuf = (char*)malloc(sizeof(char) * 1000);

    unsigned short *w_A = NULL;
    unsigned short *w_T = NULL;
    unsigned short *w_C = NULL;
    unsigned short *w_G = NULL;
    unsigned short *c_A = NULL;
    unsigned short *c_T = NULL;
    unsigned short *c_C = NULL;
    unsigned short *c_G = NULL;

    unsigned short *w_Aq = NULL;
    unsigned short *w_Tq = NULL;
    unsigned short *w_Cq = NULL;
    unsigned short *w_Gq = NULL;
    unsigned short *c_Aq = NULL;
    unsigned short *c_Tq = NULL;
    unsigned short *c_Cq = NULL;
    unsigned short *c_Gq = NULL;

    unsigned short *w_Q = NULL;
    unsigned short *c_Q = NULL;

    idx = hash_table_lookup(hashTable, processChrom);
    if(idx == -1) {
        fprintf(stderr, "%s not found in chrome name array, chr len file.\n", processChrom);
        exit(1);
    }

    int rowCnt = 0;
    int finCnt = 0;
    char curChr[100] = "chr1234567890";
    uint8_t *s, *t;
    bam_header_t *header;
    
    samFile *infp=0;
    if ((infp = sam_open(bamFileName)) == 0) { //, "r"
        fprintf(stderr, "Fail to open BAM file %s\n", bamFileName);
        exit(1);
    }

    //bam index
    hts_idx_t *bamidx = sam_index_load(infp, bamFileName); // load index
    if (bamidx == 0) { // index is unavailable
        fprintf(stderr, "Muiltype threads only works for indexed BAM files, please run samtools index %s first.\n", bamFileName);
        exit(0);
    }
    header = sam_hdr_read(infp);

    hts_itr_t *iter = sam_itr_querys(bamidx, header, processChrom);
    if (iter == NULL) { // region invalid or reference name not found
        return;
    }

    bam1_t* b = bam_init1();
    if(b == NULL) {
        fprintf(stderr, "Cannot init bam structure!\n");
        exit(1);
    }

    len = chrLen[idx];
    fprintf(stderr, "Processing chromosome: %s\n", processChrom);
    creat_mem_snp(w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, w_Aq, w_Tq, w_Cq, w_Gq, c_Aq, c_Tq, c_Cq, c_Gq, w_Q, c_Q, len);
    
    while( sam_itr_next(infp, iter, b) >= 0 ) {
        // Parse record
        if(parseBuffer(header, b, record, mapqThr) == 1)
            continue;
        // Check chrome
        if(strcmp(curChr, record->chrome) != 0) {
            // Save old chrome statics results
            
            if(rowCnt > 0) {
                // should not be here
                // Update
                fprintf(stderr, "\ncan not here!\n");
                exit(0);
                finCnt = rowCnt;
                // Print
                printSnp(methFptr, snpFptr, chrSeqArray, idx, len, minhetfreq, curChr, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, w_Aq, w_Tq, w_Cq, w_Gq, c_Aq, c_Tq, c_Cq, c_Gq,  w_Q, c_Q);
                // Memory init
                init_mem_snp(w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, w_Aq, w_Tq, w_Cq, w_Gq, c_Aq, c_Tq, c_Cq, c_Gq, w_Q, c_Q, len);
            }
            

            // Update current chrome
            strcpy(curChr, record->chrome);

            // Get chrome length
            idx = hash_table_lookup(hashTable, curChr);
            if(idx == -1) {
                fprintf(stderr, "%s not found in chrome name array.\n", curChr);
                exit(1);
            }
            len = chrLen[idx];

            // Check if chrome has been processed
            if(chrDone[idx] != 0) {
                fprintf(stderr, "%s has already been processed. The bam file is not sorted.\n", curChr);
                exit(1);
            }
            chrDone[idx] = 1;
        }

        /////////////////////////////////////////////////
        // Filter Step 1: Base quality
        /////////////////////////////////////////////////
        for(i = 0; i < record->len; i++) {
            if(record->seq[i] != 'N' && (unsigned short)(record->qual[i] - 33) < minquali)
                record->seq[i] = 'N';
        }
        /////////////////////////////////////////////////
        // Filter Step 2: SNP per base
        /////////////////////////////////////////////////
        //printf("\n%s %d %s %d\n", record->qname, record->len, record->seq, record->offset);
        off = record->offset - 1;
        iread = (int)(record->len * errorrate + 0.5);
        if(iread < 1)
            iread = 1;
        cnt = 0;
        for(i = 0; i < record->len; i++) {
            if(record->strand == '+') {
                switch(record->seq[i]) {
                    case 'A':
                        if(chrSeqArray[idx][off] != 'A' && (!(record->r12 == 2 && chrSeqArray[idx][off] == 'G')))
                            cnt++;
                        break;
                    case 'T':
                        if(chrSeqArray[idx][off] != 'T' && (!(record->r12 == 1 && chrSeqArray[idx][off] == 'C')))
                            cnt++;
                        break;
                    case 'C':
                        if(chrSeqArray[idx][off] != 'C')
                            cnt++;
                        break;
                    case 'G':
                        if(chrSeqArray[idx][off] != 'G')
                            cnt++;
                        break;
                }
            }
            else {
                switch(record->seq[i]) {
                    case 'A':
                        if(chrSeqArray[idx][off] != 'A' && (!(record->r12 == 1 && chrSeqArray[idx][off] == 'G')))
                            cnt++;
                        break;
                    case 'T':
                        if(chrSeqArray[idx][off] != 'T' && (!(record->r12 == 2 && chrSeqArray[idx][off] == 'C')))
                            cnt++;
                        break;
                    case 'C':
                        if(chrSeqArray[idx][off] != 'C')
                            cnt++;
                        break;
                    case 'G':
                        if(chrSeqArray[idx][off] != 'G')
                            cnt++;
                        break;
                }
            }

            off++;
        }

        // Snp per read
        //printf("\n%d < %d %s %f %s %d\n", cnt, iread, record->qname, errorrate, record->qual, record->offset - 1);
        if(cnt <= iread) {
            off = record->offset - 1;
            for(i = 0; i < record->len; i++) {
                if(record->strand == '+') {
                    switch(record->seq[i]) {
                        case 'A':
                            if(0 && record->r12 ==2 && chrSeqArray[idx][off] == 'G') {
                                c_G[off]++;
                                c_Gq[off] += (unsigned short)(record->qual[i] - 33);
                                c_Q[off] += record->mapq;
                            }
                            else {
                                if(record->r12 == 1) {
                                    w_A[off]++;
                                    w_Aq[off] += (unsigned short)(record->qual[i] - 33);
                                    w_Q[off] += record->mapq;
                                }
                                else {
                                    c_A[off]++;
                                    c_Aq[off] += (unsigned short)(record->qual[i] - 33);
                                    c_Q[off] += record->mapq;
                                }
                            }
                            break;
                        case 'T':
                            if(0 &&record->r12 == 1 && chrSeqArray[idx][off] == 'C') {
                                w_C[off]++;
                                w_Cq[off] += (unsigned short)(record->qual[i] - 33);
                                w_Q[off] += record->mapq;
                            }
                            else {
                                if(record->r12 == 1) {
                                    w_T[off]++;
                                    w_Tq[off] += (unsigned short)(record->qual[i] - 33);
                                    w_Q[off] += record->mapq;
                                 }
                                else {
                                    c_T[off]++;
                                    c_Tq[off] += (unsigned short)(record->qual[i] - 33);
                                    c_Q[off] += record->mapq;
                                }
                            }
                            break;
                        case 'C':
                            if(record->r12 == 1) {
                                w_C[off]++;
                                w_Cq[off] += (unsigned short)(record->qual[i] - 33);
                                w_Q[off] += record->mapq;
                            }
                            else {
                                c_C[off]++;
                                c_Cq[off] += (unsigned short)(record->qual[i] - 33);
                                c_Q[off] += record->mapq;
                            }
                            break;
                        case 'G':
                            if(record->r12 == 1) {
                                w_G[off]++;
                                w_Gq[off] += (unsigned short)(record->qual[i] - 33);
                                w_Q[off] += record->mapq;
                            }
                            else {
                                c_G[off]++;
                                c_Gq[off] += (unsigned short)(record->qual[i] - 33);
                                c_Q[off] += record->mapq;
                            }
                            break;
                    }
                }
                else {
                    switch(record->seq[i]) {
                        case 'A':
                            if(0 &&record->r12 == 1 && chrSeqArray[idx][off] == 'G') {
                                c_G[off]++;
                                c_Gq[off] += (unsigned short)(record->qual[i] - 33);
                                c_Q[off] += record->mapq;
                            }
                            else {
                                if(record->r12 == 2) {
                                    w_A[off]++;
                                    w_Aq[off] += (unsigned short)(record->qual[i] - 33);
                                    w_Q[off] += record->mapq;
                                }
                                else {
                                    c_A[off]++;
                                    c_Aq[off] += (unsigned short)(record->qual[i] - 33);
                                    c_Q[off] += record->mapq;
                                }
                            }
                            break;
                        case 'T':
                            if(0 &&record->r12 == 2 && chrSeqArray[idx][off] == 'C') {
                                w_C[off]++;
                                w_Cq[off] += (unsigned short)(record->qual[i] - 33);
                                w_Q[off] += record->mapq;
                            }
                            else {
                                if(record->r12 == 2) {
                                    w_T[off]++;
                                    w_Tq[off] += (unsigned short)(record->qual[i] - 33);
                                    w_Q[off] += record->mapq;
                                }
                                else {
                                    c_T[off]++;
                                    c_Tq[off] += (unsigned short)(record->qual[i] - 33);
                                    c_Q[off] += record->mapq;
                                }
                            }
                            break;
                        case 'C':
                            if(record->r12 == 2) {
                                w_C[off]++;
                                w_Cq[off] += (unsigned short)(record->qual[i] - 33);
                                w_Q[off] += record->mapq;
                            }
                            else {
                                c_C[off]++;
                                c_Cq[off] += (unsigned short)(record->qual[i] - 33);
                                c_Q[off] += record->mapq;
                            }
                            break;
                        case 'G':
                            if(record->r12 == 2) {
                                w_G[off]++;
                                w_Gq[off] += (unsigned short)(record->qual[i] - 33);
                                w_Q[off] += record->mapq;
                            }
                            else {
                                c_G[off]++;
                                c_Gq[off] += (unsigned short)(record->qual[i] - 33);
                                c_Q[off] += record->mapq;
                            }
                            break;
                    }
                }

                off++;
            }
        }
        else {
            #ifdef __MY_DEBUG__
            if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && record->offset <= __DEBUG_POS__ && record->offset+record->len >= __DEBUG_POS__) {
                fprintf(stderr, "read discarded because snp-ratio overflow. [threshold %d, actually %d.]\n", iread, cnt);
                #ifdef __MY_DEBUG_READ__
                dispRecord(record);
                #endif
            }
            #endif
        }
        // Update row counter
        rowCnt++;
        if(rowCnt % 100000 == 0)
            fprintf(stderr, "%s: %d reads have been dealed.\n", curChr, rowCnt);
    }

    // Last batch
    if(finCnt < rowCnt) {
        // Print
        pthread_mutex_lock(&output_mutex);
        printSnp(methFptr, snpFptr, chrSeqArray, idx, len, minhetfreq, curChr, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, w_Aq, w_Tq, w_Cq, w_Gq, c_Aq, c_Tq, c_Cq, c_Gq, w_Q, c_Q);
        pthread_mutex_unlock(&output_mutex);
        fprintf(stderr, "Memory free for SNP/meth, array len: %u\n", len);
        // Memory gathering for x_X
        free(w_A);
        free(w_T);
        free(w_C);
        free(w_G);
        free(c_A);
        free(c_T);
        free(c_C);
        free(c_G);
        // Memory gathering for x_Xq
        free(w_Aq);
        free(w_Tq);
        free(w_Cq);
        free(w_Gq);
        free(c_Aq);
        free(c_Tq);
        free(c_Cq);
        free(c_Gq);

        // Memory gathering for x_Q
        free(w_Q);
        free(c_Q);
    }

    bam_header_destroy(header);
    //bam_close(in);
    sam_close(infp);
    bam_destroy1(b);
    hts_itr_destroy(iter);
    
}

//muti-outfile
void snpProcess_multiop(FILE* methFptr, FILE* snpFptr, char* bamFileName, HashNode** hashTable, char** chrSeqArray, int* chrLen, int chrCnt, int minquali, 
    int maxcover, float minhetfreq, float errorrate, unsigned int mapqThr, char* processChrom)
{
    int i, j, off, idx, len;
    int m, n, cnt, iread;
    int* chrDone = (int*)calloc(chrCnt, sizeof(int));

    MapRecord* record = (MapRecord*)malloc(sizeof(MapRecord));
    record->qname = (char*)malloc(sizeof(char) * 1000);
    record->chrome = (char*)malloc(sizeof(char) * 1000);
    record->cigar = (char*)malloc(sizeof(char) * 1000);
    //record->tmp = (char*)malloc(sizeof(char) * 1000);
    record->seq = (char*)malloc(sizeof(char) * 1000);
    record->seqBuf = (char*)malloc(sizeof(char) * 1000);
    record->qual = (char*)malloc(sizeof(char) * 1000);
    record->qualBuf = (char*)malloc(sizeof(char) * 1000);
    record->comBuf = (char*)malloc(sizeof(char) * 1000);

    unsigned short *w_A = NULL;
    unsigned short *w_T = NULL;
    unsigned short *w_C = NULL;
    unsigned short *w_G = NULL;
    unsigned short *c_A = NULL;
    unsigned short *c_T = NULL;
    unsigned short *c_C = NULL;
    unsigned short *c_G = NULL;

    unsigned short *w_Aq = NULL;
    unsigned short *w_Tq = NULL;
    unsigned short *w_Cq = NULL;
    unsigned short *w_Gq = NULL;
    unsigned short *c_Aq = NULL;
    unsigned short *c_Tq = NULL;
    unsigned short *c_Cq = NULL;
    unsigned short *c_Gq = NULL;

    unsigned short *w_Q = NULL;
    unsigned short *c_Q = NULL;

    idx = hash_table_lookup(hashTable, processChrom);
    if(idx == -1) {
        fprintf(stderr, "%s not found in chrome name array, chr len file.\n", processChrom);
        exit(1);
    }

    int rowCnt = 0;
    int finCnt = 0;
    char curChr[100] = "chr1234567890";
    uint8_t *s, *t;
    bam_header_t *header;
    
    samFile *infp=0;
    if ((infp = sam_open(bamFileName)) == 0) { //, "r"
        fprintf(stderr, "Fail to open BAM file %s\n", bamFileName);
        exit(1);
    }

    //bam index
    hts_idx_t *bamidx = sam_index_load(infp, bamFileName); // load index
    if (bamidx == 0) { // index is unavailable
        fprintf(stderr, "Muiltype threads only works for indexed BAM files, please run samtools index %s first.\n", bamFileName);
        exit(0);
    }
    header = sam_hdr_read(infp);

    hts_itr_t *iter = sam_itr_querys(bamidx, header, processChrom);
    if (iter == NULL) { // region invalid or reference name not found
        return;
    }

    bam1_t* b = bam_init1();
    if(b == NULL) {
        fprintf(stderr, "Cannot init bam structure!\n");
        exit(1);
    }

    len = chrLen[idx];
    fprintf(stderr, "Processing chromosome: %s\n", processChrom);
    creat_mem_snp(w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, w_Aq, w_Tq, w_Cq, w_Gq, c_Aq, c_Tq, c_Cq, c_Gq, w_Q, c_Q, len);
    
    while( sam_itr_next(infp, iter, b) >= 0 ) {
        // Parse record
        if(parseBuffer(header, b, record, mapqThr) == 1)
            continue;
        // Check chrome
        if(strcmp(curChr, record->chrome) != 0) {
            // Save old chrome statics results
            
            if(rowCnt > 0) {
                // should not be here
                // Update
                fprintf(stderr, "\ncan not here!\n");
                exit(0);
                finCnt = rowCnt;
                // Print
                printSnp(methFptr, snpFptr, chrSeqArray, idx, len, minhetfreq, curChr, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, w_Aq, w_Tq, w_Cq, w_Gq, c_Aq, c_Tq, c_Cq, c_Gq,  w_Q, c_Q);
                // Memory init
                init_mem_snp(w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, w_Aq, w_Tq, w_Cq, w_Gq, c_Aq, c_Tq, c_Cq, c_Gq, w_Q, c_Q, len);
            }
            

            // Update current chrome
            strcpy(curChr, record->chrome);

            // Get chrome length
            idx = hash_table_lookup(hashTable, curChr);
            if(idx == -1) {
                fprintf(stderr, "%s not found in chrome name array.\n", curChr);
                exit(1);
            }
            len = chrLen[idx];

            // Check if chrome has been processed
            if(chrDone[idx] != 0) {
                fprintf(stderr, "%s has already been processed. The bam file is not sorted.\n", curChr);
                exit(1);
            }
            chrDone[idx] = 1;
        }

        /////////////////////////////////////////////////
        // Filter Step 1: Base quality
        /////////////////////////////////////////////////
        for(i = 0; i < record->len; i++) {
            if(record->seq[i] != 'N' && (unsigned short)(record->qual[i] - 33) < minquali)
                record->seq[i] = 'N';
        }
        /////////////////////////////////////////////////
        // Filter Step 2: SNP per base
        /////////////////////////////////////////////////
        //printf("\n%s %d %s %d\n", record->qname, record->len, record->seq, record->offset);
        off = record->offset - 1;
        iread = (int)(record->len * errorrate + 0.5);
        if(iread < 1)
            iread = 1;
        cnt = 0;
        for(i = 0; i < record->len; i++) {
            if(record->strand == '+') {
                switch(record->seq[i]) {
                    case 'A':
                        if(chrSeqArray[idx][off] != 'A' && (!(record->r12 == 2 && chrSeqArray[idx][off] == 'G')))
                            cnt++;
                        break;
                    case 'T':
                        if(chrSeqArray[idx][off] != 'T' && (!(record->r12 == 1 && chrSeqArray[idx][off] == 'C')))
                            cnt++;
                        break;
                    case 'C':
                        if(chrSeqArray[idx][off] != 'C')
                            cnt++;
                        break;
                    case 'G':
                        if(chrSeqArray[idx][off] != 'G')
                            cnt++;
                        break;
                }
            }
            else {
                switch(record->seq[i]) {
                    case 'A':
                        if(chrSeqArray[idx][off] != 'A' && (!(record->r12 == 1 && chrSeqArray[idx][off] == 'G')))
                            cnt++;
                        break;
                    case 'T':
                        if(chrSeqArray[idx][off] != 'T' && (!(record->r12 == 2 && chrSeqArray[idx][off] == 'C')))
                            cnt++;
                        break;
                    case 'C':
                        if(chrSeqArray[idx][off] != 'C')
                            cnt++;
                        break;
                    case 'G':
                        if(chrSeqArray[idx][off] != 'G')
                            cnt++;
                        break;
                }
            }

            off++;
        }

        // Snp per read
        //printf("\n%d < %d %s %f %s %d\n", cnt, iread, record->qname, errorrate, record->qual, record->offset - 1);
        if(cnt <= iread) {
            off = record->offset - 1;
            for(i = 0; i < record->len; i++) {
                if(record->strand == '+') {
                    switch(record->seq[i]) {
                        case 'A':
                            if(0 && record->r12 ==2 && chrSeqArray[idx][off] == 'G') {
                                c_G[off]++;
                                c_Gq[off] += (unsigned short)(record->qual[i] - 33);
                                c_Q[off] += record->mapq;
                            }
                            else {
                                if(record->r12 == 1) {
                                    w_A[off]++;
                                    w_Aq[off] += (unsigned short)(record->qual[i] - 33);
                                    w_Q[off] += record->mapq;
                                }
                                else {
                                    c_A[off]++;
                                    c_Aq[off] += (unsigned short)(record->qual[i] - 33);
                                    c_Q[off] += record->mapq;
                                }
                            }
                            break;
                        case 'T':
                            if(0 &&record->r12 == 1 && chrSeqArray[idx][off] == 'C') {
                                w_C[off]++;
                                w_Cq[off] += (unsigned short)(record->qual[i] - 33);
                                w_Q[off] += record->mapq;
                            }
                            else {
                                if(record->r12 == 1) {
                                    w_T[off]++;
                                    w_Tq[off] += (unsigned short)(record->qual[i] - 33);
                                    w_Q[off] += record->mapq;
                                 }
                                else {
                                    c_T[off]++;
                                    c_Tq[off] += (unsigned short)(record->qual[i] - 33);
                                    c_Q[off] += record->mapq;
                                }
                            }
                            break;
                        case 'C':
                            if(record->r12 == 1) {
                                w_C[off]++;
                                w_Cq[off] += (unsigned short)(record->qual[i] - 33);
                                w_Q[off] += record->mapq;
                            }
                            else {
                                c_C[off]++;
                                c_Cq[off] += (unsigned short)(record->qual[i] - 33);
                                c_Q[off] += record->mapq;
                            }
                            break;
                        case 'G':
                            if(record->r12 == 1) {
                                w_G[off]++;
                                w_Gq[off] += (unsigned short)(record->qual[i] - 33);
                                w_Q[off] += record->mapq;
                            }
                            else {
                                c_G[off]++;
                                c_Gq[off] += (unsigned short)(record->qual[i] - 33);
                                c_Q[off] += record->mapq;
                            }
                            break;
                    }
                }
                else {
                    switch(record->seq[i]) {
                        case 'A':
                            if(0 &&record->r12 == 1 && chrSeqArray[idx][off] == 'G') {
                                c_G[off]++;
                                c_Gq[off] += (unsigned short)(record->qual[i] - 33);
                                c_Q[off] += record->mapq;
                            }
                            else {
                                if(record->r12 == 2) {
                                    w_A[off]++;
                                    w_Aq[off] += (unsigned short)(record->qual[i] - 33);
                                    w_Q[off] += record->mapq;
                                }
                                else {
                                    c_A[off]++;
                                    c_Aq[off] += (unsigned short)(record->qual[i] - 33);
                                    c_Q[off] += record->mapq;
                                }
                            }
                            break;
                        case 'T':
                            if(0 &&record->r12 == 2 && chrSeqArray[idx][off] == 'C') {
                                w_C[off]++;
                                w_Cq[off] += (unsigned short)(record->qual[i] - 33);
                                w_Q[off] += record->mapq;
                            }
                            else {
                                if(record->r12 == 2) {
                                    w_T[off]++;
                                    w_Tq[off] += (unsigned short)(record->qual[i] - 33);
                                    w_Q[off] += record->mapq;
                                }
                                else {
                                    c_T[off]++;
                                    c_Tq[off] += (unsigned short)(record->qual[i] - 33);
                                    c_Q[off] += record->mapq;
                                }
                            }
                            break;
                        case 'C':
                            if(record->r12 == 2) {
                                w_C[off]++;
                                w_Cq[off] += (unsigned short)(record->qual[i] - 33);
                                w_Q[off] += record->mapq;
                            }
                            else {
                                c_C[off]++;
                                c_Cq[off] += (unsigned short)(record->qual[i] - 33);
                                c_Q[off] += record->mapq;
                            }
                            break;
                        case 'G':
                            if(record->r12 == 2) {
                                w_G[off]++;
                                w_Gq[off] += (unsigned short)(record->qual[i] - 33);
                                w_Q[off] += record->mapq;
                            }
                            else {
                                c_G[off]++;
                                c_Gq[off] += (unsigned short)(record->qual[i] - 33);
                                c_Q[off] += record->mapq;
                            }
                            break;
                    }
                }

                off++;
            }
        }
        else {
            #ifdef __MY_DEBUG__
            if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && record->offset <= __DEBUG_POS__ && record->offset+record->len >= __DEBUG_POS__) {
                fprintf(stderr, "read discarded because snp-ratio overflow. [threshold %d, actually %d.]\n", iread, cnt);
                #ifdef __MY_DEBUG_READ__
                dispRecord(record);
                #endif
            }
            #endif
        }
        // Update row counter
        rowCnt++;
        if(rowCnt % 100000 == 0)
            fprintf(stderr, "%s: %d reads have been dealed.\n", curChr, rowCnt);
    }

    // Last batch
    if(finCnt < rowCnt) {
        // Print
        printSnp(methFptr, snpFptr, chrSeqArray, idx, len, minhetfreq, curChr, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, w_Aq, w_Tq, w_Cq, w_Gq, c_Aq, c_Tq, c_Cq, c_Gq, w_Q, c_Q);
        // Memory gathering for x_X
        fprintf(stderr, "Memory free for SNP/meth, array len: %u\n", len);
        free(w_A);
        free(w_T);
        free(w_C);
        free(w_G);
        free(c_A);
        free(c_T);
        free(c_C);
        free(c_G);
        // Memory gathering for x_Xq
        free(w_Aq);
        free(w_Tq);
        free(w_Cq);
        free(w_Gq);
        free(c_Aq);
        free(c_Tq);
        free(c_Cq);
        free(c_Gq);

        // Memory gathering for x_Q
        free(w_Q);
        free(c_Q);
    }

    bam_header_destroy(header);
    //bam_close(in);
    sam_close(infp);
    bam_destroy1(b);
    hts_itr_destroy(iter);
    
}

//single thread
void snpProcess_singlet(FILE* methFptr, FILE* snpFptr, char* bamFileName, HashNode** hashTable, char** chrSeqArray, int* chrLen, int chrCnt, int minquali, int maxcover, float minhetfreq, float errorrate, unsigned int mapqThr)
{
    int i, j, off, idx, len;
    int m, n, cnt, iread;
    int* chrDone = (int*)calloc(chrCnt, sizeof(int));

    MapRecord* record = (MapRecord*)malloc(sizeof(MapRecord));
    record->qname = (char*)malloc(sizeof(char) * 1000);
    record->chrome = (char*)malloc(sizeof(char) * 1000);
    record->cigar = (char*)malloc(sizeof(char) * 1000);
    //record->tmp = (char*)malloc(sizeof(char) * 1000);
    record->seq = (char*)malloc(sizeof(char) * 1000);
    record->seqBuf = (char*)malloc(sizeof(char) * 1000);
    record->qual = (char*)malloc(sizeof(char) * 1000);
    record->qualBuf = (char*)malloc(sizeof(char) * 1000);
    record->comBuf = (char*)malloc(sizeof(char) * 1000);

    unsigned short *w_A = NULL;
    unsigned short *w_T = NULL;
    unsigned short *w_C = NULL;
    unsigned short *w_G = NULL;
    unsigned short *c_A = NULL;
    unsigned short *c_T = NULL;
    unsigned short *c_C = NULL;
    unsigned short *c_G = NULL;

    unsigned short *w_Aq = NULL;
    unsigned short *w_Tq = NULL;
    unsigned short *w_Cq = NULL;
    unsigned short *w_Gq = NULL;
    unsigned short *c_Aq = NULL;
    unsigned short *c_Tq = NULL;
    unsigned short *c_Cq = NULL;
    unsigned short *c_Gq = NULL;

    unsigned short *w_Q = NULL;
    unsigned short *c_Q = NULL;

    int rowCnt = 0;
    int finCnt = 0;
    char curChr[100] = "chr1234567890";
    uint8_t *s, *t;
    bam_header_t *header;

    bamFile in = bam_open(bamFileName, "r");
    if(in == NULL) {
        fprintf(stderr, "Cannot open bam file!\n");
        exit(1);
    }

    creat_mem_snp(w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, w_Aq, w_Tq, w_Cq, w_Gq, c_Aq, c_Tq, c_Cq, c_Gq, w_Q, c_Q, longestchr);

    bam1_t* b = bam_init1();
    if(b == NULL) {
        fprintf(stderr, "Cannot init bam structure!\n");
        exit(1);
    }
    header = bam_header_read(in);
    while( bam_read1(in, b) >= 0) {
        // Parse record
        if(parseBuffer(header, b, record, mapqThr) == 1)
            continue;

        // Check chrome
        if(strcmp(curChr, record->chrome) != 0) {
            // Save old chrome statics results
            if(rowCnt > 0) {
                // Update
                finCnt = rowCnt;
                // Print
                printSnp(methFptr, snpFptr, chrSeqArray, idx, len, minhetfreq, curChr, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, w_Aq, w_Tq, w_Cq, w_Gq, c_Aq, c_Tq, c_Cq, c_Gq, w_Q, c_Q);
                // Memory init
                init_mem_snp(w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, w_Aq, w_Tq, w_Cq, w_Gq, c_Aq, c_Tq, c_Cq, c_Gq, w_Q, c_Q, longestchr);
            }

            // Update current chrome
            strcpy(curChr, record->chrome);

            // Get chrome length
            idx = hash_table_lookup(hashTable, curChr);
            if(idx == -1) {
                fprintf(stderr, "%s not found in chrome name array.\n", curChr);
                exit(1);
            }
            len = chrLen[idx];

            // Check if chrome has been processed
            if(chrDone[idx] != 0) {
                fprintf(stderr, "%s has already been processed. The bam file is not sorted.\n", curChr);
                exit(1);
            }
            chrDone[idx] = 1;
        }

        /////////////////////////////////////////////////
        // Filter Step 1: Base quality
        /////////////////////////////////////////////////
        for(i = 0; i < record->len; i++) {
            if(record->seq[i] != 'N' && (unsigned short)(record->qual[i] - 33) < minquali)
                record->seq[i] = 'N';
        }
        /////////////////////////////////////////////////
        // Filter Step 2: SNP per base
        /////////////////////////////////////////////////
        //printf("\n%s %d %s %d\n", record->qname, record->len, record->seq, record->offset);
        off = record->offset - 1;
        iread = (int)(record->len * errorrate + 0.5);
        if(iread < 1)
            iread = 1;
        cnt = 0;
        for(i = 0; i < record->len; i++) {
            if(record->strand == '+') {
                switch(record->seq[i]) {
                    case 'A':
                        if(chrSeqArray[idx][off] != 'A' && (!(record->r12 == 2 && chrSeqArray[idx][off] == 'G')))
                            cnt++;
                        break;
                    case 'T':
                        if(chrSeqArray[idx][off] != 'T' && (!(record->r12 == 1 && chrSeqArray[idx][off] == 'C')))
                            cnt++;
                        break;
                    case 'C':
                        if(chrSeqArray[idx][off] != 'C')
                            cnt++;
                        break;
                    case 'G':
                        if(chrSeqArray[idx][off] != 'G')
                            cnt++;
                        break;
                }
            }
            else {
                switch(record->seq[i]) {
                    case 'A':
                        if(chrSeqArray[idx][off] != 'A' && (!(record->r12 == 1 && chrSeqArray[idx][off] == 'G')))
                            cnt++;
                        break;
                    case 'T':
                        if(chrSeqArray[idx][off] != 'T' && (!(record->r12 == 2 && chrSeqArray[idx][off] == 'C')))
                            cnt++;
                        break;
                    case 'C':
                        if(chrSeqArray[idx][off] != 'C')
                            cnt++;
                        break;
                    case 'G':
                        if(chrSeqArray[idx][off] != 'G')
                            cnt++;
                        break;
                }
            }

            off++;
        }
        // Snp per read
        //printf("\n%d < %d %s %f %s %d\n", cnt, iread, record->qname, errorrate, record->qual, record->offset - 1);
        if(cnt <= iread) {
            off = record->offset - 1;
            for(i = 0; i < record->len; i++) {
                if(record->strand == '+') {
                    switch(record->seq[i]) {
                        case 'A':
                            if(0 && record->r12 ==2 && chrSeqArray[idx][off] == 'G') {
                                c_G[off]++;
                                c_Gq[off] += (unsigned short)(record->qual[i] - 33);
                                c_Q[off] += record->mapq;
                            }
                            else {
                                if(record->r12 == 1) {
                                    w_A[off]++;
                                    w_Aq[off] += (unsigned short)(record->qual[i] - 33);
                                    w_Q[off] += record->mapq;
                                }
                                else {
                                    c_A[off]++;
                                    c_Aq[off] += (unsigned short)(record->qual[i] - 33);
                                    c_Q[off] += record->mapq;
                                }
                            }
                            break;
                        case 'T':
                            if(0 && record->r12 == 1 && chrSeqArray[idx][off] == 'C') { //0&&for methylation
                                w_C[off]++;
                                w_Cq[off] += (unsigned short)(record->qual[i] - 33);
                                w_Q[off] += record->mapq;
                            }
                            else {
                                if(record->r12 == 1) {
                                    w_T[off]++;
                                    w_Tq[off] += (unsigned short)(record->qual[i] - 33);
                                    w_Q[off] += record->mapq;
                                 }
                                else {
                                    c_T[off]++;
                                    c_Tq[off] += (unsigned short)(record->qual[i] - 33);
                                    c_Q[off] += record->mapq;
                                }
                            }
                            break;
                        case 'C':
                            if(record->r12 == 1) {
                                w_C[off]++;
                                w_Cq[off] += (unsigned short)(record->qual[i] - 33);
                                w_Q[off] += record->mapq;
                            }
                            else {
                                c_C[off]++;
                                c_Cq[off] += (unsigned short)(record->qual[i] - 33);
                                c_Q[off] += record->mapq;
                            }
                            break;
                        case 'G':
                            if(record->r12 == 1) {
                                w_G[off]++;
                                w_Gq[off] += (unsigned short)(record->qual[i] - 33);
                                w_Q[off] += record->mapq;
                            }
                            else {
                                c_G[off]++;
                                c_Gq[off] += (unsigned short)(record->qual[i] - 33);
                                c_Q[off] += record->mapq;
                            }
                            break;
                    }
                }
                else {
                    switch(record->seq[i]) {
                        case 'A':
                            if(0 && record->r12 == 1 && chrSeqArray[idx][off] == 'G') {
                                c_G[off]++;
                                c_Gq[off] += (unsigned short)(record->qual[i] - 33);
                                c_Q[off] += record->mapq;
                            }
                            else {
                                if(record->r12 == 2) {
                                    w_A[off]++;
                                    w_Aq[off] += (unsigned short)(record->qual[i] - 33);
                                    w_Q[off] += record->mapq;
                                }
                                else {
                                    c_A[off]++;
                                    c_Aq[off] += (unsigned short)(record->qual[i] - 33);
                                    c_Q[off] += record->mapq;
                                }
                            }
                            break;
                        case 'T':
                            if(0 && record->r12 == 2 && chrSeqArray[idx][off] == 'C') {
                                w_C[off]++;
                                w_Cq[off] += (unsigned short)(record->qual[i] - 33);
                                w_Q[off] += record->mapq;
                            }
                            else {
                                if(record->r12 == 2) {
                                    w_T[off]++;
                                    w_Tq[off] += (unsigned short)(record->qual[i] - 33);
                                    w_Q[off] += record->mapq;
                                }
                                else {
                                    c_T[off]++;
                                    c_Tq[off] += (unsigned short)(record->qual[i] - 33);
                                    c_Q[off] += record->mapq;
                                }
                            }
                            break;
                        case 'C':
                            if(record->r12 == 2) {
                                w_C[off]++;
                                w_Cq[off] += (unsigned short)(record->qual[i] - 33);
                                w_Q[off] += record->mapq;
                            }
                            else {
                                c_C[off]++;
                                c_Cq[off] += (unsigned short)(record->qual[i] - 33);
                                c_Q[off] += record->mapq;
                            }
                            break;
                        case 'G':
                            if(record->r12 == 2) {
                                w_G[off]++;
                                w_Gq[off] += (unsigned short)(record->qual[i] - 33);
                                w_Q[off] += record->mapq;
                            }
                            else {
                                c_G[off]++;
                                c_Gq[off] += (unsigned short)(record->qual[i] - 33);
                                c_Q[off] += record->mapq;
                            }
                            break;
                    }
                }

                off++;
            }
        }
        else {
            #ifdef __MY_DEBUG__
            if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && record->offset <= __DEBUG_POS__ && record->offset+record->len >= __DEBUG_POS__) {
                fprintf(stderr, "read discarded because snp-ratio overflow. [threshold %d, actually %d.]\n", iread, cnt);
                #ifdef __MY_DEBUG_READ__
                dispRecord(record);
                #endif
            }
            #endif
        }
        // Update row counter
        rowCnt++;
        if(rowCnt % 100000 == 0)
            fprintf(stderr, "%d records have been dealed.\n", rowCnt);
    }

    // Last batch
    if(finCnt < rowCnt) {
        // Print
        printSnp(methFptr, snpFptr, chrSeqArray, idx, len, minhetfreq, curChr, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, w_Aq, w_Tq, w_Cq, w_Gq, c_Aq, c_Tq, c_Cq, c_Gq, w_Q, c_Q);
        // Memory gathering for x_X
        fprintf(stderr, "Memory free for SNP/meth, array len: %u\n", longestchr);
        free(w_A);
        free(w_T);
        free(w_C);
        free(w_G);
        free(c_A);
        free(c_T);
        free(c_C);
        free(c_G);
        // Memory gathering for x_Xq
        free(w_Aq);
        free(w_Tq);
        free(w_Cq);
        free(w_Gq);
        free(c_Aq);
        free(c_Tq);
        free(c_Cq);
        free(c_Gq);

        // Memory gathering for x_Q
        free(w_Q);
        free(c_Q);
    }

    bam_header_destroy(header);
    bam_close(in);
    bam_destroy1(b);
}

