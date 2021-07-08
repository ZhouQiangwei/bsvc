/*
* Created on Fri Sep 21 10:29:50 2018
* @author: qwzhou
* @email: qiangwei.zhou2013@gmail.com
*/
#include <gsl/gsl_cdf.h>
#include <tr1/cmath>
#include <gsl/gsl_sf_gamma.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <string>
#include <string.h>
#include <math.h>
#include <iostream>
#include <pthread.h>

#include "hashtable.h"
#include "readgenome.h"
#include "bamprocess.h"
#include "getchrLen.h"

int multiout=0;
unsigned int longestchr = 10000;
int minquali =20;
int mincover=3; //10;
int minread2=2;
int maxcover=1000;
float minhetfreq=0.1;
float errorrate=0.02; //0.1
unsigned int mapqThr=20;
int meth=0;
int snp=1;
float minhomfreq =0.85;
float pvalue_cutoff =0.05;
int NTHREAD=0;
char bamFileName[1024];
char** chrName;
int othercover=80;
float minvarrate=0.3;
int minvarread=2;
int methmincover=4;
int printLowQ=0;
float tvalcutoff=0.01;

struct Threading
{
    pthread_t Thread;
    unsigned r;
    void *ret;
    ARGS Arg;
};

int main(int argc, char* argv[])
{
    time_t Start_Time,End_Time;
    const char* Help_String="Command Format :  BSsnpScan [options] -g GENOME -i <Bamfile> -o <SNP outfile>\n"
        "\nUsage:\n"
        "\t-g|--genome           Genome\n"
        "\t-i|--input            Sorted bam format file\n"
        "\t-o|--output           SNP output file\n"
        "\t-p|--threads          the number of threads, must have .bai file by samtools index in the bam folder.\n"
        "\t--minQ                Minimum mapping quality score to count a read, default: 20\n"
        "\t--minquali            Minimum base quality at a position to count a read, default: 15\n"
        "\t--errorrate           Minimum unmatch base allowed to count a read, default: 0.1\n"
        "\t--minread2            Minimum supporting reads at a position to call variants, default: 2\n"
        "\t--mincover            Minimum read depth at a position to make a call, default: 10\n"
        "\t--maxcover            Maximum read depth at a position to make a call, default: 1000\n"
        "\t--othercover          Maximum read depth at a position not the reference and variants, default: 80\n"
        "\t--pvalue              Default p-value threshold for calling variants, default: 0.05\n"
        "\t--minhetfreq          Minimum hetero frequency threshold, default: 0.1\n"
        "\t--minhomfreq          Minimum homo frequency threshold, default: 0.85\n"
        "\t--tvalue              Default pvalue threshold for T test in calling variants, default: 0.01\n"
        "\t--printLowQ           Print Low Quality SNP, default: 0. [0 or 1]\n"
        //"\t-m                  Report DNA methylation calling positions. default: No report\n"
        "\t-mc [filename]        Report DNA methylation calling positions, and DNA methylation file.\n"
        "\t--methmincover        DNA methylation minimum read depth at a position to make a call, default: 5\n"
        "\t--multiout            Output results write to one file, or multi files with {out.chrom} prefix. only useful when number of threads bigger than 1. [0 or 1]\n"
        //"\t-mchg                 DNA methylation CHG file.\n"
        //"\t-mchh                 DNA methylation CHH file.\n"
        "\t-h|--help             BSNPS usage.\n";
        //"\t------------          If variant position filter quality is Low, we also can set `AD>minvarread, ALFR>minvarrate, pvalue<pvalue_cutoff+0.01` as PASS.\n"
        //"\t--minvarread          default: 5\n"
        //"\t--minvarrate          default: 0.3";
        

    char refSeqFile[1024];
    char snpFileName[1024];
    char methFileName[1024];
    //char methChgFileName[1024];
    //char methChhFileName[1024];
    for(int i=1;i<argc;i++)
    {
        if(!strcmp(argv[i], "-g") ||!strcmp(argv[i], "--genome")  )
        {
            strcpy(refSeqFile, argv[++i]);
        }else if(!strcmp(argv[i], "-i") ||!strcmp(argv[i], "--input")  )
        {
            strcpy(bamFileName, argv[++i]);
        }else if(!strcmp(argv[i], "-o") ||!strcmp(argv[i], "--output")  )
        {
            strcpy(snpFileName, argv[++i]);
        }else if(!strcmp(argv[i], "-mc") )
        {
            meth=1;
            strcpy(methFileName, argv[++i]);
        }/*else if(!strcmp(argv[i], "-mchg") )
        {
            strcpy(methChgFileName, argv[++i]);
        }else if(!strcmp(argv[i], "-mchh") )
        {
            strcpy(methChhFileName, argv[++i]);
        }*/
        else if(!strcmp(argv[i], "--minquali") )
        {
            minquali=atoi(argv[++i]);
        }else if(!strcmp(argv[i], "--tvalue") )
        {
            tvalcutoff=atof(argv[++i]);
        }else if(!strcmp(argv[i], "--printLowQ") )
        {
            printLowQ=atoi(argv[++i]);
        }else if(!strcmp(argv[i], "--minQ") )
        {
            mapqThr=atoi(argv[++i]);
        }else if(!strcmp(argv[i], "--methmincover") )
        {
            methmincover=atoi(argv[++i])-1;
            if(methmincover<0) {
                fprintf(stderr, "\nmethmincover must >0\n");
                exit(0);
            }
        }else if(!strcmp(argv[i], "--minread2") )
        {
            minread2=atoi(argv[++i]);
        }else if(!strcmp(argv[i], "--errorrate") )
        {
            errorrate=atof(argv[++i]);
        }else if(!strcmp(argv[i], "--mincover") )
        {
            mincover=atoi(argv[++i]);
        }else if(!strcmp(argv[i], "--maxcover") )
        {
            maxcover=atoi(argv[++i]);
        }else if(!strcmp(argv[i], "--othercover") )
        {
            othercover=atoi(argv[++i]);
        }else if(!strcmp(argv[i], "--pvalue") )
        {
            pvalue_cutoff=atof(argv[++i]);
        }else if(!strcmp(argv[i], "--minvarread") )
        {
            minvarread=atoi(argv[++i]);
        }else if(!strcmp(argv[i], "--minvarrate") )
        {
            minvarrate=atof(argv[++i]);
        }else if(!strcmp(argv[i], "--minhetfreq") )
        {
            minhetfreq=atof(argv[++i]);
        }else if(!strcmp(argv[i], "--minhomfreq") )
        {
            minhomfreq=atof(argv[++i]);
        }else if(!strcmp(argv[i], "-p") ||!strcmp(argv[i], "--threads")  )
        {
            NTHREAD=atoi(argv[++i]);
        }else{
            fprintf(stderr, "There is no %s paramater.\n", argv[i]);
            exit(0);
        }

    }
    if(argc<3) {
        fprintf(stderr, "%s\n", Help_String);
        exit(0);
    }
    char* chrLenFile=(char*) malloc(sizeof(char)*1024);
    sprintf(chrLenFile, "%s.len", refSeqFile);
    
    fprintf(stderr, "\nBisulfite-seq Single Nucleotide Polymorphism Scan (BSNPS)\n--------------------------------------------------------------\n");

    if ( access(chrLenFile,0) ){
        fprintf(stderr, "build genome len file ...\n");
        getchrLen(refSeqFile);
    }

    int i=0; int j=0;
    long cnt;
    FILE *bamFptr, *posFptr;
    ARGS args;
    fprintf(stderr, "Chrom Len: %s\n", chrLenFile);
    fprintf(stderr, "Genome File: %s\n", refSeqFile);
    fprintf(stderr, "BAM File: %s\n", bamFileName);
    //args.bamFileName=bamFileName;

    fprintf(stderr, "SNP File: %s\n", snpFileName);
    if(meth){
        fprintf(stderr, "Meth File: %s\n", methFileName);
        fprintf(stderr, "Meth coverage: %d\n", methmincover);
        //fprintf(stderr, "MethChg File: %s\n", methChgFileName);
        //fprintf(stderr, "MethChh File: %s\n", methChhFileName);
        //args.methFileName = methFileName;
        //args.methChgFileName = methChgFileName;
        //args.methChhFileName = methChhFileName;
    }
    fprintf(stderr, "Minimum base quality: %d\n", minquali);
    fprintf(stderr, "Minimum supporting reads at a position to call variants: %d\n", minread2);
    fprintf(stderr, "Mininum read depth: %d\n", mincover);
    fprintf(stderr, "Maximum read depth: %d\n", maxcover);
    fprintf(stderr, "Minimum hetero frequency: %f\n", minhetfreq);
    fprintf(stderr, "Minimum mapping quality score: %d\n", mapqThr);
    fprintf(stderr, "Pvalue threshold for calling variants: %f\n", pvalue_cutoff);
    fprintf(stderr, "Errorate per read: %f\n", errorrate);
    fprintf(stderr, "minvarread: %d;\n", minvarread);
    fprintf(stderr, "T test pvalue threshold for calling variants: %f\n", tvalcutoff);
    
    //////////////////////////////////////////////////////////////////////////////
    // Load reference sequence & Init arrays
    //////////////////////////////////////////////////////////////////////////////
    // {{
    // Hash table
    int hash_table_size;
    args.hashTable = (HashNode**)malloc(sizeof(HashNode*) * HASH_TABLE_MAX_SIZE);
    hash_table_init(args.hashTable, &hash_table_size);
    // Init chrome name-idx hash table
    //int chrCnt;
    init_chrome_hash(args.hashTable, &hash_table_size, chrLenFile, &args.chrCnt);
    fprintf(stderr, "Init chrome name-idx hash table completed\n");
    // Init chrome name and length array
    args.chrLen = (int*)malloc(sizeof(int) * args.chrCnt);
    chrName = (char**)malloc(sizeof(char*) * args.chrCnt);
    for(i = 0; i < args.chrCnt; i++)
        chrName[i] = (char*)malloc(sizeof(char) * 50);
    init_chrome_name_len(args.hashTable, chrLenFile, args.chrCnt, chrName, args.chrLen);
    cnt = 0;
    for(i = 0; i < args.chrCnt; i++){
        if(args.chrLen[i]>longestchr) longestchr=args.chrLen[i];
        cnt += args.chrLen[i];
    }
    fprintf(stderr, "Longest chrom %d bp, total %ld bp of %d chromosomes\n", longestchr, cnt, args.chrCnt);
    // Init chrome seq array
    args.chrSeqArray = (char**)malloc(sizeof(char*) * args.chrCnt);
    for(i = 0; i < args.chrCnt; i++)
        args.chrSeqArray[i] = (char*)malloc(sizeof(char) * args.chrLen[i] + 100);
    init_chrome_seq(args.hashTable, refSeqFile, args.chrSeqArray, args.chrLen);
        for(i = 0; i < args.chrCnt; i++) {
        if(!(args.chrSeqArray[i] = (char*)realloc((void*)(args.chrSeqArray[i]), sizeof(char) * args.chrLen[i]))) {
            fprintf(stderr, "Not enough memory for chrSeqArray\n");
            exit(1);
        }
        for(j = 0; j < args.chrLen[i]; j++)
            args.chrSeqArray[i][j] = toupper(args.chrSeqArray[i][j]);
    }
    fprintf(stderr, "Build chrome seq array completed\n");

    int* chrDone = (int*)malloc(sizeof(int) * args.chrCnt);
    chrDone[i] = 0;
    // }}

    //Meth File
    if(multiout) args.methFileName = methFileName;
    if(meth == 1){
        args.methFptr = fopen(methFileName, "w");
        fprintf(args.methFptr, "chrom\tpos\tstrand\tcontext\tmethcover\t(meth+unmeth)cover\tmethylevel\tmethqual,unmethqual\twaston_cover,crick_cover\trefbase,genotype\n");
    }
    //////////////////////////////////////////////////////////////////////////////
    // SNP process
    //////////////////////////////////////////////////////////////////////////////
    fprintf(stderr, "SNP process ...\n");
    if(multiout) args.snpFileName = snpFileName;
    args.snpFptr = fopen(snpFileName, "w");
    
    //////////////////////////////////////////////////////////////////////////////
    // vcf header
    const char* vcfheader="##fileformat=VCFv4.1\n"
    "##source=BSNPS\n"
    "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">"
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">\n"
    "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Total allelic depths\">\n"
    "##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Total allelic depths on the forward strand\">\n"
    "##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Total allelic depths on the reverse strand\">\n"
    "##FILTER=<ID=Low,Description=\"Low Quality\">\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes for each ALT allele, in the same order as listed\">\n"
    "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n"
    "##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Average mapping quality\">\n"
    "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">\n"
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Quality Read Depth of bases with high-quality bases\">\n"
    "##FORMAT=<ID=PVAL,Number=1,Type=String,Description=\"P-value from Fisher's Exact Test\">\n"
    "##FORMAT=<ID=SP,Number=1,Type=Integer,Description=\"Phred-scaled strand bias P-value\">\n"
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Depth of variant-supporting bases\">\n"
    "##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths on the forward strand\">\n"
    "##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths on the reverse strand\">\n"
    "##FORMAT=<ID=BSD,Number=8,Type=Integer,Description=\"Depth of A,T,C,G in watson strand and crick strand\">\n"
    "##FORMAT=<ID=BSQ,Number=8,Type=Integer,Description=\"Avarage Base Quality of A,T,C,G in watson strand and crick strand\">\n"
    "##FORMAT=<ID=ALFR,Number=R,Type=Float,Description=\"Allele frequency\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1";

    //##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">
    //##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">
    //##INFO=<ID=IDV,Number=1,Type=Integer,Description=\"Maximum number of reads supporting an indel\">
    //##INFO=<ID=IMF,Number=1,Type=Float,Description=\"Maximum fraction of reads supporting an indel\">
    //////////////////////////////////////////////////////////////////////////////

    fprintf(args.snpFptr, "%s\n", vcfheader);

    if (NTHREAD)
    {
        Threading* Thread_Info=(Threading*) malloc(sizeof(Threading)*NTHREAD);
        pthread_attr_t Attrib;
        pthread_attr_init(&Attrib);
        pthread_attr_setdetachstate(&Attrib, PTHREAD_CREATE_JOINABLE);

        for (int i=0;i<NTHREAD;i++)
        {
            args.ThreadID=i;
            Thread_Info[i].Arg=args;
            Thread_Info[i].r=pthread_create(&Thread_Info[i].Thread,&Attrib,npsnpAnalysis,(void*) &Thread_Info[i].Arg);
            if(Thread_Info[i].r) {printf("Launch_Threads():Cannot create thread..\n");exit(-1);}
        }
        pthread_attr_destroy(&Attrib);
        for (int i=0;i<NTHREAD;i++)
        {
            pthread_join(Thread_Info[i].Thread,NULL);
        }
        free(Thread_Info);
    }else{
        args.ThreadID=0;
        npsnpAnalysis(&args);
    }
    /*   
    if(multiout){
        char tempoutfile[1000];
        for(int i=0;i<chrCnt;i++){
            sprintf(tempoutfile, "%s.%s.vcf", args->snpFileName, chrName[processchrom]);
            FILE* snptempfp = fopen(tempoutfile, "r");

            fclose(snptempfp);
        }
    }
    */
    fclose(args.snpFptr);
    if(meth==1) fclose(args.methFptr);

    fprintf(stderr, "SNP process done!\n");

    return 0;
}


