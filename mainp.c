#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>


#include "hash_funcs.h"
#include "chrome_funcs.h"
#include "sam_funcs.h"

int vQualMin =20;
int nLayerMin=2;
int nLayerMax=1000;
float vSnpRate=0.1;
float vSnpPerBase=0.4;
unsigned int mapqThr=20;
int meth=0;
int snp=1;

int main(int argc, char **argv)
{
	int i, j;
	long cnt;
	FILE *bamFptr, *posFptr;
	//////////////////////////////////////////////////////////////////////////////
	// Parse command args
	//////////////////////////////////////////////////////////////////////////////
	// {{
	if(argc != 14) {
		printf("Not enough args! Argc = %d.\n", argc);
		exit(1);
	}
        ARGS args;
        THREAD=0;

	char* chrLenFile = argv[1];
	char* refSeqFile = argv[2];
	args.bamFileName = argv[3];
	args.snpFileName = argv[4];
	args.methCgFileName = argv[5];
        args.methChgFileName = argv[6];
        args.methChhFileName = argv[7];

        vQualMin = atoi(argv[8]);
	nLayerMin = atoi(argv[9]);
	nLayerMax = atoi(argv[10]);
	vSnpRate = atof(argv[11]);
	vSnpPerBase = atof(argv[12]);
	mapqThr = atoi(argv[13]);
	
	fprintf(stderr, "chrLenFile = %s.\n", chrLenFile);
	fprintf(stderr, "refSeqFile = %s.\n", refSeqFile);
	fprintf(stderr, "bamFileName = %s.\n", bamFileName);
	fprintf(stderr, "snpFileName = %s.\n", snpFileName);
	fprintf(stderr, "methCgFileName = %s.\n", methCgFileName);
    fprintf(stderr, "methChgFileName = %s.\n", methChgFileName);
    fprintf(stderr, "methChhFileName = %s.\n", methChhFileName);

	fprintf(stderr, "vQualMin = %d.\n", vQualMin);
	fprintf(stderr, "nLayerMax = %d.\n", nLayerMax);
	fprintf(stderr, "vSnpRate = %f.\n", vSnpRate);
	fprintf(stderr, "vSnpPerBase = %f.\n", vSnpPerBase);
	fprintf(stderr, "mapqThr = %d.\n", mapqThr);
	// }}

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
	fprintf(stderr, "Init chrome name-idx hash table completed.\n");
	// Init chrome name and length array
	args.chrLen = (int*)malloc(sizeof(int) * args.chrCnt);
	char** chrName = (char**)malloc(sizeof(char*) * args.chrCnt);
	for(i = 0; i < args.chrCnt; i++)
		chrName[i] = (char*)malloc(sizeof(char) * 50);
	init_chrome_name_len(args.hashTable, chrLenFile, args.chrCnt, chrName, args.chrLen);
	cnt = 0;
	for(i = 0; i < args.chrCnt; i++)
		cnt += args.chrLen[i];
	fprintf(stderr, "Init chrome name-len array completed, total %ld bp of %d chromosomes.\n", cnt, chrCnt);
	// Init chrome seq array
	args.chrSeqArray = (char**)malloc(sizeof(char*) * args.chrCnt);
	for(i = 0; i < args.chrCnt; i++)
		args.chrSeqArray[i] = (char*)malloc(sizeof(char) * args.chrLen[i] + 100);
	init_chrome_seq(args.hashTable, refSeqFile, args.chrSeqArray, args.chrLen);
        for(i = 0; i < args.chrCnt; i++) {
		if(!(args.chrSeqArray[i] = (char*)realloc((void*)(args.chrSeqArray[i]), sizeof(char) * args.chrLen[i]))) {
			fprintf(stderr, "Not enough memory for chrSeqArray.\n");
			exit(1);
		}
		for(j = 0; j < args.chrLen[i]; j++)
			args.chrSeqArray[i][j] = toupper(args.chrSeqArray[i][j]);
	}
	fprintf(stderr, "Init chrome seq array completed.\n");

	int* chrDone = (int*)malloc(sizeof(int) * args.chrCnt);
	chrDone[i] = 0;
	// }}

	//////////////////////////////////////////////////////////////////////////////
	// SNP process
	//////////////////////////////////////////////////////////////////////////////
	Thread_Arg T;
        if (THREAD)
        {
            Threading* Thread_Info=(Threading*) malloc(sizeof(Threading)*NTHREAD);
            pthread_attr_t Attrib;
            pthread_attr_init(&Attrib);
            pthread_attr_setdetachstate(&Attrib, PTHREAD_CREATE_JOINABLE);
            for (int i=0;i<NTHREAD;i++)
            {
                args.ThreadID=i;
                Thread_Info[i].r=pthread_create(&Thread_Info[i].Thread,&Attrib,npsnpAnalysis,(void*) &args);
                if(Thread_Info[i].r) {printf("Launch_Threads():Cannot create thread..\n");exit(-1);}
            }
            pthread_attr_destroy(&Attrib);
            for (int i=0;i<NTHREAD;i++)
            {
                pthread_join(Thread_Info[i].Thread,NULL);
            }
            free(Thread_Info);

//snpAnalysis(bamFileName, snpFileName, methCgFileName, methChgFileName, methChhFileName, hashTable, chrSeqArray, chrLen, chrCnt);
        }else
            npsnpAnalysis(&args);
      //   snpAnalysis(bamFileName, snpFileName, methCgFileName, methChgFileName, methChhFileName, hashTable, chrSeqArray, chrLen, chrCnt, vQualMin, nLayerMin, nLayerMax, vSnpRate, vSnpPerBase, mapqThr, meth, snp);
	fprintf(stderr, "SNP process completed.\n");

	return 0;
}
