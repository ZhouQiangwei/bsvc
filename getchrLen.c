#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define BUFSIZE 1000000

int getchrLen(char* refSeqFile)
{
	FILE* iFPtr;
	FILE* oFPtr;
	char chrLenFile[500];
	
	strcpy(chrLenFile, refSeqFile);
	strcat(chrLenFile, ".len");
	
	fprintf(stderr, "Process genome and get refseq file: %s.\n", refSeqFile);
	//fprintf(stderr, "chrLenFile = %s.\n", chrLenFile);
	
	char* readBuffer = (char*)malloc(sizeof(char) * BUFSIZE);
	int len = 0;
	int lines = 0;
	
	char* token;
	char seps[] = " \t\n\r";
	char chrName[100];
	
	if(!(iFPtr = fopen(refSeqFile, "r"))) {
		printf("File %s open error.\n", refSeqFile);
		exit(1);
	}
	
	if(!(oFPtr = fopen(chrLenFile, "w"))) {
		printf("File %s open error.\n", chrLenFile);
		exit(1);
	}
	
	while(fgets(readBuffer, BUFSIZE, iFPtr)) {
		if(strlen(readBuffer) >= BUFSIZE - 1) {
            fprintf(stderr, "Too many characters in one row! Try to split the long row into several short rows (fewer than %d characters per row).\n", BUFSIZE);
            exit(1);
        }
		
		if(readBuffer[0] == '>') {
			if(lines > 0) {
				fprintf(oFPtr, "%s\t%d\n", chrName, len);
			}
			// Save name
			token = strtok(readBuffer + 1, seps);	
			strcpy(chrName, token);
			len = 0;	
		}
		else {
			// Substract \n
			len += strlen(readBuffer) - 1;
		}
		
		lines++;
	}
	
	if(lines > 0) {
		fprintf(oFPtr, "%s\t%d\n", chrName, len);
	}
	
	fclose(iFPtr);
	fclose(oFPtr);
	
	return 0;
}
