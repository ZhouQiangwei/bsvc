CC = g++
FLAGS = -O2 
LDFLAGS = -O2
HTSLIB=submodules/htslib
SAMTOOLS=submodules/samtools
LIBS = -m64 -I./$(SAMTOOLS) -I$(HTSLIB) -L$(SAMTOOLS) -L$(HTSLIB) -lbam -lz -lpthread -lm -lgsl -lgslcblas -lhts

TARGET = bsnps
TOCOMPILE = bsnps.o bamprocess.o hashtable.o readgenome.o fisher.o bayes.o getchrLen.o genotype.o

all: ${TOCOMPILE}
	${CC} $(LDFLAGS) -o $(TARGET) ${TOCOMPILE} ${LIBS}

.c.o:
	$(CC) ${FLAGS} ${LIBS} -c $*.c

clean:
	rm -f *.o $(TARGET)

