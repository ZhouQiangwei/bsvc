CC = g++
FLAGS = -O2 
LDFLAGS = -O2
HTSLIB=submodules/htslib
SAMTOOLS=./  #submodules/samtools
BAMFLAG= #-lbam
LIBS = -m64 -I$(SAMTOOLS) -I$(HTSLIB) -L$(SAMTOOLS) -L$(HTSLIB) $(BAMFLAG) -lz -lpthread -lm -lgsl -lgslcblas -lhts

TARGET = bsvc
TOCOMPILE = bsnps.o bamprocess.o hashtable.o readgenome.o fisher.o bayes.o getchrLen.o genotype.o tstudenttest.o

all: ${TOCOMPILE}
	make -C $(HTSLIB) lib-static
	#make clean-so -C $(HTSLIB)
	#make -C $(SAMTOOLS) lib
	${CC} $(LDFLAGS) -o $(TARGET) ${TOCOMPILE} ${LIBS}

.c.o:
	$(CC) ${FLAGS} ${LIBS} -c $*.c

deepclean: clean
	make clean -C $(HTSLIB)
	#make clean -C $(SAMTOOLS)

clean:
	rm -f *.o $(TARGET)

