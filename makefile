CC=g++

FLAGS=-g -O2 -fopenmp
CFLAGS=-g -O2 -fopenmp

#FLAGS=-g -fopenmp
#CFLAGS=-g -fopenmp

SAMTOOLS_ROOT=/home/jry/zhaodan_practice/samtools-0.1.19/
MLPACK_LIB=/home/jry/workspace/mlpack-3.0.0/build/lib/
MLPACK_INCLUDE=/home/jry/workspace/mlpack-3.0.0/build/include/
FLAGS+=-I${SAMTOOLS_ROOT}
FLAGS+=-I${MLPACK_INCLUDE}
LFLAGS=-lm -L${SAMTOOLS_ROOT} -L${MLPACK_LIB} -lbam -lz -lpthread -lmlpack -lboost_program_options -larmadillo -g
SOURCE = cmds scan distribution refseq polyscan param utilities homo window bamreader sample chi somatic
OBJS= $(patsubst %,%.o,$(SOURCE))

all: check-samtools msisensor

%.o:%.cpp
	$(CC) $(FLAGS) -g -c $< -o $@

check-samtools:
    ifndef SAMTOOLS_ROOT
        $(error SAMTOOLS_ROOT is undefined)
    endif

msisensor: $(OBJS)
	$(CC) $^ $(CFLAGS) $(LFLAGS) -o $@ 
clean:
	rm -f *.o msisensor


