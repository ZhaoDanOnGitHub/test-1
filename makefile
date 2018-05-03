CC=g++

FLAGS=-g -O2 -fopenmp
CFLAGS=-g -O2 -fopenmp

#FLAGS=-g -fopenmp
#CFLAGS=-g -fopenmp

SAMTOOLS_ROOT=/home/zhaodan/workspace/samtools-0.1.19/
FLAGS+=-I${SAMTOOLS_ROOT}
LFLAGS=-lm -L${SAMTOOLS_ROOT} -lbam -lz -lpthread
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


