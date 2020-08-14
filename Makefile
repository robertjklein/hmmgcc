# Makefile for hmmgcc
# Modified from SRE's Makefile for SQUID

## Compiler
## Only specifiy gcc if different from cc (i.e. not on Linux)
CC = cc
#CC = gcc

## Where to install
#
INSTALLDIR = /usr/local/bin

## Compiler flags 
#
CFLAGS =   
#CFLAGS = -g # GCC debugging
#CFLAGS = -O2 # GCC optimization 

#######
## should not need to modify below this line
#######
LIBS      = -lm


SRC =   hmm.c viterbi.c hmmgcc.c

OBJ =   hmm.o viterbi.o hmmgcc.o

all: hmmgcc
hmmgcc: $(OBJ)
	$(CC) $(CFLAGS) -o hmmgcc $(OBJ) $(LIBS)

install: hmmgcc
	test -d $(INSTALLDIR) || mkdir -p $(INSTALLDIR)
	cp hmmgcc $(INSTALLDIR)

depend:
	makedepend $(SRC)

clean:
	-rm -f *.o *~ core TAGS  
clobber:
	-rm -f *.o *~ core TAGS $(PROGS)


# DO NOT DELETE

