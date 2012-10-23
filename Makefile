SHELL=bash
.PHONY: gitstatus.txt help clean
BITS=
#BITS=-m32
#BITS=-m64

MAIN=vcsbm

all: ${MAIN}

clean:
	-rm ${MAIN} *.o */*.o


CXXFLAGS=       \
          -std=gnu++0x   \
          -Wmissing-field-initializers   \
          -Wsign-compare   \
          -Wuninitialized   \
          -Wunused-parameter    \
          -Wunused             \
          -Wnon-virtual-dtor            \
          -Wall -Wformat -Werror -Wextra #-Wconversion # scf.cpp doesn't like -Wconversion

CXX=g++
CC=g++
#CXXFLAGS= ${BITS}     -g
LDFLAGS+= -lrt
LDFLAGS+= `gsl-config --libs`
CXXFLAGS:= ${BITS} -O3        ${CXXFLAGS} -std=gnu++0x # -DNDEBUG
#CXXFLAGS+= -p -pg
#CXXFLAGS=              -O2                 

#${MAIN}: CXXFLAGS += -DNDEBUG
${MAIN}: gitstatus.o cmdline.o ${MAIN}.o
	${CXX} ${CXXFLAGS} $^ ${LDFLAGS} -o $@

#lineGraph: lineGraph.o shmGraphRaw.o Range.o
gitstatus.txt: 
	@echo '"\n"' > gitstatus.txt
	-@type git > /dev/null 2>&1 && { git log | head -n 1 ; git status ; } | head -n 20 | sed -re 's/"/\\"/g ; s/^/"/g; s/$$/\\n"/g; ' > gitstatus.txt

gitstatus.o: comment.txt  gitstatus.txt

cmdline.c.FORCE:
	# remake cmdline.c . But it's OK unless you change the .ggo file. You'll need gengetopt(1) to be able to run this.
	gengetopt  --unamed-opts < cmdline.ggo
