CC=clang++
INCFLAG=$(shell rivet-config --cppflags)
LINKFLAG=$(shell rivet-config --ldflags)
LIBS=-lRivet -lHepMC -lfastjet

all: test.cc
	${CC} -o test test.cc ${INCFLAG} ${LINKFLAG} ${LIBS}
