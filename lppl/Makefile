CC = gcc

CFLAGS = -g -DDEBUG -Wall -fopenmp
# CFLAGS = -O3 -Wall -ftree-vectorize -funroll-loops -fopenmp
# CFLAGS = -O3 -DNDEBUG -Wall -ftree-vectorize -funroll-loops -fopenmp

INC = -I . -I /home/vincenzo/Freeware/levmar-2.4 -I ../shared
OBJ = args.o print.o parms.o func-lppl.o levmar.o guess.o weight.o \
      iterations.o noise.o \
      ../shared/smalloc.o ../shared/debug.o ../shared/slock.o \
      ../shared/error.o ../shared/vector.o ../shared/smmap.o ../shared/sfork.o \
      ../shared/ssignal.o ../shared/time.o ../shared/random.o \
      ../shared/minmax.o \
      read.o correlation.o
LIBDIR = -L /home/vincenzo/Freeware/levmar-2.4
LIBS = -lm -lpthread -lpopt -llapack -lblas -llevmar

.PHONY: all
all:	lppl explore-lppl guess-lppl gen-lppl

lppl:	lppl.c $(OBJ) $(LMAR)
	$(CC) $(INC) $(CFLAGS) -o $@ $(OBJ) lppl.c $(LIBDIR) $(LIBS)

explore-lppl:	explore-lppl.c $(OBJ) $(LMAR)
	$(CC) $(INC) $(CFLAGS) -o $@ $(OBJ) explore-lppl.c $(LIBDIR) $(LIBS)

guess-lppl:	guess-lppl.c $(OBJ) $(LMAR)
	$(CC) $(INC) $(CFLAGS) -o $@ $(OBJ) guess-lppl.c $(LIBDIR) $(LIBS)

gen-lppl:	gen-lppl.c $(OBJ) $(LMAR)
	$(CC) $(INC) $(CFLAGS) -o $@ $(OBJ) gen-lppl.c $(LIBDIR) $(LIBS) 

gen-noise:	gen-noise.c $(OBJ)
	$(CC) $(INC) $(CFLAGS) -o $@ $(OBJ) gen-noise.c $(LIBDIR) $(LIBS) 

utest:	utest.c $(OBJ) $(LMAR)
	$(CC) $(INC) $(CFLAGS) -o $@ $(OBJ) utest.c $(LIBDIR) $(LIBS) -lcunit

.PHONY: clean
clean:
	rm -f $(OBJ) lppl.o explore-lppl.o lppl explore-lppl guess-lppl gen-lppl utest

.PHONY: fromscratch
fromscratch: clean all

%.o:    %.c Makefile lppl.h 
	$(CC) $(INC) $(CFLAGS) -c -o $@ $<
