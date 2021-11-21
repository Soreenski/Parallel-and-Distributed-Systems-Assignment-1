CC=gcc
MPICC=mpicc
CILKCC=/usr/local/OpenCilk-9.0.1-Linux/bin/clang
CFLAGS=-O3
PTHREADSFLAGS = -O3 -pthread -std=c99


default: all

sequential: mmio.o coo2csc.o sequential.o 
	$(CC) $(CFLAGS) -o sequential mmio.c coo2csc.c sequential.c

pthreads: mmio.o coo2csc.o pthreads.c
	$(CC) $(PTHREADSFLAGS) -o pthreads mmio.c coo2csc.c pthreads.c

openck: mmio.o coo2csc.o openck.c
	$(CILKCC) $(CFLAGS) -o openck mmio.c coo2csc.c openck.c -fcilkplus

openmp: mmio.o coo2csc.o openmp.c
	$(CC) $(CFLAGS) -o openmp mmio.c coo2csc.c openmp.c -fopenmp

test: mmio.o coo2csc.o test.c
	$(CC) $(PTHREADSFLAGS) -o test mmio.c coo2csc.c test.c

.PHONY: clean

all: sequential pthreads openck openmp

clean: 
	rm -f sequential pthreads openck openmp sequential.o pthreads.o openck.o openmp.o example mmio.o coo2csc.o