# GNU Makefile for MEX compilation for MPSPack MATLAB toolbox's @utils directory
#
# Currently linux OS only. Needs GNU Scientific Library (gsl) installed
#
#
# Copyright (C) 2008, 2009, Alex Barnett, Timo Betcke

include ../make.inc

default: all

all: gslbesselj gslbesseljnu greengardrokhlinhank103 inpolyc
special: greengardrokhlinhank106

.PHONY: all special clean

gslbesselj: gslbesselj.c
	mex gslbesselj.c $(GSL) $(BLAS)

gslbesseljnu: gslbesseljnu.c
	mex gslbesseljnu.c $(GSL) $(BLAS)

greengardrokhlinhank103: greengardrokhlinhank103.c hank103.o
	mex greengardrokhlinhank103.c hank103.o LD=$(FC)

greengardrokhlinhank106: greengardrokhlinhank106.c hank103.o hank106.o
	mex greengardrokhlinhank106.c hank103.o hank106.o LD=$(FC)

hank103.o: hank103.f
	$(FC) -fPIC -c hank103.f -O3

hank106.o: hank106.f
	$(FC) -fPIC -c hank106.f -O3

inpolyc: inpolyc.c
	mex inpolyc.c

clean:
	rm -f *.o *.mexglx *.mexa64
