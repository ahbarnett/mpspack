# GNU Makefile for MEX compilation for MPSPack MATLAB toolbox's @utils directory
#
# Currently linux OS (glnx86) only. Needs GNU Scientific Library (gsl) installed
#
# Barnett 9/5/08

CC = gcc
FC = gfortran

default: all

all: gslbesselj.mexglx gslbesseljnu.mexglx greengardrokhlinhank103.mexglx
special: greengardrokhlinhank106.mexglx

.PHONY: all special clean

gslbesselj.mexglx: gslbesselj.c
	mex gslbesselj.c -lgsl -lblas

gslbesseljnu.mexglx: gslbesseljnu.c
	mex gslbesseljnu.c -lgsl -lblas

greengardrokhlinhank103.mexglx: greengardrokhlinhank103.c hank103.o
	mex greengardrokhlinhank103.c hank103.o LD=$(FC)

greengardrokhlinhank106.mexglx: greengardrokhlinhank106.c hank103.o hank106.o
	mex greengardrokhlinhank106.c hank103.o hank106.o LD=$(FC)

hank103.o: hank103.f
	$(FC) -c hank103.f -O3

hank106.o: hank106.f
	$(FC) -c hank106.f -O3

clean:
	rm *.o *.mexglx
