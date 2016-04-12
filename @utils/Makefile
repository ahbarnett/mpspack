# GNU Makefile for MEX compilation for MPSpack MATLAB toolbox's @utils directory
#
# These are only needed for tweaks to accelerate math libraries; MPSpack
# will run without them. Tested on linux only.
#
# Simplified version 4/12/16.
# (C) 2008-2016 Alex Barnett

include ../make.inc

default: all

all: gslbesselj gslbesseljnu greengardrokhlinhank103 greengardrokhlinhank106 inpolyc

.PHONY: all clean

gslbesselj: gslbesselj.c
	mex gslbesselj.c $(GSL) $(BLAS)

gslbesseljnu: gslbesseljnu.c
	mex gslbesseljnu.c $(GSL) $(BLAS)

greengardrokhlinhank103: greengardrokhlinhank103.c hank103.o
	mex greengardrokhlinhank103.c hank103.o LD=$(FC)

greengardrokhlinhank106: greengardrokhlinhank106.c hank103.o hank106.o
	mex greengardrokhlinhank106.c hank103.o hank106.o LD=$(FC)

hank103.o: hank103.f
	$(FC) $(FFLAGS) -c hank103.f

hank106.o: hank106.f
	$(FC) $(FFLAGS) -c hank106.f

inpolyc: inpolyc.c
	mex inpolyc.c

# remember not to delete the shipped Luong MEX executables (!):
clean:
	rm -f *.o greengardrokhlinhank103.mex* greengardrokhlinhank106.mex* gslbesseljnu.mex* gslbesselj.mex* inpolyc.mex*
