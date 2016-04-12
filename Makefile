# GNU Makefile for MPSpack setup.
# Compiles fortran/MEX codes, FMM interface, etc.
# Note that MPSpack can work without them.
#
# Tested on linux only.
# Simplified to single makefile for git 4/12/16
# (C) 2008 - 2016 Alex Barnett

CC = gcc
FLIBS = -fPIC -O3          # PIC needed for MEX to link against .o
FC = gfortran
U = @utils                 # our utils directory

######## library locations: adjust as necessary (shown for ubuntu 14.04 linux)
BLAS = -lblas
GSL = -lgsl                      # GNU Scientific Library
FMM2D = /usr/local/fmmlib2d/     # point to your FMM installation
# (get from http://www.cims.nyu.edu/cmcl/fmm2dlib/fmm2dlib.html )
LP2D = /home/alex/physics/leslie/gimbutas/lp2d   # your LP2D installation
# (get from https://math.dartmouth.edu/~ahb/software/lp2d.tgz )

# decide if Rokhlin's hank106 available...
ifneq ("$(wildcard $(U)/hank106.f)","")
HANKELS = $(U)/hank103.o
else
HAVEHANK106 = 0
endif

# ***** finish this up!

all: hankel bessel inpolyc
.PHONY: all

hankel:

special:
	(cd @utils; make special)

clean:
	(cd @utils; make clean)
