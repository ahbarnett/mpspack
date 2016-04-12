# GNU Makefile for MPSpack tweaks only.
# Compiles fortran/MEX codes, FMM interface, etc.
# Note that MPSpack can work without them.
#
# Tested on linux only.
# Bare-bones for git 4/12/16
# (C) 2008 - 2016 Alex Barnett

include make.inc

.PHONY: all tar clean

default: all

all:
	(cd @utils; make)

tar:
	(cd ..; tar zcvf mpspack.tgz mpspack)

clean:
	(cd @utils; make clean)
