# GNU Makefile for MEX compilation for MPSPack MATLAB toolbox
#
# Currently linux OS (glnx86) only.
#
# Barnett 9/5/08

default: all

all:
	(cd @utils; make)
special:
	(cd @utils; make special)

.PHONY: all special clean

tar:
	(cd ..; tar zcvf mpspack.tar.gz mpspack)

clean:
	(cd @utils; make clean)
