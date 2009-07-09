# GNU Makefile for fortran/MEX compilation for MPSpack MATLAB toolbox
#
# Currently linux OS (glnx86) only.
#
# Barnett 9/5/08, updated creation of archive & getrevisionnumber 1/21/09

rev := $(shell ./getrevisionnumber)
pkg := mpspack-r$(rev)

default: all

all:
	(cd @utils; make)
special:
	(cd @utils; make special)

.PHONY: all special clean

tar:
#	echo $(rev); echo $(pkg)
	(cd ..; svn export mpspack $(pkg); tar zcvf $(pkg).tar.gz $(pkg); rm -Rf $(pkg))

clean:
	(cd @utils; make clean)
