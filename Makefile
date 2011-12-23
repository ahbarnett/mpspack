# GNU Makefile for fortran/MEX compilation for MPSpack MATLAB toolbox
#
# Currently linux OS (glnx86) only.
#
# Barnett 9/5/08, updated creation of archive & getrevisionnumber 1/21/09
# simplified 8/5/09

# (C) 2008 - 2011 Alex Barnett and Timo Betcke

#rev := $(shell ./getrevisionnumber)
#pkg := mpspack-r$(rev)
pkg := mpspack-1.2beta

default: all

all:
	(cd @utils; make)
special:
	(cd @utils; make special)

.PHONY: all special clean

tar:
#	echo $(rev); echo $(pkg)
	(cd ..; /opt/CollabNet_Subversion/bin/svn export mpspack $(pkg); tar zcvf $(pkg).tar.gz $(pkg); rm -Rf $(pkg))

clean:
	(cd @utils; make clean)
