# GNU Makefile for fortran/MEX compilation for MPSpack MATLAB toolbox
#
# Currently linux OS (glnx86) only.
#
# Barnett 9/5/08, updated creation of archive & getrevisionnumber 1/21/09
# simplified 8/5/09

# (C) 2008 - 2012 Alex Barnett and Timo Betcke

#rev := $(shell ./getrevisionnumber)
#pkg := mpspack-r$(rev)
pkg := mpspack-1.31

default: all

all:
	(cd @utils; make)
special:
	(cd @utils; make special)

.PHONY: all special clean

tar:
#	echo $(rev); echo $(pkg)
#       note sequential here to avoid race conditions...
	(cd ..; svn export mpspack $(pkg); tar zcvf $(pkg).tar.gz $(pkg))
	(cd ..; zip -r $(pkg).zip $(pkg))
	rm -Rf ../$(pkg)

clean:
	(cd @utils; make clean)
