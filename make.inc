# Make.inc: common compiler and library locations for MPSpack makefiles
# Barnett, probably 2008.

CC = gcc
FC = gfortran
FFLAGS = -fPIC -O3          # PIC needed for MEX to link against .o

# LIBRARY LOCATIONS...
GSL = -lgsl              # Adapt to the location of your GSL library

# Set the location of the BLAS libraries on your system
# If you have ATLAS installed uncomment one of the following and
# change path if necessary
BLAS = -lblas                  # Generic blas package (if installed)
# Multithreaded version
# BLAS = -L/usr/lib64/atlas -lptf77blas -lptcblas -latlas
# Singlethreaded version
# BLAS = -L/usr/lib64/atlas -lf77blas -lcblas -latlas

# Fast algorithm libs:
FMM2D = /usr/local/fmmlib2d/     # point to your FMM installation
# (get from http://www.cims.nyu.edu/cmcl/fmm2dlib/fmm2dlib.html )
LP2D = /home/alex/physics/leslie/gimbutas/lp2d   # your LP2D installation
# (get from https://math.dartmouth.edu/~ahb/software/lp2d.tgz )
