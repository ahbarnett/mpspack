# Make.inc

CC = gcc
FC = gfortran
GSL = -lgsl # Adapt to the location of your GSL library

# Set the location of the blas libraries on your system
# If you have ATLAS installed uncomment one of the following and
# change path if necessary

# Multithreaded version
# BLAS = -L/usr/lib64/atlas -lptf77blas -lptcblas -latlas
# Singlethreaded version
# BLAS = -L/usr/lib64/atlas -lf77blas -lcblas -latlas

# Generic blas package (if installed)
BLAS = -lblas
