FC = ifort
CC = cc

MKLDIR = /mn/stornext/u2/sigurdkn/local/intel_mkl/mkl/lib/intel64

# The include and linking commands for LAPACK.
LAPACK_LINK = -shared-intel -Wl,-rpath,$(MKLDIR) -L$(MKLDIR)  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread

HOMELIB = -L$(HOME)/.local/lib
HOMEINC = -I$(HOME)/.local/include
HKELIB = -L/mn/stornext/u2/hke/owl/local/lib
HKEINC = -I/mn/stornext/u2/hke/owl/local/include

FITSLIB = $(HOMELIB) -lcfitsio
HEALPIXLIB = $(HOMELIB) -lhealpix
QUIETLIB = $(HOMELIB) -lquiet
HDFLIB = $(HKELIB) -lhdf5_fortran -lhdf5
#MISCLIBS = -llapack -lblas -lhdf5_fortran -lhdf5
#MISCLIBS = $(HKELIB) -lblas

#LIBDIRS =  -L/mn/stornext/u2/hke/owl/local/lib
#LIBDIRS = -L$(HOME)/.local/lib \
#	  -L/mn/stornext/u2/sigurdkn/local/lib \
#	  -L/mn/stornext/u2/hke/local/lib
#	  -L$(HOME)/src/quiet/quiet_svn/oslo/src/f90/include \
#	  -L/mn/regulus/u1/sigurdkn/local/intel_mkl/mkl/lib/intel64/
#	  -L/mn/regulus/u1/sigurdkn/local/lib \

INCS = $(HOMEINC)
#INCS = 	-I/mn/stornext/u2/hke/owl/local/include
#INCS = -I$(HOME)/.local/include
#	-I/mn/stornext/u2/sigurdkn/local/include
#	-I/mn/stornext/u2/hke/local/include
#	-I$(HOME)/src/quiet/quiet_svn/oslo/src/f90/include \
#	-I/mn/regulus/u1/sigurdkn/local/include \
#

#LIBS = $(LAPACK_LINK) $(LIBDIRS) $(QUIETLIB) $(HEALPIXLIB) $(MISCLIBS) $(FITSLIB) -openmp -cxxlib
LIBS = $(HDFLIB) $(LAPACK_LINK) $(QUIETLIB) $(HEALPIXLIB) $(MISCLIBS) $(FITSLIB) -openmp -cxxlib

FFLAGS = -g -C -traceback -check all
#FFLAGS = -g -C
#FFLAGS = -O3
#FFLAGS = -O3

SIMMAPOBJ = simmap.o

simmap: $(SIMMAPOBJ)
	$(FC) -o $@ $^ $(LIBS)

%.o: %.f90
	$(FC) $(FFLAGS) $(INCS) -c $*.f90

.PHONY: clean
clean:
	rm -f *.mod *.o *~
