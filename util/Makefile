# $RCSfile$

# Compiler vendor name, to select machine_$(CVN)
CVN = lf95

FC= lf95
FOPTS= -g --chk --info
PLATFORM= LF95.Linux.slow
PGSLIB= /van1/vsnyder/mlspgs/SDPToolkit/LF95
CC=cc
HDF5INC=/users/pwagner/apps/test/vrh9.0-lf95/hdf5-1.4.4/include
HDF5LIB=/users/pwagner/apps/test/vrh9.0-lf95/hdf5-1.4.4/lib

compare.o: compare.f90 machine.mod
	$(FC) $(FOPTS) -c compare.f90

machine_$(CVN).o machine.mod: machine_$(CVN).f90
	$(FC) $(FOPTS) -c machine_$(CVN).f90

compare: compare.o machine_$(CVN).o
	$(FC) $(FOPTS) -o compare compare.o machine_$(CVN).o

CondenseLeakLog: CondenseLeakLog.f90
	$(FC) $(FOPTS) -o CondenseLeakLog -I $(HOME)/mlspgs/lib/$(PLATFORM) \
          CondenseLeakLog.f90 -L $(HOME)/mlspgs/lib/$(PLATFORM) -lmls

ConvertPFA: ConvertPFA.f90
	$(FC) $(FOPTS) -o ConvertPFA -I $(HOME)/mlspgs/lib/$(PLATFORM) \
          ConvertPFA.f90 -L $(HOME)/mlspgs/lib/$(PLATFORM) -lmls

convert_spectroscopy: convert_spectroscopy.f90
	$(FC) $(FOPTS) -o convert_spectroscopy -I $(HOME)/mlspgs/lib/$(PLATFORM) \
          convert_spectroscopy.f90 -L $(HOME)/mlspgs/lib/$(PLATFORM) -lmls \
          -L $(PGSLIB) -lPGSTK

DumpPFAData: DumpPFAData.f90
	$(FC) $(FOPTS) -o DumpPFAData -I $(HOME)/mlspgs/lib/$(PLATFORM) \
          DumpPFAData.f90 -L $(HOME)/mlspgs/lib/$(PLATFORM) -lmls

Filter_Fit: Filter_Fit.f90
	$(FC) $(FOPTS) -o Filter_Fit Filter_Fit.f90 \
          -L/van1/vsnyder/math77/fortran/LF95/opt/ -lmath77 \
          -L/van1/vsnyder/math77/fortran/blas/LF95/opt/ -lblas

f90tex: f90tex.f90
	$(FC) $(FOPTS) -o f90tex f90tex.f90 \
          -I $(HOME)/mlspgs/lib/$(PLATFORM) \
          -L $(HOME)/mlspgs/lib/$(PLATFORM) -lmls -static

init_gen: init_gen.f90
	$(FC) $(FOPTS) -o init_gen init_gen.f90 \
        -I $(HOME)/mlspgs/lib/$(PLATFORM) \
          -L $(HOME)/mlspgs/lib/$(PLATFORM) -lmls

h5cat: h5cat.c
	$(CC) -g -o h5cat -I $(HDF5INC) h5cat.c \
          -L $(HDF5LIB) -lz -lm -lhdf5

# $Log$
# Revision 1.4  2003/08/13 18:16:45  livesey
# Typo
#
# Revision 1.3  2003/08/13 18:16:36  livesey
# New HDF5 library.
#
# Revision 1.2  2001/07/25 02:07:24  vsnyder
# Account for init_gen
#
# Revision 1.1  2001/06/01 21:52:45  vsnyder
# Initial commit
#
