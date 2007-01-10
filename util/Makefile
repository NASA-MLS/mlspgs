# $RCSfile$

# Compiler vendor name, to select machine_$(CVN)
CVN = lf95

FC= lf95
FOPTS= -g --chk
LDOPTS=
PLATFORM= LF95.Linux.slow
#PGSLIB= /van1/vsnyder/mlspgs/SDPToolkit/LF95
CC=cc
#HDF5INC=/users/pwagner/apps/test/vrh9.0-lf95/hdf5-1.4.4/include
#HDF5LIB=/users/pwagner/apps/test/vrh9.0-lf95/hdf5-1.4.4/lib

FAFTER=2>&1|hl
LAFTER=2>&1|hl
MLSF95=LF95
MLSPLAT=Linux
MLSCONFG=LF95.Linux

PGSINC=/software/toolkit/E-LF6.1-tk5.2.11.A/toolkit5.2.11.A/include
PGSTK=/software/toolkit/E-LF6.1-tk5.2.11.A/toolkit5.2.11.A/lib/linux
HDFEOS=/software/toolkit/E-LF6.1-tk5.2.11.A/hdfeos/lib/linux
HDFEOS5=/software/toolkit/E-LF6.1-tk5.2.11.A/hdfeos5/lib/linux
HDF=/software/toolkit/E-LF6.1-tk5.2.11.A/HDF4.1r5/lib
HDF5=/software/toolkit/E-LF6.1-tk5.2.11.A/hdf5-1.6.2/lib
GCTP=/software/toolkit/E-LF6.1-tk5.2.11.A/hdfeos/lib/linux
PGSLIB=/software/toolkit/E-LF6.1-tk5.2.11.A/toolkit5.2.11.A/lib/linux
HDFEOS_LIB=/software/toolkit/E-LF6.1-tk5.2.11.A/hdfeos/lib/linux
HDFLIB=/software/toolkit/E-LF6.1-tk5.2.11.A/HDF4.1r5/lib

machine_$(CVN).o machine.mod: machine_$(CVN).f90
	$(FC) $(FOPTS) -c machine_$(CVN).f90

compare: compare.f90 machine.mod machine_$(CVN).o
	$(FC) $(FOPTS) -o compare compare.f90 machine_$(CVN).o

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

CopyElevMap: CopyElevMap.f90
	$(FC) $(FOPTS) -o CopyElevMap CopyElevMap.f90

Dump_HDF5_Structure.o: Dump_HDF5_Structure.f90 MyMLSDataInfo.mod
	$(FC) -c $(FOPTS) -I . -I ../lib/$(MLSCONFG) -I $(HDF5)/../fortran/src \
          -I $(HDF5)/../include Dump_HDF5_Structure.f90 $(FAFTER)

Dump_HDF5_Structure: Dump_HDF5_Structure.o MyMLSDataInfo.o
	$(FC) $(LDOPTS) -o Dump_HDF5_Structure \
         Dump_HDF5_Structure.o MyMLSDataInfo.o \
         -L ../lib/$(MLSCONFG) -lmls \
         -L $(PGSLIB) -lPGSTK \
         -L $(HDF5) -lhdf5_fortran -lhdf5  \
         -lz \
         -ljpeg \
         -lz -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.3 -lgcc $(LAFTER)

DumpPFAData: DumpPFAData.f90
	$(FC) $(FOPTS) -o DumpPFAData -I $(HOME)/mlspgs/lib/$(PLATFORM) \
          DumpPFAData.f90 -L $(HOME)/mlspgs/lib/$(PLATFORM) -lmls

Filter_Fit: Filter_Fit.f90
	$(FC) $(FOPTS) -o Filter_Fit Filter_Fit.f90 \
          -L/van1/vsnyder/math77/fortran/LF95/opt/ -lmath77 \
          -L/van1/vsnyder/math77/fortran/blas/LF95/opt/ -lblas

f90tex: f90tex.f90 machine.mod machine_$(CVN).o
	$(FC) $(FOPTS) -o f90tex f90tex.f90 machine_$(CVN).o -static

init_gen: init_gen.f90
	$(FC) $(FOPTS) -o init_gen init_gen.f90 \
        -I $(HOME)/mlspgs/lib/$(PLATFORM) \
          -L $(HOME)/mlspgs/lib/$(PLATFORM) -lmls

h5cat: h5cat.c
	$(CC) -g -o h5cat -I $(HDF5INC) h5cat.c \
          -L $(HDF5LIB) -lz -lm -lhdf5

MyMLSDataInfo.mod: MyMLSDataInfo.o

MyMLSDataInfo.o: MyMLSDataInfo.f90
	$(FC) -c $(FOPTS) -I . -I ../lib/$(MLSCONFG) -I $(HDF5)/../fortran/src \
          -I $(HDF5)/../include MyMLSDataInfo.f90 $(FAFTER)

# $Log$
# Revision 1.7  2004/09/23 23:16:11  vsnyder
# Simplify make for compare
#
# Revision 1.6  2004/09/22 19:22:50  vsnyder
# Repair make for f90tex
#
# Revision 1.5  2004/06/17 01:05:26  vsnyder
# Account for ConvertPFA, DumpPFAData
#
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
