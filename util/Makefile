# $RCSfile$

FC= lf95
FOPTS= -g --chk --info
PLATFORM= LF95.Linux
PGSLIB= /van1/vsnyder/mlspgs/SDPToolkit/LF95

CondenseLeakLog: CondenseLeakLog.f90
	$(FC) $(FOPTS) -o CondenseLeakLog -I $(HOME)/mlspgs/lib/$(PLATFORM) \
          CondenseLeakLog.f90 -L $(HOME)/mlspgs/lib/$(PLATFORM) -lmls

convert_spectroscopy: convert_spectroscopy.f90
	$(FC) $(FOPTS) -o convert_spectroscopy -I $(HOME)/mlspgs/lib/$(PLATFORM) \
          convert_spectroscopy.f90 -L $(HOME)/mlspgs/lib/$(PLATFORM) -lmls \
          -L $(PGSLIB) -lPGSTK

Filter_Fit: Filter_Fit.f90
	$(FC) $(FOPTS) -o Filter_Fit Filter_Fit.f90 \
          -L/van1/vsnyder/math77/fortran/LF95/opt/ -lmath77 \
          -L/van1/vsnyder/math77/fortran/blas/LF95/opt/ -lblas

f90tex: f90tex.f90
	$(FC) $(FOPTS) -o f90tex f90tex.f90 \
          -I $(HOME)/mlspgs/lib/$(PLATFORM) \
          -L $(HOME)/mlspgs/lib/$(PLATFORM) -lmls

init_gen: init_gen.f90
	$(FC) $(FOPTS) -o init_gen init_gen.f90 \
        -I $(HOME)/mlspgs/lib/$(PLATFORM) \
          -L $(HOME)/mlspgs/lib/$(PLATFORM) -lmls

# $Log$
# Revision 1.1  2001/06/01 21:52:45  vsnyder
# Initial commit
#
