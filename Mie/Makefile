# Copyright 2005, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# "$Id$"

FC=     ifort
FOPTS= 
LDOPTS=
MATH77 = math77
PLAT=   IFC.Linux

.SUFFIXES:
.SUFFIXES: .o .dbg.o .f90 .f .c

vpath %.f math77/ $(MATH77)/

.f.o:
	$(FC) -c -o $@ $(FOPTS) $<

all: Mie_Efficiencies Mie_Least_Squares Mie_Tables \
     Test_Derivs Test_Ice_Distribution Test_Phase

# amach has a test that is subverted by optimization
amach.o: math77/amach.f
	$(FC) -c -g math77/amach.f

ice_size_distribution.mod Ice_Size_Distribution.o: Ice_Size_Distribution.f90
	$(FC) -c $(FOPTS) Ice_Size_Distribution.f90 -I $(HOME)/mlspgs/lib/$(PLAT)

mie_efficiencies_m.mod Mie_Efficiencies_m.o: Mie_Efficiencies_m.f90 sphbess.mod
	$(FC) -c $(FOPTS) Mie_Efficiencies_m.f90 -I $(HOME)/mlspgs/lib/$(PLAT)

Mie_Efficiencies_p.o: Mie_Efficiencies_p.f90 mie_efficiencies_m.mod
	$(FC) -c $(FOPTS) Mie_Efficiencies_p.f90 -I $(HOME)/mlspgs/fwdmdl/$(PLAT) \
              -I $(HOME)/mlspgs/lib/$(PLAT)

Mie_Least_Squares.o : Mie_Least_Squares.f90
	$(FC) -c $(FOPTS) Mie_Least_Squares.f90 -I $(HOME)/mlspgs/lib/$(PLAT)

Mie_Tables.o: Mie_Tables.f90 ice_size_distribution.mod \
              mie_efficiencies_m.mod phase_m.mod p1_m.mod
	$(FC) -c $(FOPTS) Mie_Tables.f90 -I $(HOME)/mlspgs/fwdmdl/$(PLAT) \
              -I $(HOME)/mlspgs/lib/$(PLAT)

p1_m.mod P1_m.o: P1_m.f90
	$(FC) -c $(FOPTS) P1_m.f90 -I $(HOME)/mlspgs/lib/$(PLAT)

phase_m.mod Phase_m.o: Phase_m.f90 p1_m.mod
	$(FC) -c $(FOPTS) Phase_m.f90 -I $(HOME)/mlspgs/lib/$(PLAT)

sphbes.o sphbess.mod: sphbes.f90
	$(FC) -c $(FOPTS) sphbes.f90 -I $(HOME)/mlspgs/lib/$(PLAT)

Mie_Efficiencies: Mie_Efficiencies_p.o Mie_Efficiencies_m.o sphbes.o
	$(FC) -o Mie_Efficiencies $(LDOPTS) \
              Mie_Efficiencies_p.o Mie_Efficiencies_m.o \
              sphbes.o \
              $(HOME)/mlspgs/lib/$(PLAT)/Constants.o \
              $(HOME)/mlspgs/lib/$(PLAT)/MLSKinds.o \
              $(HOME)/mlspgs/lib/$(PLAT)/physics.o \
              $(HOME)/mlspgs/fwdmdl/$(PLAT)/RefractiveIndex.o 

Mie_Least_Squares: Mie_Least_Squares.o \
            amach.o daxpy.o dcopy.o dcov2.o ddot.o dhfti.o dhtcc.o dhtgen.o \
            dnrm2.o dswap.o erfin.o ermor.o ermsg.o ierm1.o ierv1.o
	$(FC) -o Mie_Least_Squares  Mie_Least_Squares.o \
              amach.o daxpy.o dcopy.o dcov2.o ddot.o dhfti.o dhtcc.o dhtgen.o \
              dnrm2.o dswap.o erfin.o ermor.o ermsg.o ierm1.o ierv1.o \
              $(HOME)/mlspgs/lib/$(PLAT)/Constants.o \
              $(HOME)/mlspgs/lib/$(PLAT)/MLSKinds.o

Mie_Tables: Mie_Tables.o Mie_Efficiencies_m.o \
            amach.o dcopy.o \
            dint1.o dinta.o dintdl.o dintdu.o dintf.o \
            dintm.o dintma.o  dintns.o dinto.o dintop.o dintsm.o \
            dmess.o mess.o \
            Ice_Size_Distribution.o Phase_m.o P1_m.o sphbes.o
	$(FC) -o Mie_Tables $(LDOPTS) \
              Mie_Tables.o Mie_Efficiencies_m.o \
              amach.o dcopy.o \
              dint1.o dinta.o dintdl.o dintdu.o dintf.o \
              dintm.o dintma.o  dintns.o dinto.o dintop.o dintsm.o \
              dmess.o mess.o \
              Ice_Size_Distribution.o Phase_m.o P1_m.o sphbes.o \
              $(HOME)/mlspgs/lib/$(PLAT)/Constants.o \
              $(HOME)/mlspgs/lib/$(PLAT)/MLSKinds.o \
              $(HOME)/mlspgs/lib/$(PLAT)/physics.o \
              $(HOME)/mlspgs/fwdmdl/$(PLAT)/RefractiveIndex.o 

Test_Derivs.o: Test_Derivs.f90 ice_size_distribution.mod \
               mie_efficiencies_m.mod
	$(FC) -c $(FOPTS) Test_Derivs.f90 -I $(HOME)/mlspgs/lib/$(PLAT)

Test_Derivs: Test_Derivs.o Ice_Size_Distribution.o Mie_Efficiencies_m.o \
             sphbes.o
	$(FC) -o Test_Derivs $(LDOPTS) Test_Derivs.o Ice_Size_Distribution.o \
              Mie_Efficiencies_m.o sphbes.o \
              $(HOME)/mlspgs/lib/$(PLAT)/Constants.o \
              $(HOME)/mlspgs/lib/$(PLAT)/MLSKinds.o \
              $(HOME)/mlspgs/lib/$(PLAT)/physics.o \
              $(HOME)/mlspgs/fwdmdl/$(PLAT)/RefractiveIndex.o 

Test_Ice_Distribution.o: Test_Ice_Distribution.f90 ice_size_distribution.mod
	$(FC) -c $(FOPTS) Test_Ice_Distribution.f90 -I $(HOME)/mlspgs/lib/$(PLAT)

Test_Ice_Distribution: Test_Ice_Distribution.o Ice_Size_Distribution.o
	$(FC) -o Test_Ice_Distribution $(LDOPTS) Test_Ice_Distribution.o \
              Ice_Size_Distribution.o \
              $(HOME)/mlspgs/lib/$(PLAT)/Constants.o \
              $(HOME)/mlspgs/lib/$(PLAT)/MLSKinds.o \
              $(HOME)/mlspgs/lib/$(PLAT)/physics.o \
              $(HOME)/mlspgs/fwdmdl/$(PLAT)/RefractiveIndex.o 

Test_Phase.o: Test_Phase.f90 mie_efficiencies_m.mod phase_m.mod p1_m.mod
	$(FC) -c $(FOPTS) Test_Phase.f90 -I $(HOME)/mlspgs/lib/$(PLAT)

Test_Phase: Test_Phase.o
	$(FC) -o Test_Phase Test_Phase.o Mie_Efficiencies_m.o Phase_m.o \
              P1_m.o sphbes.o \
              amach.o dcopy.o dint1.o dinta.o dintdl.o dintdu.o dintf.o \
              dintm.o dintma.o dintns.o dinto.o dintop.o dintsm.o \
              dmess.o mess.o \
              $(HOME)/mlspgs/lib/$(PLAT)/Constants.o \
              $(HOME)/mlspgs/lib/$(PLAT)/MLSKinds.o \
              $(HOME)/mlspgs/lib/$(PLAT)/physics.o \
              $(HOME)/mlspgs/fwdmdl/$(PLAT)/RefractiveIndex.o 

clean:
	rm -f Mie_Efficiencies_m.o Mie_Efficiencies_m.o Mie_Tables.o \
              Test_Derivs.o Test_Derivs Test_Phase.o Test_Phase \
              Test_Ice_Distribution.o Test_Ice_Distribution \
              amach.o dint1.o dinta.o dintdl.o dintdu.o dintf.o dintm.o \
              dintma.o dintns.o dinto.o dintop.o dintsm.o dmess.o mess.o \
              Ice_Size_Distribution.o ice_size_distribution.mod \
              mie_efficiencies_m.mod Phase_m.o phase_m.mod P1_m.o p1_m.mod \
              sphbes.o sphbess.mod

daxpy.o: $(MATH77)/daxpy.f
dcopy.o: $(MATH77)/dcopy.f
dcov2.o: $(MATH77)/dcov2.f
ddot.o: $(MATH77)/ddot.f
dhfti.o: $(MATH77)/dhfti.f
dhtcc.o: $(MATH77)/dhtcc.f
dhtgen.o: $(MATH77)/dhtgen.f
dint1.o: $(MATH77)/dint1.f
dinta.o: $(MATH77)/dinta.f
dintdl.o: $(MATH77)/dintdl.f
dintdu.o: $(MATH77)/dintdu.f
dintf.o: $(MATH77)/dintf.f
dintm.o: $(MATH77)/dintm.f
dintma.o: $(MATH77)/dintma.f
dintns.o: $(MATH77)/dintns.f
dinto.o: $(MATH77)/dinto.f
dintop.o: $(MATH77)/dintop.f
dintsm.o: $(MATH77)/dintsm.f
dmess.o: $(MATH77)/dmess.f
dnrm2.o: $(MATH77)/dnrm2.f
dswap.o: $(MATH77)/dswap.f
erfin.o: $(MATH77)/erfin.o
ermor.o: $(MATH77)/ermor.o
ermsg.o: $(MATH77)/ermsg.o
ierm1.o: $(MATH77)/ierm1.o
ierv1.o: $(MATH77)/ierv1.o
mess.o: $(MATH77)/mess.f

# $Log$
