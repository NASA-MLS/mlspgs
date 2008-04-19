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

vpath %.f $(MATH77)/

%.o: %.f
	$(FC) -c -o $@ $(FOPTS) $<

all: Mie_Tables

# amach has a test that is subverted by optimization
amach.o: $(MATH77)/amach.f
	$(FC) -c -g $(MATH77)/amach.f

ice_size_distribution.mod Ice_Size_Distribution.o: Ice_Size_Distribution.f90
	$(FC) -c $(FOPTS) Ice_Size_Distribution.f90 -I $(HOME)/mlspgs/lib/$(PLAT)

mie_efficiencies_m.mod Mie_Efficiencies_m.o: Mie_Efficiencies_m.f90 sphbess.mod
	$(FC) -c $(FOPTS) Mie_Efficiencies_m.f90 -I $(HOME)/mlspgs/lib/$(PLAT)

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

clean:
	rm -f Mie_Efficiencies_m.o Mie_Tables.o \
              amach.o dcopy.o \
              dint1.o dinta.o dintdl.o dintdu.o dintf.o dintm.o \
              dintma.o dintns.o dinto.o dintop.o dintsm.o dmess.o mess.o \
              Ice_Size_Distribution.o ice_size_distribution.mod \
              mie_efficiencies_m.mod Phase_m.o phase_m.mod P1_m.o p1_m.mod \
              sphbes.o sphbess.mod

dcopy.o: $(MATH77)/dcopy.f
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
mess.o: $(MATH77)/mess.f

# $Log$
# Revision 1.2  2008/04/19 01:37:27  vsnyder
# Removed test harness makes
#
# Revision 1.1  2008/04/19 01:15:26  vsnyder
# Initial commit
#
