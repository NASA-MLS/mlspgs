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

math77 libmath77.a: amach.o dcopy.o dcpval.o derm1.o derv1.o \
        dint1.o dinta.o dintdl.o dintdu.o dintf.o \
        dintm.o dintma.o  dintns.o dinto.o dintop.o dintsm.o \
        dmess.o dpfit.o drotm.o drotmg.o \
        erfin.o ermsg.o ierm1.o ierv1.o mess.o
	ar r libmath77.a $^

%.o: %.f
	$(FC) -c -o $@ $(FOPTS) $<

# amach has a test that is subverted by optimization
amach.o: amach.f
	$(FC) -c -g $(AMACHOPTS) amach.f

clean:
	rm -f *.o *.a

# $Log$
