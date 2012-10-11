# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# "$Id$"

.PHONY: Mie_Tables Mie_Tables_NoHDF clean tar AllDVI

Mie_Tables:
	rm -f Write*HDF_m.o writehdf_m.mod
	make -f Makefile.second $@

Mie_Tables_NoHDF:
	rm -f Write*HDF_m.o writehdf_m.mod
	make -f Makefile.second $@ NO=No

clean tar AllDVI:
	make -f Makefile.second $@

# $Log$
# Revision 1.5  2010/05/07 00:52:32  vsnyder
# New top-level makefile
#
