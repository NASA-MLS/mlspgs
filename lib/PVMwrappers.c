/* This is a set of cfortran.h wrappers for the pvm stuff.  It is
   needed to make the f90 bindings work right */

/* $Id$ */

#include "cfortran.h"

FCALLSCFUN1(INT, pvm_pkstr, PVM_PKSTR, pvm_pkstr, STRING);
     
FCALLSCFUN3(INT, pvm_pkdouble, PVM_PKDOUBLE, pvm_pkdouble, PDOUBLE,
	    INT, INT);

FCALLSCFUN3(INT, pvm_pkint, PVM_PKINT, pvm_pkint, PINT, INT, INT);

FCALLSCFUN1(INT, pvm_upkstr, PVM_UPKSTR, pvm_upkstr, PSTRING);
     
FCALLSCFUN3(INT, pvm_upkdouble, PVM_UPKDOUBLE, pvm_upkdouble, PDOUBLE,
	    INT, INT);

FCALLSCFUN3(INT, pvm_upkint, PVM_UPKINT, pvm_upkint, PINT, INT, INT);



