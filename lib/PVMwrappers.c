/* This is a set of cfortran.h wrappers for the pvm stuff.  It is
   needed to make the f90 bindings work right */

/* $Id$ */

#include "cfortran.h"
#include "string.h"

#define MAXPVMARGS 64
#define ARGLEN 64

     
/* This horrible collection of code is needed so I can get pvmspawn
   complete with command line arguments */

static int noArgs;
static char *args[MAXPVMARGS];

void clearpvmargs()
{
  noArgs = 0;
}

void nextpvmarg( char *str )
{
  int myLen;

  myLen = strlen(str);
  args[noArgs] = (char*) malloc(myLen+1);
  strncpy ( args[noArgs], str, myLen );
  noArgs++;
}

void freepvmargs()
{
  int i;
  for ( i = 0; i < noArgs; i++ ) free(args[i]);
  noArgs=0;
}

int mypvmspawn ( char *task, int flag, char *where, int ntask, int *tids )
{
  int i;

  args[noArgs]=NULL;
  return ( pvm_spawn ( task, args, flag, where, ntask, tids ) );
}

/* Now other wrappers */ 

FCALLSCFUN1(INT, pvm_pkstr, PVM_PKSTR, pvm_pkstr, STRING);
     
FCALLSCFUN3(INT, pvm_pkdouble, PVM_PKDOUBLE, pvm_pkdouble, PDOUBLE,
	    INT, INT);

FCALLSCFUN3(INT, pvm_pkint, PVM_PKINT, pvm_pkint, PINT, INT, INT);

FCALLSCFUN1(INT, pvm_upkstr, PVM_UPKSTR, pvm_upkstr, PSTRING);
     
FCALLSCFUN3(INT, pvm_upkdouble, PVM_UPKDOUBLE, pvm_upkdouble, PDOUBLE,
	    INT, INT);

FCALLSCFUN3(INT, pvm_upkint, PVM_UPKINT, pvm_upkint, PINT, INT, INT);

FCALLSCSUB1(usleep, USLEEP, usleep, INT)

FCALLSCSUB0(clearpvmargs, CLEARPVMARGS, clearpvmargs )

FCALLSCSUB0(freepvmargs, FREEPVMARGS, freepvmargs )

FCALLSCSUB1(nextpvmarg, NEXTPVMARG, nextpvmarg, PSTRING )

FCALLSCFUN5(INT, mypvmspawn, MYPVMSPAWN, mypvmspawn, 
	    STRING, INT, STRING, INT, PINT )
