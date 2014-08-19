/* Copyright 2005, by the California Institute of Technology. ALL
   RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
   commercial use must be negotiated with the Office of Technology Transfer
   at the California Institute of Technology.

   This software may be subject to U.S. export control laws. By accepting this
   software, the user agrees to comply with all applicable U.S. export laws and
   regulations. User has the responsibility to obtain export licenses, or other
   export authority as may be required before exporting such information to
   foreign countries or providing access to foreign persons.
*/

#include <sys/time.h>
#include <sys/resource.h>

static char RCSfile[] = "$RCSfile$";

static char Id[] = "$Id$";

int GetResourceUsage ( struct rusage *usage )
{

/* A wrapper for getrusage with "who" == RUSAGE_SELF */

  return getrusage ( RUSAGE_SELF, usage ) ;

}

/* $Log$
 * Revision 2.1  2014/08/19 16:54:42  vsnyder
 * Initial commit
 * */
