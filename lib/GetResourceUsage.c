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
#include <unistd.h>

static char RCSfile[] = "$RCSfile$";

static char Id[] = "$Id$";

long GetPID ()
{

/* A wrapper for getpid, with a known result kind
   because we don't know the typedef for pid_t */

  long MyPID;
  MyPID = getpid();
  return MyPID;
}

int GetResourceUsage ( struct rusage *usage )
{

/* A wrapper for getrusage with "who" == RUSAGE_SELF
   because we don't know the value of RUSAGE_SELF */

  return getrusage ( RUSAGE_SELF, usage ) ;

}

long TicksPerSecond ()
{

/* A wrapper for sysconf(_SC_CLK_TCK)
   because we don't know the value of _SC_CLK_TCK */

  return sysconf ( _SC_CLK_TCK ) ;

}

/* $Log$
 * Revision 2.2  2014/08/19 23:13:29  vsnyder
 * Add GetPID and TicksPerSecond
 *
/* Revision 2.1  2014/08/19 16:54:42  vsnyder
/* Initial commit
/* */
