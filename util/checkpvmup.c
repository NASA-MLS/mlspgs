/*
 * Copyright 2005, by the California Institute of Technology. ALL
 * RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
 * commercial use must be negotiated with the Office of Technology Transfer
 * at the California Institute of Technology.

 * This software may be subject to U.S. export control laws. By accepting this
 * software, the user agrees to comply with all applicable U.S. export laws and
 * regulations. User has the responsibility to obtain export licenses, or other
 * export authority as may be required before exporting such information to
 * foreign countries or providing access to foreign persons.
 */
#include <stdlib.h>
#include <pvm3.h>

/* Compile me with:
cc -o checkpvmup checkpvmup.c -L$PVM_ROOT/lib/$PVM_ARCH -lpvm3
*/

main ( int argc, void *argv[] )
{
  int mytid;
  mytid = pvm_mytid ();
  exit ( mytid < 0 );
}
