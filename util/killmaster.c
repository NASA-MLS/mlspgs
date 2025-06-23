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
cc -o killmaster killmaster.c -L/usr/share/pvm3/lib/LINUX -L/usr/share/pvm3/lib/LINUXI386 -lpvm3
invoke as 'killmaster 0x<tid>' where the master <tid> is obtained from the pvm "ps -a" command 
*/

main ( int argc, void *argv[] )
{
  int tid;
  int bufID;
  int info;

  tid = (int) strtod ( argv[1], NULL );
  printf ( "Task id is:%x\n", tid );
  bufID = pvm_initsend ( PvmDataDefault );
  info = pvm_pkint ( &tid, 1, 1);
  info = pvm_send ( tid, 999 );
}
