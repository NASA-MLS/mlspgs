#include "/usr/share/pvm3/include/pvm3.h"

/* Compile me with:
cc -o killmaster killmaster.c -L/usr/share/pvm3/lib/LINUX -lpvm3
*/

main ( int argc, void *argv[] )
{
  int tid;
  int bufID;
  int info;

  tid = atoi ( argv[1] );
  bufID = pvm_initsend ( PvmDataDefault );
  info = pvm_pkbyte ( &tid, 1, 1);
  info = pvm_send ( tid, 999 );
}
