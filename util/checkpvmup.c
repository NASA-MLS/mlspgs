#include <stdlib.h>
#include "/usr/share/pvm3/include/pvm3.h"

/* Compile me with:
cc -o checkpvmup checkpvmup.c -L/usr/share/pvm3/lib/$PVM_ARCH -lpvm3
*/

main ( int argc, void *argv[] )
{
  int mytid;
  mytid = pvm_mytid ();
  exit ( mytid < 0 );
}
