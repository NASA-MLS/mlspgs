/* $Id */

#include <stdio.h>
#include <stdlib.h>
#include "/usr/share/pvm3/include/pvm3.h"

/* Compile me with:
cc -o machineok machineok.c -L/usr/share/pvm3/lib/LINUX -lpvm3
*/

main ( int argc, void *argv[] )
{
  int tid;
  int bufID;
  int info;
  const int zeroseven[2] = { 0, 7 };

  if ( argc != 3 ) {
    printf ( "Usage: machineok <tid> <machine-name>\n" );
    return;
  }

  tid = (int) strtod ( argv[1], NULL );
  printf ( "Telling task 0x%x that %s is fixed.\n", tid, argv[2] );
  bufID = pvm_initsend ( PvmDataDefault );
  info = pvm_pkint ( &zeroseven, 2, 1);
  info = pvm_pkstr ( argv[2] );
  info = pvm_send ( tid, 800 );
}

/*
$Log$
/*
