#include <stdio.h>
#include <regex.h>

/* Highlight a pattern in the input by coloring it in the output. */
/* "regex" is used to find the patterns.                          */

  static char Id[] = "$Id$";

main ( int argc, char* argv[] )
{ char b[256];                     /* input buffer */
  char after[] = "\033[0m";        /* black on transparent background */
  char blue[] = "\033[00;34;1m";   /* bold blue on transparent background */
  char puce[] = "\033[00;35;1m";   /* bold puce on transparent background */
  char red[] = "\033[00;31;1m";    /* bold red on transparent background */
  int i;                           /* Subscript/loop inductor */

#define PSIZ 256         /* Compiled pattern buffer size */

  typedef struct
  { char* find;          /* Pattern to find */
    char* color;         /* String to but before it */
    regex_t preg[PSIZ];  /* Compiled pattern -- from regcomp */
  } pat;

  regmatch_t match;      /* Where is the string that matches the pattern? */

  pat pats[] = { "[0-9]+-[SU]", red, {},     /* lf95 errors */
                 "[0-9]+-W", blue, {},       /* lf95 warnings */
                 "[0-9]+-I", puce, {},       /* lf95 informative */
                 /* NAG f95 patterns: */
                 "Deleted", red, {},       "Error", red, {},
                 "Fatal", red, {},         "Panic", red, {},
                 "Warning", blue, {},      "Extension", blue, {},
                 "Obsolescent", blue, {},  "Info", puce, {},
                 "undefined", red, {} };     /* during linking */

#define NPAT sizeof pats / sizeof pats[0]

  /* Compile the patterns */
  for ( i=0; i<NPAT ; i++ ) regcomp ( pats[i].preg, pats[i].find, REG_EXTENDED );

  setvbuf ( stdout, (char*)NULL, _IONBF, 1 ); /* unbuffer the output */

  while (1)
  { /* Get a line: */
    if ( fgets ( b, sizeof b, stdin ) == NULL ) return(0);
    /* Look for a pattern match */
    for ( i=0; i<NPAT; i++ )
    { if ( regexec ( pats[i].preg, b, 1, &match, 0 ) == 0 )
      /* Got a match; put the desired color around it */
      { fwrite ( b, sizeof(char), match.rm_so, stdout );
        printf ( "%s", pats[i].color );
        fwrite ( &b[match.rm_so], sizeof(char), match.rm_eo - match.rm_so, stdout );
        printf ( "%s%s", after, &b[match.rm_eo] );
        goto cycle;
      }
    }
    /* No match; just echo the input */
    printf ( "%s", b );
  cycle:
  }
}

/*
$Log$
*/
