// Check processor time spent by code between calls; and print it 
// flag: 's' - start; '0' - reset; 'f' - finish; 'p' - print->stderr
//
// Rubin, September 2003

#include <time.h>
#include <stdio.h>


unsigned long
clocks( char flag )
{
  static unsigned long nn, totalt, thiscall;
  switch( flag )
    {
    case 's': 					// start
      thiscall = clock();
      break;
    case 'f':					// finish
      totalt += clock() - thiscall;
      nn++;
      break;
    case '0':
      nn = totalt = 0L;				// reset
      break;
    case 'p':					//print
      fprintf( stderr, "#clocks: %g /ms/; n=%ld\n", 1000.*totalt/CLOCKS_PER_SEC, nn );
      break;
    }//switch 

  return (unsigned long)totalt;
}
