/* icap as in Demir etc.
 */
//Rubin, December 2004
// icapmax -> variable
//Rubin, October 2007


#include "sa.h"


REAL
icap(struct Cpar *C, struct Caintra_state *Ca )
{
  return C->icapmax/(1+(0.0004/Ca->casub));
}/** icap **/
