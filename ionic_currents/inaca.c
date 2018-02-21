// calc Inaca

#include "sa.h"

REAL
inaca(struct State *S, struct State *Sn, REAL ht, struct Cpar *C, struct Caintra_state *Ca )
{
  static int first=1; 
  static REAL  gnaca, dnaca; 
  static REAL nai3, nao3;

  // first time
  for( ; first; first=0 )
    {
      // consts
      gnaca = 0.5;
      dnaca = 1e-4;

    }//forfirst


  nai3 = S->nai*S->nai*S->nai;
  nao3 = nao*nao*nao;

  return C->knaca*( nai3*cao  * exp(S->E*0.03743*gnaca) - nao3*Ca->casub * exp(S->E*0.03743f*(gnaca-1)))  / (1+ dnaca*(nao3*Ca->casub+nai3*cao));


} /** inaca **/
