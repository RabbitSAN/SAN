// Sustained component of 4-AP sensitive current.
#include "sa.h"


////////// tau and ss 
void
tss_r( REAL ee, REAL *tau, REAL *ss )
{
  *ss  = 1./(exp(-(ee-10.93)/19.7) +1);
  *tau = 15.59e-3/( 1.037*exp((ee+30.61)*0.09) + 0.369*exp(-0.12*(ee+23.84)) ) +2.98e-3;
}///tss_r

///////////////////////////////////////////////////////////////////////
REAL
isus_t(struct State *S, struct State *Sn, REAL ht, struct Table T, struct Cpar *C, struct Nernst *nernst)
{
  static int first=1; 
  static REAL tr[NTABLE],rss[NTABLE];
  REAL tr_t, rss_t;

  // first time
  for( ; first; first=0 )
    {
      //fill tables
      { int i; REAL ee;
      for( i=-1; ++i<NTABLE; )
	{
	  ee = ETMIN+ (ETMAX-ETMIN)*i/(NTABLE-1);
	  // tabs
	  tss_r( ee, tr+i, rss+i );
	}//fori
      }//fill
    }//forfirst

   // estimate table values
  {
    register int iii=T.n1;
    tr_t  = tr[iii] + T.frac*(tr[iii+1]-tr[iii]);
    rss_t = rss[iii] + T.frac*(rss[iii+1]-rss[iii]);
  }
 
    Sn->r = rss_t - (rss_t - S->r)*exp(-ht/tr_t);	// solve eq


  return C->gsus*S->r*(S->E-nernst->EK);

} /** isus_t **/
