// calc if, Hyperpolarization-activated current


#include "sa.h"

////////// tau and ss 
void
tss_y( REAL ee, REAL *tau, REAL *ss )
{
  REAL ay, by;

  ay = exp(-(ee + 78.91)/26.62);
  by = exp((ee + 75.13)/21.25);

  *tau = 1./ (ay + by);
  *ss = *tau * ay;
}///tss_y

///////////////////////////////////////////////////////////////////////
REAL
if_t(struct State *S, struct State *Sn, REAL ht, struct Table T, struct Cpar *C, struct Is *I, struct Nernst *nernst)
{
  static int first=1; 
  static REAL ty[NTABLE],yss[NTABLE]; 
  REAL ty_t, yss_t, ifna, ifk;

  // first time
  for( ; first; first=0 )
    {

      //fill tables
      { int i; REAL ee;
      for( i=-1; ++i<NTABLE; )
	{
	  ee = ETMIN+ (ETMAX-ETMIN)*i/(NTABLE-1);
	  // tabs
	  tss_y( ee, ty+i, yss+i );
	}//fori
      }//filltables
    }//forfirst
  ///////////////////////////////////////////////////


  // estimate table values
  {
    register int iii=T.n1;
    ty_t  = ty[iii] + T.frac*(ty[iii+1]-ty[iii]);
    yss_t = yss[iii] + T.frac*(yss[iii+1]-yss[iii]);
  }
 
    
    Sn->y = yss_t - (yss_t - S->y)*exp(-ht/ty_t);	// solve eq

    /// fifna, fifk -- extern vars to keep track of these currents separately
    I->fifna = ifna = C->gfna*S->y*(S->E- nernst->ENa);
    I->fifk  = ifk  = C->gfk*S->y*(S->E- nernst->EK);

    return ifna+ifk;

} /** if_t **/
