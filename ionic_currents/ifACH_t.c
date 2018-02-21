// calc if, Hyperpolarization-activated current
// Shift activation curve due to ACH

// use nonconst nai, ki

#include "sa.h"
#include <math.h>

////////// tau and ss 
void
tssACH_y( REAL ee, REAL *tau, REAL *ss, REAL ach )
{
  REAL ay, by;
  REAL achshft;

  // ACH consts for If
#define smax 	(-7.5) 	// /mv/  (Zhang's .f)
#define nf	0.69	// exponent
#define k05	1.26e-8	// /M/ ACh
  
  achshft = smax* pow(ach,nf) /( pow(k05,nf) + pow(ach,nf) );


  ay = exp(-(ee + 78.91 -achshft)/26.62);
  by = exp( (ee + 75.13 -achshft)/21.25);

  *tau = 1./ (ay + by);
  *ss = *tau * ay;
}///tss_y

///////////////////////////////////////////////////////////////////////
REAL
ifach_t(struct State *S, struct State *Sn, REAL ht, struct Table T, struct Cpar *C , struct Is *I, struct Nernst *nernst, REAL ach)
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
	  tssACH_y( ee, ty+i, yss+i, ach );
	}//fori
      }//filltables
    }//forfirst
  ///////////////////////////////////////////////////

  //The table above is not used now to account for non constant ACh, tssACH_y is called at each timestep instead.
    tssACH_y( S->E, &ty_t, &yss_t, ach );
    
    Sn->y = yss_t - (yss_t - S->y)*exp(-ht/ty_t);	// solve eq

    /// fifna, fifk -- extern vars to keep track of these currents separately
    I->fifna = ifna = C->gfna*S->y*(S->E - nernst->ENa);
    I->fifk  = ifk  = C->gfk*S->y*(S->E - nernst->EK);

    return ifna+ifk;

} /** ifach_t **/
