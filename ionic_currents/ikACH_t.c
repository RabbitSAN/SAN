// ACh dependent iK

#include "sa.h"

////////// tau and ss 
void
tss_jach( REAL ee, REAL *tau, REAL *ss )
{
  REAL a, b;

  a = 73.1; // /1/s/
  b = 120./(1.+exp(-(ee + 50.)/15.));
    
  *tau 	= 1./(a+b);
  *ss 	= a* *tau;
}///tss_jach

////////// tau and ss 
void
tss_kach( REAL ee, REAL *tau, REAL *ss )
{
  REAL a, b;

  a = 3.7; // /1/s/
  b = 5.82/(1.+exp(-(ee + 50.)/15.));
    
  *tau 	= 1./(a+b);
  *ss 	= a* *tau;
}///tss_kach



///////////////////////////////////////////////////////////////////////
REAL
ikach_t(struct State *S, struct State *Sn, REAL ht, struct Table T, struct Cpar *C, struct Nernst *nernst, REAL ach)
{
  static int first=1; 
  static REAL tj[NTABLE],jss[NTABLE], tk[NTABLE],kss[NTABLE]; 
  REAL tj_t, jss_t, tk_t, kss_t; 
  static REAL coefach;

  // first time
  for( ; first; first=0 )
    {
      int i;


      //fill tables
      for( i=-1; ++i<NTABLE; )
	{
	  REAL ee;
	  
	  ee = ETMIN+ (ETMAX-ETMIN)*i/(NTABLE-1);

	  // fill tabs
	  tss_jach( ee, tj+i, jss+i );
	  tss_kach( ee, tk+i, kss+i );
	}//fori

#define gkachmax (3*0.0198) // muS
#define nkach	1.5
#define k05kach	2.8e-7 // M
    }//forfirst

  coefach = C->kachikach *(ach/(2e-7 +ach))*ko/(10+ko);

  // estimate table values
  {
    register int iii=T.n1;
    tj_t  = tj[iii] + T.frac*(tj[iii+1]-tj[iii]);
    jss_t = jss[iii] + T.frac*(jss[iii+1]-jss[iii]);
    tk_t  = tk[iii] + T.frac*(tk[iii+1]-tk[iii]);
    kss_t = kss[iii] + T.frac*(kss[iii+1]-kss[iii]);
  }

  Sn->jach = jss_t - (jss_t - S->jach)*exp(-ht/tj_t);	// solve eq
  Sn->kach = kss_t - (kss_t - S->kach)*exp(-ht/tk_t);	// solve eq

  return coefach*S->jach*S->kach*(S->E-nernst->EK)/(1+exp( (S->E -nernst->EK-140)/(2.5*RTF) ) );

} /** ikach_t **/
