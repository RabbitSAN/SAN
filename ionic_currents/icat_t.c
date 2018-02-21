// calc T-type iCa

#include "sa.h"

////////// tau and ss 
void
tss_dt( REAL ee, REAL *tau, REAL *ss )
{
  REAL adt, bdt;

    adt = exp((ee + 26.3f) /30) *1068;
    bdt = exp(-(ee + 26.3f) /30) *1068;
    
    *tau = 1./(adt+bdt);
    *ss  = 1./(exp(-(ee+37)/6.8) +1);
}///tss_dt

////////// tau and ss 
void
tss_ft( REAL ee, REAL *tau, REAL *ss )
{
  REAL aft, bft;

    aft = exp(-(ee + 71.7)/83.3) *15.3;
    bft = exp((ee + 71.7) /15.38) *15;

    *tau = 1./(aft + bft);
    *ss  = 1./(exp((ee + 71.)/9.) +1.);

}////tss_ft

///////////////////////////////////////////////////////////////////////
REAL
icat_t(struct State *S, struct State *Sn, REAL ht, struct Table T, struct Cpar *C )
{
  static int first=1; 
  static REAL Ecat; 
  static REAL tdt[NTABLE],dtss[NTABLE], tft[NTABLE],ftss[NTABLE]; 
  REAL tdt_t, dtss_t, tft_t, ftss_t; 


  // first time
  for( ; first; first=0 )
    {
      int i;

      Ecat = 45;

      //fill tables
      for( i=-1; ++i<NTABLE; )
	{
	  REAL ee;
	  
	  ee = ETMIN+ (ETMAX-ETMIN)*i/(NTABLE-1);

	  // fill tabs
	  tss_dt( ee, tdt+i, dtss+i );
	  tss_ft( ee, tft+i, ftss+i );
	}//fori
    }//forfirst


  // estimate table values
  {
    register int iii=T.n1;
    tdt_t  = tdt[iii] + T.frac*(tdt[iii+1]-tdt[iii]);
    dtss_t = dtss[iii] + T.frac*(dtss[iii+1]-dtss[iii]);
    tft_t  = tft[iii] + T.frac*(tft[iii+1]-tft[iii]);
    ftss_t = ftss[iii] + T.frac*(ftss[iii+1]-ftss[iii]);
  }

  Sn->dt = dtss_t - (dtss_t - S->dt)*exp(-ht/tdt_t);	// solve eq
  Sn->ft = ftss_t - (ftss_t - S->ft)*exp(-ht/tft_t);	// solve eq

  return  C->gcat*S->ft*S->dt*(S->E-Ecat);

} /** icat **/
