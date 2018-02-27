// calc L-type iCa

#include "sa.h"

////////// tau and ss 
void
tssACH_dl( REAL ee, REAL *tau, REAL *ss )
{
  REAL adl, bdl;

#define ETOL 1e-4 	// tolerance for E

  if( fabs(ee+35) <ETOL ) 
    adl = 14.19*2.5-42.45*ee /(exp(-0.208*ee )-1);
  else if( fabs(ee) <ETOL )
    adl=-14.19*(ee +35)/(exp(-(ee +35)/2.5)-1)+42.45*0.208;
  else	// normal
    adl = -14.19*(ee+35)/(exp(-(ee+35)/2.5)-1.) - 42.45*ee/(exp(-0.208*ee)-1.); 
  // check for zeros
  if( fabs(ee-5) <ETOL )
    bdl = 5.71/0.4;
  else		//normal
    bdl = 5.71*(ee-5.)/(exp(0.4*(ee-5))-1); 
  
  *tau = 1./(adl+bdl);
  *ss  = 1./(1+exp(-(ee+23.1)/6));

}////tssACH_dl


////////// tau and ss 
void
tssACH_fl( REAL ee, REAL *tau, REAL *ss )
{
  REAL afl, bfl;


  //check for zeros
  if( fabs(ee +28) <ETOL )
    afl = 3.12*4;
  else		// normal
    afl = 3.12*(ee+28)/(exp((ee+28)*0.25) -1);

  bfl = 25./(1+exp(-(ee+28)*0.25)); //
	  
  *tau = 1./(afl+bfl);
  *ss  = 1./(1+exp((ee+45)*0.2));

}////tssACH_fl





////////////////////////////////////////////////////////////////////////
REAL
icalach_t(struct State *S, struct State *Sn, REAL ht, struct Table T,  struct Cpar *C, REAL ach )
{
  static int first=1; 
  static REAL Ecal; 
  static REAL tdl[NTABLE],dlss[NTABLE], tfl[NTABLE],flss[NTABLE]; 
  REAL tdl_t, dlss_t, tfl_t, flss_t; 
  // gca should be here, but it proved to be varied from center out

  static REAL gach; 	// factor [0..1] to show deprression Ical by ACH

  // first time
  for(; first; first=0 )
    {
      int i;

      Ecal = 46.4; //mV

      //fill tables
      for( i=-1; ++i<NTABLE; )
	{
	  REAL ee;
	  
	  ee = ETMIN+ (ETMAX-ETMIN)*i/(NTABLE-1);

	  // fill tabs
	  tssACH_dl( ee, tdl+i, dlss+i );
	  tssACH_fl( ee, tfl+i, flss+i );

	}//fori      
    }//forfirst

      // ACH depression

	 gach = 1. - C->kachical* ach/(1.2e-7 +ach);

  // estimate table values
  {
    register int iii=T.n1;
    tdl_t  = tdl[iii] + T.frac*(tdl[iii+1]-tdl[iii]);
    dlss_t = dlss[iii] + T.frac*(dlss[iii+1]-dlss[iii]);
    tfl_t  = tfl[iii] + T.frac*(tfl[iii+1]-tfl[iii]);
    flss_t = flss[iii] + T.frac*(flss[iii+1]-flss[iii]);
  }

  Sn->dl = dlss_t - (dlss_t - S->dl)*exp(-ht/tdl_t);	// solve eq
  Sn->fl = flss_t - (flss_t - S->fl)*exp(-ht/tfl_t);	// solve eq

  // check for reality
//  if( tfl_t<ht || tdl_t<ht ) printf( "!!!ERROR: ht, tdl, tfl = %g %g %g\n", ht, tdl_t, tfl_t );


  // mintau
  mintau = MIN(mintau, tdl_t );
  mintau = MIN(mintau, tfl_t );


  return gach*C->gcal*( S->fl*S->dl + 0.006/(1+exp(-(S->E+14.1)/6)) )*(S->E-Ecal);
} /** icalach_t **/
