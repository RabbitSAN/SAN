// calc ina

#include "sa.h"

////////// tau and ss 
void
tss_m( REAL ee, REAL *tau, REAL *ss )
{

  *ss  = pow(1+ exp(-(ee+30.32)/5.46), -0.333333333 );  // f-code
  *tau = 6.247e-4/ (0.832*exp(-0.335*(ee+56.7))+0.627*exp(0.082*(ee+65.01)))+4e-5;

}///tss_m

////////// tau and ss 
void
tss_h1( REAL ee, REAL *tau, REAL *ss )
{
  *ss  = 1./(1+exp((ee+66.1)/6.4));
  *tau = 3.717e-6*exp(-0.2815*(ee+17.1)) / (1+3.732e-3*exp(-0.3426*(ee+37.76))) + 5.977e-4;

}///tss_h1

////////// tau and ss 
void
tss_h2( REAL ee, REAL *tau, REAL *ss )
{
  *ss  = 1./(1+exp((ee+66.1)/6.4));
  *tau = 3.186e-8*exp(-0.6219*(ee+18.8)) / (1+7.189e-5*exp(-0.6683*(ee+34.07))) + 3.556e-3;

}///tss_h2

///////////////////////////////////////////////////////////////////////
REAL
ina_t(struct State *S, struct State *Sn, REAL ht, struct Table T, struct Cpar *C, struct Nernst *nernst)
{
  static int first=1; 
  static REAL Ena; 
  static REAL tm[NTABLE],mss[NTABLE];
  static REAL th1[NTABLE],h1ss[NTABLE];
  static REAL th2[NTABLE],h2ss[NTABLE];
  REAL mss_t, tm_t, fna, h, h1ss_t, th1_t, h2ss_t, th2_t;

  // first time
  for( ; first; first=0 )
    {
      //fill tables
      { int i; REAL ee;
      for( i=-1; ++i<NTABLE; )
	{
	  ee = ETMIN+ (ETMAX-ETMIN)*i/(NTABLE-1);
	  // tabs
	  tss_m( ee, tm+i, mss+i );
	  tss_h1( ee, th1+i, h1ss+i );
	  tss_h2( ee, th2+i, h2ss+i );
	}//fori
      }//fill
    }//forfirst

   // estimate table values
  {
    register int iii=T.n1;
    tm_t  = tm[iii] + T.frac*(tm[iii+1]-tm[iii]);
    mss_t = mss[iii] + T.frac*(mss[iii+1]-mss[iii]);
    th1_t  = th1[iii] + T.frac*(th1[iii+1]-th1[iii]);
    h1ss_t = h1ss[iii] + T.frac*(h1ss[iii+1]-h1ss[iii]);
    th2_t  = th2[iii] + T.frac*(th2[iii+1]-th2[iii]);
    h2ss_t = h2ss[iii] + T.frac*(h2ss[iii+1]-h2ss[iii]);
  }
  

  Sn->m = mss_t - (mss_t-S->m)*exp(-ht/tm_t);	// solve eq
  Sn->h1 = h1ss_t - (h1ss_t-S->h1)*exp(-ht/th1_t);	// solve eq
  Sn->h2 = h2ss_t - (h2ss_t-S->h2)*exp(-ht/th2_t);	// solve eq

  fna = 9.52e-2 *exp(-6.3e-2*(S->E+34.4)) / (1.+1.66*exp(-0.225*(S->E+63.7))) + 8.69e-2;

  h = (1-fna)*S->h1 + fna*S->h2;


  return C->gna*S->m*S->m*S->m*h*nao*FRD/RTF* (exp((S->E-nernst->ENa)/RTF)-1)/(exp(S->E/RTF)-1) *S->E;

} /** ina **/
