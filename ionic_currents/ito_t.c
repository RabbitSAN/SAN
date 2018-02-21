// calc Ito - transient
#include "sa.h"



////////// tau and ss 
void
tss_q( REAL ee, REAL *tau, REAL *ss )
{
  *ss = 1./ (exp((ee + 59.37)/13.1) + 1);
  *tau = 10.1e-3 + 65.17e-3/( 0.5686*exp(-0.08161*(ee+49.0))+0.7174*exp(0.2719*(ee+50.93)) );

}///tss_q


///////////////////////////////////////////////////////////////////////
REAL
ito_t(struct State *S, struct State *Sn, REAL ht, struct Table T, struct Cpar *C, struct Nernst *nernst)
{
  static int first=1; 
  static REAL tq[NTABLE],qss[NTABLE];
  REAL tq_t, qss_t;


  // first time
  for( ; first; first=0 )
    {

      //fill tables
      { int i; REAL ee;
      for( i=-1; ++i<NTABLE; )
	{
	  ee = ETMIN+ (ETMAX-ETMIN)*i/(NTABLE-1);
	  // tabs	
	  tss_q( ee, tq+i, qss+i );
	}//fori
      }//fill
    }//forfirst

  // estimate table values
  {
    register int iii=T.n1;
    tq_t  = tq[iii] + T.frac*(tq[iii+1]-tq[iii]);
    qss_t = qss[iii] + T.frac*(qss[iii+1]-qss[iii]);
  }

  Sn->q = qss_t - (qss_t - S->q)*exp(-ht/tq_t);	// solve eq



  // r vriable is already changed in isus
  //   Sn->r = rss_t - (rss_t - S->r)*exp(-ht/tr_t);	// solve eq

  return C->gto*S->q*S->r*(S->E-nernst->EK);

} /** ito_t **/
