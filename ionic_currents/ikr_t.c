// calc i K r, rapid delayed rectifying

#include "sa.h"

////////// tau and ss 
void
tss_paf( REAL ee, REAL *tau, REAL *ss )
{
    *ss = 1./(exp(-(ee+14.2)/10.6) + 1);
    *tau = 1./( exp((ee-9)/15.9) *37.2 + exp(-(ee-9)/22.5) *0.96 );

}///tss_paf

////////// tau and ss 
void
tss_pas( REAL ee, REAL *tau, REAL *ss )
{
    *ss  = 1./(exp(-(ee+14.2)/10.6) + 1);
    *tau = 1./( exp((ee-9)/17) *4.2 + exp(-(ee-9)/21.6) *0.15);
}////tss_pass

////////// tau and ss 
void
tss_pii( REAL ee, REAL *tau, REAL *ss )
{
    *ss  = 1./(exp((ee+18.6)/10.1) + 1);
    *tau = 0.002;
}////tss_pii



///////////////////////////////////////////////////////////////////////
REAL
ikr_t(struct State *S, struct State *Sn, REAL ht, struct Table T, struct Cpar *C, struct Nernst *nernst )
{
  static int first=1; 
  static REAL tpaf[NTABLE],pafss[NTABLE], tpas[NTABLE],passs[NTABLE], tpii[NTABLE],piiss[NTABLE]; 
  REAL pafss_t, tpaf_t, passs_t, tpas_t, piiss_t, tpii_t, fkr, pa;

  // first time
  for( ; first; first=0 )
    {
      int i;


      //fill tables
      for( i=-1; ++i<NTABLE; )
	{
	  REAL ee;
	  ee = ETMIN+ (ETMAX-ETMIN)*i/(NTABLE-1);
	  // tabs
	  tss_paf( ee, tpaf+i, pafss+i );
	  tss_pas( ee, tpas+i, passs+i );
	  tss_pii( ee, tpii+i, piiss+i );
	}//fori
    }//forfirst

  // estimate table values
  {
    register int iii=T.n1;
    tpaf_t  = tpaf[iii] + T.frac*(tpaf[iii+1]-tpaf[iii]);
    pafss_t = pafss[iii] + T.frac*(pafss[iii+1]-pafss[iii]);
    tpas_t  = tpas[iii] + T.frac*(tpas[iii+1]-tpas[iii]);
    passs_t = passs[iii] + T.frac*(passs[iii+1]-passs[iii]);
    tpii_t  = tpii[iii] + T.frac*(tpii[iii+1]-tpii[iii]);
    piiss_t = piiss[iii] + T.frac*(piiss[iii+1]-piiss[iii]);
  }


    Sn->paf = pafss_t - (pafss_t - S->paf)*exp(-ht/tpaf_t);	// solve eq
    Sn->pas = passs_t - (passs_t - S->pas)*exp(-ht/tpas_t);	// solve eq
    Sn->pii = piiss_t - (piiss_t - S->pii)*exp(-ht/tpii_t);	// solve eq

    fkr = 0.4;
    pa = (1-fkr)*S->paf + fkr*S->pas;

 
    return C->gkr*pa*S->pii*(S->E-nernst->EK);

} /** ikr_t **/
