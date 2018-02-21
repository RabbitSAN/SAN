/* Intracellular Ca dynamics as in: 
   [Dynamical description of sinoatrial node pacemaking: 
    improved mathematical model for primary pacemaker cell.
    Yasutaka Kurata, Ichiro Hisatome, Sunao Imanishi, and Toshishige Shibamoto
    Am J Physiol Heart Circ Physiol 283: H2074­H2101, 2002;]
*/
// Rubin, September 2003
#include <stdio.h>

#include "sa.h"
#include "ca.h"



REAL ca_intra(REAL ht, REAL icatotal, struct Cpar *C, struct Caintra_state *Ca )
{
  REAL dftc, dftmc, dftmm, dfcmi, dfcms, dfcq;	// time derivatives
  REAL jcadiff, jrel, jup, jtr; 		// fluxes


  ht *= 1000; 	// time in ms (was s)


  // Ca++ buffering -- derivatives on t
  dftc 	= kftc*Ca->cai*(1.-Ca->ftc)-kbtc*Ca->ftc;
  dfcmi = kfcm*Ca->cai*(1.-Ca->fcmi)-kbcm*Ca->fcmi;
  dfcms = kfcm*Ca->casub*(1.-Ca->fcms)-kbcm*Ca->fcms;
  dfcq 	= kfcq*Ca->carel*(1.-Ca->fcq)-kbcq*Ca->fcq;
  dftmc = kftmc*Ca->cai*(1.-Ca->ftmc-Ca->ftmm)-kbtmc*Ca->ftmc;
  dftmm = kftmm*mgi*(1.-Ca->ftmc-Ca->ftmm)-kbtmm*Ca->ftmm; 


  // Ca++ fluxes: diffusion and because of SR
  jcadiff = (Ca->casub-Ca->cai)/tdiffca;
  jrel 	= prel*(Ca->carel-Ca->casub)/(1.+(krel*krel/(Ca->casub*Ca->casub)));
  jup	= pup/(1.+kup/Ca->cai);
  jtr	= (Ca->caup-Ca->carel)/ttr;

  // Now, Ca concentrations

  Ca->cai = Ca->cai + ht*( (jcadiff*C->vsub-jup*C->vup)/C->vi - (cmtot*dfcmi + tctot*dftc + tmctot*dftmc) );

  Ca->casub	= Ca->casub + ht*( (   -1000*(icatotal)/(2.*FRD)  + jrel*C->vrel)/C->vsub -jcadiff - cmtot*dfcms );

  Ca->carel 	= Ca->carel + ht*( jtr - jrel  -cqtot*dfcq  );

  Ca->caup 	= Ca->caup + ht*( jup - jtr*C->vrel/C->vup );


  // fractional occupancies
  Ca->ftc  	= Ca->ftc  + ht* dftc;
  Ca->fcmi  	= Ca->fcmi + ht* dfcmi;
  Ca->fcms  	= Ca->fcms + ht* dfcms;
  Ca->fcq  	= Ca->fcq  + ht* dfcq;
  Ca->ftmc  	= Ca->ftmc + ht* dftmc;
  Ca->ftmm  	= Ca->ftmm + ht* dftmm;


  //  return casub_new;
  return Ca->cai;
} /** ca_intra **/


