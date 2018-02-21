// Main module to cal SA (central)

// The full set of equations can be found in [Aliev RR. Conceptual and Detailed Models of Electrical Activity of the Myocardium. 
// Lambert Academic Publishing: Saarbrucken, Germany. 2012. ISBN:978-3-8465-3943-9 (in Russian).]


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "sa.h"

#define READIC 0	//If set to 1 read initial state from state.dat file. If set to 0 use default initial conditions.
#define ACH 1		//Turn on/off ACh.
#define WRITEIC 0	//Write new IC to the state.dat file after simulations.
#define IONS 1		//Turn on/off intracellular ionic concentrations changes.

// extern subs
extern REAL ibgna(struct State *, struct State *, REAL, struct Cpar *, struct Nernst *);
extern REAL ibgca(struct State *, struct State *, REAL, struct Cpar *, struct Caintra_state *, struct Nernst *);
extern REAL ibgk(struct State *, struct State *, REAL, struct Cpar *, struct Nernst *);
extern REAL ical_t(struct State *, struct State *, REAL, struct Table, struct Cpar *);
extern REAL icalach_t(struct State *, struct State *, REAL, struct Table, struct Cpar *, REAL ach);
extern REAL icat_t(struct State *, struct State *, REAL, struct Table, struct Cpar *);
extern REAL ikr_t(struct State *, struct State *, REAL, struct Table, struct Cpar *, struct Nernst *);
extern REAL iks_t(struct State *, struct State *, REAL, struct Table, struct Cpar *);
extern REAL ikach_t(struct State *, struct State *, REAL, struct Table, struct Cpar *, struct Nernst *, REAL ach);
extern REAL if_t(struct State *, struct State *, REAL, struct Table, struct Cpar *, struct Is *, struct Nernst *);
extern REAL ifach_t(struct State *, struct State *, REAL, struct Table, struct Cpar *, struct Is *, struct Nernst *, REAL ach);
extern REAL ina_t(struct State *, struct State *, REAL, struct Table, struct Cpar *, struct Nernst *);
extern REAL isus_t(struct State *, struct State *, REAL, struct Table, struct Cpar *, struct Nernst *);
extern REAL ito_t(struct State *, struct State *, REAL, struct Table, struct Cpar *, struct Nernst *);
extern REAL inaca(struct State *, struct State *, REAL, struct Cpar *, struct Caintra_state *);
extern REAL ip(struct State *, struct State *, REAL, struct Cpar *);

extern REAL icap(struct Cpar *, struct Caintra_state *);
extern REAL ca_intra( REAL, REAL, struct Cpar *, struct Caintra_state *);
extern unsigned long clocks( char );

REAL basic_pars( REAL vthr, REAL vo1, REAL vn1, REAL t, REAL ht );
void settype(struct Cpar *);

REAL t, ht=0.01e-3;		// [s]

REAL itotal=0.;			// total current = sum all currents




REAL istim;			// [na] stimul current -- variable
 
REAL trun=50;//30;		// [s] run time 

//////////////////////////////////////////////////////////////////////////
// Set parameters corresponding to particular cell type.
// C=0 corresponds to central SAN.
// C=1 corresponds to peripheral SAN.
//////////////////////////////////////////////////////////////////////////
void settype(struct Cpar *C){
  C->cm=cmC+C->ctype*(cmP-cmC);
  C->gna=gnaC+C->ctype*(gnaP-gnaC);
  C->gto=gtoC+C->ctype*(gtoP-gtoC);
  C->gsus=gsusC+C->ctype*(gsusP-gsusC);
  C->gkr=gkrC+C->ctype*(gkrP-gkrC);
  C->gks=gksC+C->ctype*(gksP-gksC);
  C->gfna=gfnaC+C->ctype*(gfnaP-gfnaC);
  C->gfk=gfkC+C->ctype*(gfkP-gfkC);
  C->gbna=gbnaC+C->ctype*(gbnaP-gbnaC);
  C->gbca=gbcaC+C->ctype*(gbcaP-gbcaC);
  C->gbk=gbkC+C->ctype*(gbkP-gbkC);
  C->ipss=(ipssC+C->ctype*(ipssP-ipssC));
  C->knaca=knacaC+C->ctype*(knacaP-knacaC);
  C->icapmax=icapmaxC+C->ctype*(icapmaxP-icapmaxC);

  C->kachical=kachicalC+C->ctype*(kachicalP-kachicalC);
  C->kachicat=kachicatC+C->ctype*(kachicatP-kachicatC);
  C->kachikach=kachikachC+C->ctype*(kachikachP-kachikachC);

  C->vc=(0.11*C->cm*1e6);
  C->vrel=(0.0012*C->vc);
  C->vup=(0.0116*C->vc);
  C->vsub=(0.01*C->vc);
  C->vi=(0.46*C->vc-C->vsub);
	
  C->gcal=gcalC+C->ctype*(gcalP-gcalC);
  C->gcat=gcatC+C->ctype*(gcatP-gcatC);
  
}
//////////////////////////////////////////////////////////////////////////
// Total Current
//////////////////////////////////////////////////////////////////////////
REAL
fitotal(struct State *St, struct State *Stn,REAL ht, struct Cpar *Cp, struct Is *I, struct Table T, struct Caintra_state *Ca, struct Nernst *nernst, REAL ach)
{
  I->fibgna = ibgna(St, Stn, ht, Cp, nernst);
  I->fibgca = ibgca(St, Stn, ht, Cp,Ca, nernst);
  I->fibgk = ibgk(St, Stn, ht, Cp, nernst);
  if(1/* t>0 */) clocks('s');


#if ACH
  I->fical = icalach_t(St, Stn, ht, T, Cp, ach);
  I->ficat = icat_t(St, Stn, ht, T, Cp);
  I->fif = ifach_t(St, Stn, ht, T, Cp, I, nernst, ach);
  I->fikach = ikach_t(St, Stn, ht, T, Cp, nernst, ach);
#else //noACH
  I->fical = ical_t(St, Stn, ht, T, Cp);
  I->ficat = icat_t(St, Stn, ht, T, Cp);
  I->fif = if_t(St, Stn, ht, T, Cp, I, nernst);
#endif //ACH
  I->fikr = ikr_t(St, Stn, ht, T, Cp, nernst);
  I->fiks = iks_t(St, Stn, ht, T, Cp);
  I->fisus = isus_t(St, Stn, ht, T, Cp, nernst); //Keep isus and ito in this order for proper gating variables calculation.
  I->fito = ito_t(St, Stn, ht, T, Cp, nernst);
  I->fina = ina_t(St, Stn, ht, T, Cp, nernst);


  I->finaca = inaca(St, Stn, ht, Cp, Ca);
  I->fip = ip(St, Stn, ht, Cp);
  I->ficap = icap(Cp,Ca);
  return I->fibgna+I->fibgca+I->fibgk+I->fical+I->ficat+I->fif+I->fikr+I->fiks+I->fina+I->finaca+I->fip+I->fisus+I->fito +I->fikach +I->ficap;
}
///////////////////////////////////////////////////////////////////////////
// Func( t, E ) is a right hand side for eq dE/dt = Func
///////////////////////////////////////////////////////////////////////////
REAL Func( REAL t, struct State *Stn, struct State *St, struct Cpar *Cp, struct Is *I, struct Table T,  struct Caintra_state *Ca, struct Nernst *nernst, REAL ach){
  return -(istim+ fitotal(St, Stn, t, Cp, I, T,Ca, nernst, ach))/Cp->cm;  	// RH part +istim
}/** Func **/

///////////////////////////////////////////////////////////////////////////
// calculate Nernst potentials using intra- and extra- concentrations 
///////////////////////////////////////////////////////////////////////////
void	potentials(struct State *Stn, struct Caintra_state *Ca, struct Nernst *nernst)
{
  nernst->ENa 	= RTF* log(nao/Stn->nai);
  nernst->EK 		= RTF* log(ko/Stn->ki);
  nernst->ECa =	0.5*RTF*log(cao/Ca->cai); // zca=2;
} /** potentials **/

///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// Default initial conditions 
///////////////////////////////////////////////////////////////////////////

void initial_conditions(struct State *Stn,struct Caintra_state *Ca)
{
	Stn->E=-70.5; 
	Stn->m=0.086;
	Stn->h1=0.74; 
	Stn->h2=0.085;
	Stn->fl=0.995;
	Stn->dl=0.00034; 
	Stn->ft=0.55; 
	Stn->dt=0.0071; 
	Stn->q=0.74;
	Stn->r=0.015;
	Stn->paf=0.021;
	Stn->pas=0.28;
	Stn->pii=0.994;
	Stn->xs=0.097;
	Stn->y=0.077;
	Stn->jach=0.77;
	Stn->kach=0.59;
	Stn->qa=0.93;
	Stn->qi=0.011;
	Stn->nai=8.24;
	Stn->ki=140.;
	Ca->ftc=0.057; 
	Ca->ftmc=0.64;
	Ca->ftmm=0.32; 
	Ca->fcmi=0.11;
	Ca->fcms=0.033;
	Ca->fcq=0.21;
	Ca->cai=0.002; 
	Ca->casub=8.2e-5; 
	Ca->caup=1.66;
	Ca->carel=0.22;

}

int main(int argc, char **argv){

//Initialize structures, cell type, ach concentration.
    struct State St, Stn;
    struct Caintra_state Ca;
    struct Cpar Cp;
    struct Nernst nernst;
    struct Is I;

    Cp.ctype=atof(argv[1]);	//Read cell type from command line. 0 - central SAN, 1 - peripheral SAN.
    settype(&Cp);    
    REAL ach=0e-8; //2.5e-8;

#if READIC
  {
    FILE *fin = fopen( "state.dat", "r" );
    if(!fin) exit(puts("!!! Cannot open IC file"));
    fread( &Stn, sizeof(struct State), 1, fin );
    fread( &Ca, sizeof(struct Caintra_state), 1, fin );
    fclose(fin);
  }
#else
  {
	  initial_conditions(&Stn,&Ca);
  }
#endif
	FILE *fout_second = fopen( "out.txt", "w+" );
	fprintf(fout_second, "#E\tnai\tki\tcai\tcasub\tcaup\tcarel\tibgna\tibgca\tibgk\tical\ticat\tif\tikr\tiks\tina\tinaca\tip\tisus\tito\tikach\ticap\tifna\tifk\n");
  	int counter = 0;

  	// time-stepping 
  	for( t=0; t<trun; t+=ht )
  	{
      	potentials(&Stn, &Ca, &nernst);	// Nernst potentials


   	{
   		static float tgt=0.;
		if( t>=tgt )
		{
	       		tgt += 0.001;	// [s].  Change output frequency here. 
	       		fprintf(fout_second, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", t, Stn.E, Stn.nai, Stn.ki, Ca.cai, Ca.casub, Ca.caup, Ca.carel, I.fibgna, I.fibgca, I.fibgk, I.fical, I.ficat, I.fif, I.fikr, I.fiks, I.fina, I.finaca, I.fip, I.fisus, I.fito, I.fikach, I.ficap, I.fifna, I.fifk);

		}
   }//output
	    
	{St = Stn; mintau = 1e33;}

#if WRITEIC
	{static float twrite=0.;
		if( t>=twrite ){
			twrite+=100;
			FILE *fout = fopen( "state.dat", "w" );
			fwrite( &Stn, sizeof(struct State), 1, fout);
			fwrite( &Ca, sizeof(struct Caintra_state), 1, fout);
			fclose(fout);
		}
	}
#endif //WRITEIC


//T structure is used for lookup tables when currents are calculated.
      {
	int ii; double dd;
	dd = (St.E-ETMIN)/(ETMAX-ETMIN)*(NTABLE-1); 
	ii = (int)dd;
	
	T.n1 = ii;
	if( T.n1<0 || T.n1>=NTABLE ) fprintf(stderr, "!!!T.n1=%d\n", T.n1);
	T.frac = dd-ii;
      }
    

      // Forward Euler step
       Stn.E = St.E + ht*(Func(ht,&Stn, &St, &Cp, &I, T, &Ca, &nernst, ach));
       Ca.cai = ca_intra(ht, I.fical+I.ficat-2*I.finaca+I.ficap+I.fibgca, &Cp, &Ca);
#if IONS
       Stn.nai = St.nai+ht* -1e6*(I.fina+3.*I.finaca+3.*I.fip+I.fibgna+I.fifna)/(FRD*Cp.vi); 
       Stn.ki  = St.ki+ht* -1e6*(-2.*I.fip+I.fikr+I.fiks+I.fikach+I.fito+I.fisus+I.fifk+I.fibgk)/(FRD*Cp.vi);
#endif
       
       // tracking of the progress
       counter += 1;
       if ((counter%100000) == 0)
       {
       		printf("Step #%d\n", counter);
       }

    }//for_t
	
#if WRITEIC
	FILE *fout = fopen( "state.dat", "w" );
	fwrite( &Stn, sizeof(struct State), 1, fout);
	fwrite( &Ca, sizeof(struct Caintra_state), 1, fout);
	fclose(fout);
#endif //WRITEIC
  	
  	clocks('p');
  	fclose(fout_second);
	return 0;
}//main
