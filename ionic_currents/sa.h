// basic variables and descriptions
// !!!!!!!!!!!! test for UNIVERSAL cell
#ifndef SA_H
#define SA_H 1

#include <math.h>

#define MIN(_x,_y)  ( ((_x)<(_y)) ? (_x) : (_y) )
#define SQR(_x)     ( (_x)*(_x) )

typedef double REAL;		// our type

REAL mintau;

struct State { REAL E, m, h1, h2, fl, dl, ft, dt, q, r, paf, pas, pii, xs, y, jach, kach,  qa, qi, nai, ki; };
struct Caintra_state {REAL ftc, ftmc, ftmm, fcmi, fcms, fcq, cai, casub, caup, carel; };
struct Cpar {REAL ctype,  cm, gna, gto, gsus, gkr, gks, gfna, gfk, gbna, gbca, gbk, ipss, knaca, icapmax, kachical, kachicat, kachikach, gcal, gcat, vc, vrel, vup, vsub, vi; };
struct Nernst {REAL ENa, EK, ECa;};
struct Is {REAL fibgna, fibgca, fibgk, fical, ficat, fif, fikr, fiks, fina, finaca, fip, fisus, fito, fikach, fist, ficap, fifna, fifk;};

#define NTABLE 1000
#define ETMIN (-100.)
#define ETMAX 100.
struct Table { int n1; float frac; } T;
struct Tlocal { REAL to; REAL uo; REAL tn; REAL un; };

// intracell concentrations; const for this model
#define nao	140.		// mM
#define cao	2.		// mM
#define ko		5.4		// mM	//normal K

#define cmC 	20e-6		// muF
#define gnaC	0.		// muS
#define gtoC	4.91e-3		// muS
#define gsusC	6.65e-5		// muS
#define gkrC	7.97e-4		// muS
#define gksC	5.18e-4		// muS
#define gfnaC	0.0548e-2		// muS
#define gfkC	0.0548e-2		// muS
#define gbnaC	5.8e-5		// muS
#define gbcaC	1.32e-5		// muS  //SAV
#define gbkC	2.52e-5		// muS
#define ipssC	0.0520064		// nA //adjusted!
#define knacaC	0.25e-5		// nA     //default
#define gcalC    0.7604e-2*1
#define gcatC    0.214e-2 		// gcal+gcat=const
#define icapmaxC 0.0132		   /*nA*/ //test

#define  kachicalC	0.14
#define  kachicatC	0.14
#define  kachikachC  	0.00705

#define cmP 	50e-6			// muF   
#define gnaP	9.23077e-07//*0.05		// muS
#define gtoP	0.0202722			// muS
#define gsusP	0.00541935		// muS
#define gkrP	0.0104523			// muS
#define gksP	0.00509804		// muS
#define gfnaP	0.00530769		// muS //orig
#define gfkP	0.00530769		// muS //orig
#define gbnaP	0.000145385		// muS
#define gbcaP	3.30769e-05		// muS
#define gbkP	6.3e-05			// muS
#define ipssP	0.88				// nA  //INCREASED   
#define knacaP	6.76923e-05		// nA  //INCREASED   

#define gcalP    0.0506923*0.5  		// muS
#define gcatP    0.0106923   		// muS
#define icapmaxP 0.066			   /*nA*/ //test


#define  kachicalP	0.56
#define  kachicatP	0.56
#define  kachikachP  	0.0792

/* Consts for Solution of Conductance */
#define R 	8314.      	// Universal Gas Constant (J/kmol*K)
#define FRD 	96485. 		// Faraday's Constant (C/mol)
#define TEMP 	310.   		// Temperature (K)
#define RTF (R*TEMP/FRD)	// const

extern REAL ht;      					// Time step (ms)
extern REAL t;       					// Time (ms)

extern REAL fifna, fifk;
//extern REAL cm,gna,gto,gsus,gkr,gks,gfna,gfk,gbna,gbca,gbk,ipss,knaca, kachical, kachicat, kachikach, gcal, gcat, icapmax, vc,vrel,vup,vsub,vi;

#endif //ifndef SA_H
