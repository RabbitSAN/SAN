/* consts for Intracellular Ca dynamics as in: 
   [Dynamical description of sinoatrial node pacemaking: 
    improved mathematical model for primary pacemaker cell.
    Yasutaka Kurata, Ichiro Hisatome, Sunao Imanishi, and Toshishige Shibamoto
    Am J Physiol Heart Circ Physiol 283: H2074­H2101, 2002;]
*/
// Rubin, December 2004

//struct Caintra_state {float ftc, ftmc, ftmm, fcmi, fcms, fcq, cai, casub, caup, carel; } Ca;   

#define	tdiffca	0.04	//
#define	ttr	60.0	//
//#define	ttr	30.0	//
#define	krel	0.0012	//
#define	pup	0.005	//
#define	prel	5.0	//
#define	kup	0.0006	//

#define	mgi	2.5	// Mg -- assume const

//#define	vi	1.5939	//--> in sa.h
//#define	vrel	0.0042	//
//#define	vup	0.0406	//
//#define	vsub	0.035	//
#define	tctot	0.031	//
#define	tmctot	0.062	//
#define	cmtot	0.045	//
#define	cqtot	10.0	//
#define	kftc	88.8	//
#define	kftmm	2.277	//
#define	kftmc	227.7	//
#define	kfcm	227.7	//
#define	kfcq	0.534	//
#define	kbtc	0.446	//
#define	kbtmc	0.00751 //
#define	kbtmm	0.751	//
#define	kbcm	0.542	//
#define	kbcq	0.445	//
