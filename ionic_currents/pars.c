/* Measure basic pars of stimul, print to stdout

Rubin, Jan 2004
*/

#include <stdio.h>
#include "sa.h"		// we only need REAL from there 


REAL
basic_pars( REAL vthr, REAL vo1, REAL vn1, REAL t, REAL ht )
{
 static int first=1;
 static REAL pa1, mdp1, dvmax1, vdvmax1, apd1, cl1, tcl1; 
 REAL ret=-1.;

 for( ;first;first=0)
 {	
   // PA = Peak Amplitude
	fprintf( stderr, "### CL\t APD\t PA\t MDP\t dvmax\t vdvmax\t t\n" );
	tcl1=pa1=dvmax1=vdvmax1= -3e33; 
	mdp1=3e33;
 }//first

 //front
 if( vn1>=vthr && vo1<vthr ) 
 {
   if( tcl1 >0 ) // not first front
   {	
     ret = cl1 =t-tcl1;			// return the CL
     fprintf( stderr, "%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\n", 
	      cl1, apd1, pa1, mdp1, dvmax1, vdvmax1, t);
   }
   tcl1 = t;
   pa1=dvmax1=vdvmax1= -3e33; 
   mdp1=3e33;
 }//front

 //tail
 if( vn1<vthr && vo1>=vthr ) 
 {
   apd1 = t-tcl1;
 }//tail

 // every step
 if( vn1>pa1 ) pa1 = vn1;
 else if( vn1<mdp1 ) mdp1 = vn1;      //with "else" there is a little economy
 if( (vn1-vo1)/ht >dvmax1 ) 
 {
   dvmax1 = (vn1-vo1)/ht;
   vdvmax1 = vo1;
 }
 

 return ret; 		// return value (now CL)
}////basic_pars
