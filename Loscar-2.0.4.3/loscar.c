/****************************************************************
 
   file: loscar.c VERSION: LOSCAR V 2.0.4.3      09/2018

   LOSCAR Model: Long-term Ocean-atmosphere-Sediment 
                 CArbon cycle Reservoir Model

   *** LOSCAR comes with ABSOLUTELY NO WARRANTY ***
   *** Use at your own risk. DO NOT DISTRIBUTE  ***
 
   When using LOSCAR, cite as:

   Zeebe, R. E., Zachos, J. C., and Dickens, G. R. Carbon dioxide
   forcing alone insufficient to explain Paleocene-Eocene Thermal
   Maximum warming. Nature Geoscience, doi:10.1038/NGEO578 (2009) 
   
   Richard E. Zeebe
   School of Ocean and Earth 
   Science and Technology 
   Department of Oceanography 
   University of Hawaii at Manoa
   1000 Pope Road, MSB 504
   Honolulu, HI 96822, USA
   email: loscar.model@gmail.com
	   
   /home/zeebe/GCC/Loscar

   MAKE (modern version): make loscar 
   MAKE (paleo  version): make loscar PALEO=1

   COMPILER:	   
   gcc -Wall -lm -o ...
   other options:
   -Wall -Wextra -Wundef -pedantic 
   Optimize: -O2 -O3

   CODE: should be ANSI C compatible

   USE (linux example):  
   modern version: $./loscar.x preind.inp  >& loscar.log &
   paleo  version: $./loscar.x prepetm.inp >& loscar.log &
		   
   SOLVER: NR
   void stiff(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)  
   Fourth-order Rosenbrock
  
   TO DO: 

       @ mass balances
	   @ splint, SatHor?, removecr?
 	   @ check ERROR SCALING in case of problems. see solver.c, 
	     odeint(). normalize all variables to order 1.

   COUNTERS: require that tfinal > tstart. then:
			 
   int     what                            set by	  range 
   ==================================================================
   nstep   total steps taken (=nok+nbad)   solver 	 1 to MAXSTP 
   kmax    desired No. output values       user      2 to MAXSTP+1
   kount   actual  No. output values       solver    2 to ?
   lt      = kount

   vars are allocated using kmax and include value at tstart.
   thus if min(nstep)=1 => min(kmax)=2.

   NOTE: generally the following convention is used: 
	     TRUE  = ON  = 1 
	     FALSE = OFF = 0
 	   
   NOTE: fcsys(). a2 for kb should be -2.608, not +2.608 
   (thanks to James Rae). use plus when checking old results.
   
   NOTE: vectors initialized using NR start at index 1, not 0!

   NOTE: gcc 4.4.3 linux defines __STDC__   
  	   
   DEBUG/CHECK: 

   printf("\n>>> %.15e <<<\n",x); exit(0);
   grep -n bla *
	   
   test cases 
   V2.0.4.3 modern: [368.33502623] paleo: [1067.37290006] new K*(sij)
   V2.0.4.2 modern: [368.33502623] paleo: [1070.03120666]
   V2.0.3 FTYS, CINP = 1000., TFINAL = 200.e3, O2 mol/m3
	      Final Atm CO2:  1070.02463727 (gcc 4.6.1: same. omit -ansi)
   V2.0.2 FTYS, CINP = 1000., TFINAL = 200.e3
	      Final Atm CO2:  1070.03221410
   V2.0: FSEDCC, kmax=1000, ldrestart=1, ffflag=1,  tsnsflag=0/1.
	     Final Atm CO2:  368.33323766 TSCAL=10 no ts
	     (including empty Tethys equations)
   V1.9: FSEDCC, kmax=1000, ldrestart=1, ffflag=1,  tsnsflag=0/1.
	     Final Atm CO2:  368.33502623 TSCAL=10 no ts
         (including SHLF)
	     Final Atm CO2:  368.33502927 TSCAL=10 no ts
	     (including dox, d13C, new total boron)
	     Final Atm CO2:  368.90552696 TSCAL=10 no ts
	     (including dox, y0preind-1.9.1.dat)
   V1.8: FSED, kmax=1000, ldrestart=1, ffflag=1,  tsnsflag=0/1.
	     Final Atm CO2:  368.92057262 TSCAL=10 no ts
	                     391.11286579 TSCAL=10    ts
	                     391.10400401 TSCAL=1     ts
                         391.04926531 ftsns()
         diff pco2a.dat check/pco2a.dat-1.8ts
   V1.7: FSED, kmax=1000, ldrestart=1, ffflag=1, tsnsflag=1.
	     Final Atm CO2: 390.97608672 
         diff pco2a.dat check/pco2a.dat1.7
   V1.6: FSED
	     Final Atm CO2: 368.90458443 
         diff pco2.dat check/pco2.dat1.6
   V1.3: FSEDU, NOCT=3,NB=10,NCATM=1,t=[0,5000],yp: +dic, pco2,
	     epslvr=1e-4, kmax=500, no optimization. 
	     Final Atm CO2: 286.05596348 a2-kb=-2.608
                        286.04585891 a2-kb=+2.608
         diff dic.dat check/dic.dat1.3
   V1.1: NOCT=3,NB=10,tfinal=5000,PO40,yp: bio+thmfun,
	     epslvr=1e-6, kmax=100, no optimization.
         diff po4.dat check/po4.dat1.1 
   V1.0: NOCT=2,NB=10,tfinal=5,y0=i,yp=-y.
         diff dic.dat check/dic.dat1.0 

	   
   updates:

   09/03/18 new K* corrections (ZeebeTyrrell18)
   12/06/17 updated sediment routines. see initstart, derivs
   11/30/11 revised double comparison
            see Makefile. gcc 4.6.1 -ansi sets constant=long double
   10/31/11 CLIMSENS -> sclim (as input)
   10/28/11 csys(): O2 solubility: ml/kg => ml/l	   
   10/26/11 PE 13C input values. new y0prepetm-2.0.3.dat
   10/25/11 co3satm[aipt] all basins
   10/23/11 Tethys done. first steady-state tests OK
   10/17/11 include Tethys V 2.0
   10/15/11 include shelf/deep rain
   10/08/11 add Ca,Mg to fcsys()
   09/23/11 include sediment C13	   
   09/20/11 fkrg: 9>10. new y0preind.dat
   09/18/11 Total boron: 416.0 DOE94 => 432.5 Lee10.
		    new y0preind.dat
   09/16/11 include 13C (dicc).
   09/15/11 include dissolved O2 (dox) + MM kinetics. V 1.9
   04/20/11 added MM kinetics HighLat PO4.	   
   04/19/11 added fdapi (change initial dic, alk, po4).
            added warning for high ALK/DIC ratio 
            (Follow's CO2 routine gives negative [H+])
   04/12/11 include temperature as ocean tracer (tcb). temp change
		    (CO2-sensitivity) much easier and numerically stable.
			LOSCAR V 1.8		    
   04/10/11 checked and cleaned up global vars (a few renamed)
   04/05/11 included temperature change (CO2-sensitivity) via
	        temperature update. => a few issues with numerics.
			LOSCAR V 1.7
   04/03/11 calc and write ccd. matches Loscar.m
		    (interp. lin., fcsml=0.20, dz=10 m)
		    replaced femiss() by finterp().
			LOSCAR V 1.6
   04/02/11	calc and write pco2ocn, pH, CO3, omega.		
   03/26/11 added README.txt. a2 for kb set to -2.608,
            new y0preind.dat
   03/23/11 read file names from restart file
            (restart, save, emissions).
   03/20/11 read and save restart values OK.
   03/19/11 read emissions OK. fcsml=0.05.
  	        LOSCAR V 1.5
   03/16/11 sediments OK. incl. dic,alk,po4,pco2.
   03/13/11 1st sediment test OK. dic,alk,po4=const.
		    matches Loscar.m
  	        LOSCAR V 1.4
   03/06/11 dic, pco2 OK. output formatted. set epsj=epslvr. 
  	        LOSCAR V 1.3
   03/04/11 readparms() OK.
  	        LOSCAR V 1.2
   03/02/11 function prototypes and defs => ANSI.  
   03/01/11 code check: splint -mustfreefresh 
		    -nullpass (+bounds)
			'+= +(...)' gives *** Internal Bug at 
			constraintGeneration.c:1441. (minus OK)
   02/28/11 V 1.1 tested on LNX32/64, MAC, XP. 
   02/27/11 csys() OK. ALK 10-box done. OK.
  	        LOSCAR V 1.1	   
   02/21/11 PO4 10-box done. matches Loscar.m. OK.
   02/19/11 V 1.0 tested on
            Linux,  i486,    gcc 4.4.3
            Linux,   x86_64, gcc 4.4.5 -m32 compiler option
            Mac,    i686_64, gcc 4.0.1
            Win XP,  x86,    gcc 3.4.2, +GNU make for Win
   02/19/11 common.h: global vars.
   02/18/11 thmfun() OK.
   02/16/11 TO DO: PO4 10-box, csys(), thmfun().
   02/15/11 writedat() works. PO4 3-box OK.
	        LOSCAR V 1.0
   02/13/11 new file. (based on loscar.m, m-file started in 2005)

*****************************************************************/
#include <stdio.h>
#include <time.h>
#include "defs.h"
#include "utils.h"

#include "common.h"

#ifndef TIMESTEP
  #define TIMESTEP 5000
#endif 

int index;

/*============================================================*/
/*====== function declarations and global variables ==========*/
/*============================================================*/
/* see defs.h, utils.h, and common.h                          */


/*============================================================*/
/*============================================================*/
/*==================== int main() ============================*/
/*============================================================*/
/*============================================================*/
int main(int argc, char **argv)
{
 int nok,nbad,lt;

 /* time vars */
 clock_t tick_start, tick_end;
 double wall_time;

 /* default: initialize parameters/variables */
 initfree(1);

 /* optional: read control parameters/variables */
 if(argc == 2){
    cntrfflag = 1;
    readparms(argv[1]);
 }
 else{
    cntrfflag = 0;
    printf("\n@ Setting up run. Using default parameters\n");
 }

 /* check variables and allocate variables depending on */ 
 /* parameters that change if control file is read      */
 initfree(2);
	
 /* optional: run fossil fuel scenario, load emissions */ 
 if(ffflag == 1){
   printf("\n@ Loading emissions: '%s'\n",ffldstr);
   reademiss();
 }
 else{
   printf("\n@ Emissions loaded: none\n");
 }	

 /* initialize y-start values (default or load) */
 initstart();
	
 /* print start parameters/variables */
 initfree(3);

 /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
 /*%%%%%%%%%%%%%%%%%%%%%%%% solver %%%%%%%%%%%%%%%%%%%%%%%%%%%*/
 /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

 printf("\n@ Starting integration");
 printf("[tstart tfinal]=[%.2e %.2e]\n",t0,tfinal);
 printf("  Dum spiro spero.\n");		
 tick_start = clock(); /* start clock */
	
 odeint(ystart,NEQ,t0,tfinal,epslvr,h1slv,hminslv,&nok,&nbad,derivs,stiff);

 printf("\n@ Done integrating\n");
 printf("  Alea iacta est.\n");		
 tick_end  = clock(); /* end clock */
 wall_time = (double)(tick_end-tick_start)/CLOCKS_PER_SEC;
 printf("\n@ Wallclock Time = %f sec \n",wall_time);
	
 /* print step info */
 printf("\n@ Step Info:\n");
 printf("     %6d steps OK\n",nok);
 printf("     %6d steps FIXED\n",nbad);
 printf("     %6d steps TOTAL\n",nok+nbad);
 printf("  kmax=%d, kount=%d\n\n",kmax,kount);
 if(kount < 2)
    ferrwrt("main(): kount < 2. increase tfinal?");
 lt = kount;
	
 printf("@ [tfinal      solution[1][lt]]\n");
 printf("  [%.2e  ",tmv[lt]);
 printf("%12.8f] \n\n",yy[1][lt]);
 if(NCATM == 1){
   printf("@ Final Atm CO2: ");
   printf("%12.8f \n",yy[NOCT*NB+1][lt]/PPMTOMMSQ/CNTI);
 }

 /* write data to files */
 printf("\n@ Checking and writing output (.dat) files\n");
 writedat(tmv,yy,lt,tcb0,spm,cac,mgc,s4c,dsv,zv,nz,svrestart,fpsvstr);

 /* free parameters/variables */
 printf("\n@ Freeing memory\n");
 initfree(4);

 return(0);	
}
/*============================================================*/
/*==================== int main END ==========================*/
/*============================================================*/

