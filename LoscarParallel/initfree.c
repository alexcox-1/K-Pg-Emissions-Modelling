/* ac32 changed nsi to 0.60 */

/****************************************************************
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
	   
*****************************************************************/
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "defs.h"
#include "utils.h"

#include "common.h"


/*============================================================*/
/*==================== initfree() ============================*/
/*============================================================*/
/* initialize, allocate, print, and free global parameter/variables      
   flag = 1,2,3,4

   1: initialize and allocate
   2: alloc vars that may change if control file is read
   3: print
   4: free

 NOTE: helper arrays (say x) holding a list of numbers may be 
       intialized starting at index 0. x[0] is set to 0 and 
	   is not used.
   updates:

   09/03/18 new K* corrections (ZeebeTyrrell18)
   11/30/11 revised double comparison (excess precision)		   
   10/21/11 init Tethys variables
   10/18/11 free_dvector(gp,1,NOC=>NB)
   10/17/11 include Tethys V 2.0   
   09/20/11 fkrg: 9>10. new y0preind.dat		   
   04/19/11 added fdapi (change initial dic, alk, po4).		   
   04/15/11 temperature: added tcb as tracer (rm tcv). start priorities:
            1. control file 2. restart file 3. intern default
            NOTE: restart is read after cntrl => remember tcb0c
   04/09/11 added flag: 
            check variables and allocate variables depending on 
            parameters that change if control file is read  
   02/18/11 new file
			   
 */
void initfree(int flag)
{
 int i,k;
 char mssg[BUFSIZ];
 double d13cvc,d13cin,d13ckrg,kspv[3];
	
 FILE *fpout;

 time_t timer;

 /*==== defs: area, height of ocean boxes. ==========*/	 
#ifdef FTYS /* TETHYS */
 /* area fractions         A    I   P   T  */
 double fanoc[NOC+1] = {0,0.15,.14,.52,0.09};
 /* deep volume fractions */	
 double fdvol[3+1]   = {0,0.16,.16,.68};
 double vres;	
#else
 /* area fractions         A    I   P      */
 double fanoc[NOC+1] = {0,0.26,.18,.46};
#endif
 double fhl = 0.1; /* H */	
	
 /* height (m)          L    I     D */
 double hlid[3+1] = {0,100.,900.,(HAV-1000.)};
	
 /*==== defs: temperature ======================*/	 
  /* temp (deg C)     L   I   D  */
#ifdef FTYS /* TETHYS */
 double tc3[3+1] = {0,25.,16.,12.}; /* AIP    */
 double tct[3+1] = {0,18.,14.,12.}; /* Tethys */
#else
 double tc3[3+1] = {0,20.,10.,2.};  /* AIP    */
#endif	
	
 /*==== defs: mixing parameters ====================*/	 
 /* mixing              A     I     P    TLI    TII */
#ifdef FTYS /* TETHYS */
 double  mv0[KOC+1]={0,3.5e6,3.5e6,7.0e6,3.2e6,2.e6}; /* 3.5 3.5 7.0 3.2 2.*/
 double mhd0[NOC+1]={0,4.0e6,4.0e6,6.e6,0.7e6};       /* 4.0 4.0 6.0 0.7   */	 
#else
 double  mv0[KOC+1]={0,5.5e6,4.5e6,6.5e6}; /* 5.5 4.5 6.5 */
 double mhd0[NOC+1]={0,3.0e6,2.0e6,8.0e6}; /* 3.0 2.0 8.0 */	 
#endif
	
 /*==== defs:fraction EPH, remineralized in deep A,I,P boxes */
 /*                    A   I   P                             */
 double gp0[NOC+1]={0,0.3,0.3,0.4};

#ifdef FSED
	
 int *indxv;
 double gam,dz;
	
 /*==== sediment box depths                                  */
 /* (km) => m, see below                                     */
 double dsv0[NSD+1]={0,.1,.6,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.5};

#ifdef FTYS /* TETHYS */
 /*==== area fractions A,I,P,T  (2x2 deg, Karen Bice)        */
 /* percent => frac, see below                               */
 double asva0[NSD+1]={0,1.1407,10.0216,8.7160,6.3724,4.7915,3.4973,
	          12.2935,10.7438,8.4983,11.4134,11.0948,11.4167,1.e-6};
 double asvi0[NSD+1]={0,0.2501,5.5345,5.5145,8.9550,4.6283,6.5361,
	           11.7221,12.6050,14.6295,13.8384,8.0807,7.7057,1.e-6};
 double asvp0[NSD+1]={0,0.1673,2.8333,3.0599,2.5389,1.4218,5.0153,
	          9.8023,14.0117,10.1975,20.0019,13.7155,17.2346,1.e-6};
 double asvt0[NSD+1]={0,7.0534,46.5363,22.4068,7.1501,2.4261,3.9946,
	               2.7063,0.8532,0.9814,3.0802,2.3583,0.4533,1.e-6};
#else	
 /*==== area fractions A,I,P   (Menard & Smith, 1966)        */
 /* percent => frac, see below                               */
 double asva0[NSD+1]={0,7.0297,5.1729,1.9106,2.3882,4.2988,4.2988,  
	         9.6711,9.6711,16.2389,16.2389,11.1712,11.1712,0.7388};
 double asvi0[NSD+1]={0,3.5710,2.6844,1.5907,1.9884,5.0146,5.0146,
	         12.6299,12.6299,18.3221,18.3221,8.4957,8.4957,1.2407};
 double asvp0[NSD+1]={0,1.6358,2.5901,1.4484,1.8105,3.4372,3.4372,
	       10.9275,10.9275,17.5411,17.5411,13.4784,13.4784,1.7468};
#endif /* FTYS */
#endif /* FSED */
	
 /*============================================================*/
 /*============================================================*/
 /* initialize */	
 if(flag == 1){

 timer=time(NULL);
 printf("%s\n",asctime(localtime(&timer)));
 printf("  Heus!\n\n");	 
 printf("   ************************************************\n");	
 printf("   * LOSCAR V %s                  %2d/%4d      *\n",
        VLOSCAR,VMONTH,VYEAR);	
 printf("   *                                              *\n");		 
 printf("   * Richard E. Zeebe                             *\n");		 
 printf("   * SOEST, University of Hawaii                  *\n");		 
 printf("   * Honolulu, HI 96822, USA                      *\n");		 
 printf("   * loscar.model@gmail.com                       *\n");		 
 printf("   *                                              *\n");		 
 printf("   * LOSCAR: Long-term Ocean-atmosphere-Sediment  *\n");		 
 printf("   * CArbon cycle Reservoir Model                 *\n");		 
 printf("   *                                              *\n");		 
 printf("   * LOSCAR comes with ABSOLUTELY NO WARRANTY.    *\n");		 
 printf("   * Use at your own risk. DO NOT DISTRIBUTE.     *\n");		 
 printf("   *                                              *\n");		 
 printf("   ************************************************\n");	

#ifdef FTYS
 printf("\n@ This is the PALEO setup including Tethys\n");
#else
 printf("\n@ This is the MODERN setup\n");
#endif
	 
 printf("\n@ No. ocean basins         : %3d",NOC);
 printf("\n@ No. ocean boxes          : %3d",NB);
 printf("\n@ No. ocean tracers        : %3d",NOCT);
 printf("\n@ Atmospheric Carbon       : %3d",NCATM);
 printf("\n@ Atmospheric Carbon-13    : %3d",NCCATM);
 printf("\n@ No. sediment depth levels: %3d",NSD);
 printf("\n@ No. sediment C-13 levels : %3d",NSDCC);
 printf("\n@ No. equations            : %3d",NEQ);
 printf("\n");
	 
	 
 /* @@@@@@@@@@@@@@@@@@@@@@ CAUTION @@@@@@@@@@@@@@@@@@@@@@@@
    when using global variables (see common.h), no warnings
    are issued if these are used uninititalized. 
	So the order matters! check values in files:
    parms.out, ystart.out
    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */ 


 /*** default flags: superseded by parameter input file!  ***/

 /* load restart (ystart) ON/OFF = 1/0 */
 ldrestart = 1;
 /* save restart (yfinal) ON/OFF = 1/0 */
 svrestart = 0;
 /* fossil fuel flag, load emissions ON/OFF = 1/0 */
 ffflag    = 0; /* 0 default (runtest) */
 /* OCN temperature change (CO2-sensitivity) ON/OFF = 1/0  */
 tsnsflag  = 0;
 /* carbon input ON/OFF = 1/0 */	
 cinpflag  = 0;

	 
 /* default file names: superseded by parameter input file!  */
 strcpy(fpldstr,"dat/y0preind-2.0.4.dat");   /* load restart   */	 
 strcpy(fpsvstr,"dat/y0preind-2.0.x.dat");   /* save restart   */
 strcpy(ffldstr,"dat/Emss/1000_0500.dat"); /* load emissions */

#ifdef FTYS
 ldrestart = 1;	
 ffflag    = 0;	 
 cinpflag  = 0; /* 0 default (runtest) */
 strcpy(fpldstr,"dat/y0prepetm-2.0.4.dat"); /* load restart   */
 strcpy(fpsvstr,"dat/y0prepetm-2.0.x.dat"); /* save restart   */
#endif	 

 /* default OCN temperature change (CO2-sensitivity) */
 sclim = 3.0;
	 
 /* conveyor transport */
 /* fconv: 1 = NADW, 2 = NPDW, 3 = SO */	 
#ifdef FTYS /* TETHYS */
 fconv = 3;            /* 3 */
 thc0  = 25.e6*YTOSEC; /* (Sv -> m3/y) 25 */	 
 tto   = 02.e6*YTOSEC; /* (Sv -> m3/y)  2 */;	
#else
 fconv = 1;            /* 1 */
 thc0  = 20.e6*YTOSEC; /* (Sv -> m3/y) 20 */
 tto   = 0.0;	 
#endif
 thc = thc0;
 /* PETM: add SO to NPDW, see loscarderivs() */ 
 tso = 0.0; /* default = 0.0 */
	 
 /* TH branches */
 thbra = 0.20;  /* 0.20 upwelled into intermdt Atl */
 thbri = 0.20;  /* 0.20 upwelled into intermdt Ind */

#ifdef FTYS /* TETHYS */
 t0     = 0.0   ; /* (y) interval: tstart  0.     */ 
 tfinal = 200.e3; /* (y) interval: tend    200.e3 */	 
#else	 
 t0     = 1700.00; /* (y) interval: tstart  1700. */ 
 tfinal = 3000.00; /* (y) interval: tend    3000. */	 
#endif	 
 /* Solver options */
 epslvr  = 1.e-04;            /* accuracy 1.e-04 (sloppy) */ 
 kmax    = 1000;              /* max # output values      */
 dxsav   = (tfinal-t0)/kmax;  /* approx output interval   */
 h1slv   = 0.01;              /* 1st guess step size      */
 hminslv = 0.00;              /* min       step size      */
	 
 /* fcsml, numerics: linear f_calcite drop  */
 /*	during dissolution if fc < fcsml. 0.05  */
 fcsml = 0.05;	 

 /* allocate globals (more allocation, see flag == 2)  */	 
 ystart = dvector(1,NEQ);        /* ystart             */
 vb     = dvector(1,NB);         /* volume boxes       */
 ab     = dvector(1,NB);         /* areas              */
 hb     = dvector(1,NB);         /* height boxes       */
 kkv    = ivector(1,NTS);        /* Surface: Low-lat(NLS) + High(1) */ 
 kiv    = ivector(1,NOC);        /* Interm             */ 
 mxv    = dvector(1,KOC);        /* mixing AIP TLI TII */
 mhd    = dvector(1,NOC);        /* mixing H-AIP T     */
 gp     = dvector(1,NB);         /* H-frac remineralized in AIP */
 hgssv  = dvector(1,NB);         /* csys: H+ guess     */
 tcb0   = dvector(1,NB);	     /* temperature boxes init       */
 tcb0c  = dvector(1,NB);	     /* temperature boxes init cntrl */
 salv   = dvector(1,NB);         /* salinity    boxes  */
 prsv   = dvector(1,NB);         /* pressure    boxes  */
 spm    = dmatrix(1,NB,1,2);     /* SP matrix   boxes  */
 fdapi  = dvector(1,3);          /* factor initial dic, alk, po4 */
	 

 /* Boxes: Atlantic, Indic, Pacific, High, Tethys

	  A I P           A I P         A I P          H         L   I  D
 Low (1,2,3), Interm (4,5,6), Deep (7,8,9), High (10), Teth (11,12,13)

 
 AL   1      
 IL   2    
 PL   3     
 AI   4    
 II   5     
 PI   6     
 AD   7   
 ID   8   
 PD   9   
 H   10
 TL  11
 TI  12
 TD  13
	 
 */


 /*==== area, height, volume of ocean boxes. =======*/	 
 /* (a bit complicated but gives exact ocean volume
	for given area fractions and LI box heights)    */

#ifdef FTYS	 
 if(NB != 13){
    fprintf(stderr,"\nNB = %d",NB); 
    ferrx("intifree(): set NB (#boxes) to 13");
 }	
#else
 if(NB != 10){
    fprintf(stderr,"\nNB = %d",NB); 
    ferrx("intifree(): set NB (#boxes) to 10");
 }	
#endif
	 
 /* area fractions          A    I   P   T
 double fanoc[NOC+1] = {0,0.26,.18,.46}; 
 double fanoc[NOC+1] = {0,0.15,.14,.52,0.09};
 */
 for(k=1;k<=3;k++){
       ab[k]   = fanoc[k]*AOC;	
       ab[k+3] = fanoc[k]*AOC;	
       ab[k+6] = fanoc[k]*AOC;	
 } 
 ab[10] = fhl*AOC;	/* HighLat fraction 0.1 */

#ifdef FTYS /* TETHYS */
 for(k=1;k<=3;k++)
       ab[k+10] = fanoc[4]*AOC;	 /* Tethys */		 
#endif	 
	 
 /* height (m)          L    I     D 
 double hlid[3+1] = {0,100.,900.,(HAV-1000.)}; */
 for(k=1;k<=3;k++){
       hb[k]   = hlid[1];	
       hb[k+3] = hlid[2];	
       hb[k+6] = hlid[3];	
 } 
 hb[10] = 250.;	

#ifdef FTYS /* TETHYS */
 for(k=1;k<=2;k++)
       hb[k+10] = hlid[k];	/* Tethys L I */		 
 hb[13] = 200.;	            /* Tethys D   */
#endif	

 /*=== calculate final volume and height ===*/
#ifdef FTYS /* TETHYS */
 /* residual volume */ 
 vres = VOC-(hlid[1]+hlid[2])*(1.-fhl)*AOC
		   -hb[10]*ab[10]-hb[13]*ab[13];
	 
 /* distribute into deep AIP */	 
 for(k=1;k<=3;k++)
     hb[k+6] = vres*fdvol[k]/ab[k+6];
	 
 /* volume */	 
 vvelm(vb,ab,hb,NB); 
#else	 
 /* volume */
 vvelm(vb,ab,hb,NB); 
	 
 /* now add to deep boxes: the
    volume below H box = ab[10]*(HAV-hb[10]) */
 for(k=1;k<=3;k++)
	 vb[k+6] += ab[10]*(HAV-hb[10])/3.;

 /* recalculate height */
 for(k=1;k<=NB;k++)
       hb[k] = vb[k]/ab[k];	
#endif /* FTYS */

 /* set box indices */	
 /* kkv: surface    */	 
 /* kiv: interm     */	 
 for(k=1;k<=3;k++){
      kkv[k] = k;   /* AIP */
	  kiv[k] = k+3; /* AIP */
 }
 kkv[4] = 10;	    /* H   */
#ifdef FTYS
 kkv[5] = 11;	 /* Tethys */
 kiv[4] = 12;	 /* Tethys */
#endif 	 
	 

 /*============ temp, sal, pressure =============*/	 
 /* set internal default                         */
 /* temp (deg C)      L   I   D 
 double tc3[3+1] = {0,20.,10.,2.};	*/ 
 for(k=1;k<=3;k++){
       tcb0[k]   = tc3[1];	
       tcb0[k+3] = tc3[2];	
       tcb0[k+6] = tc3[3];	
 } 
 tcb0[10] = 2.;	

#ifdef FTYS /* TETHYS */
 for(k=1;k<=3;k++)
       tcb0[k+10] = tct[k];			 
 tcb0[10] = 12.;	
#endif
	 
 /* NOTE: start temperature (ystart) is set in initstart() */

 /* salinity */	 
 for(k=1;k<=NB;k++)
       salv[k] = SALOCN;

 /* pressure (bars) */	 
 for(k=1;k<=3;k++){ /* mean box height (*0.5), => bars (*0.1) */
       prsv[k]   =  0.5*hb[k]                 *0.1;  /* surf  */
       prsv[k+3] = (0.5*hb[k+3]+hb[k]        )*0.1;  /* intrm */
       prsv[k+6] = (0.5*hb[k+6]+hb[k]+hb[k+3])*0.1;  /* deep  */
 }
 prsv[10] = 0.5*hb[10]*0.1; /* high */

#ifdef FTYS /* TETHYS */
 k = 11;	 
 prsv[k]   =  0.5*hb[k]                 *0.1;  /* surf  */
 prsv[k+1] = (0.5*hb[k+1]+hb[k]        )*0.1;  /* intrm */
 prsv[k+2] = (0.5*hb[k+2]+hb[k]+hb[k+1])*0.1;  /* deep  */
#endif	 
	 
 /* store SP in matrix */
 fvec2arr(spm,salv,prsv,NB);

 /* set Ca,Mg */	 
#ifdef FTYS /* TETHYS */
 cac = 20.0e-3; /* mol/kg 20.0e-3 */
 mgc = 30.0e-3; /* mol/kg 30.0e-3 */
 s4c = 14.0e-3; /* mol/kg 14.0e-3 */	 
#else
 cac = CAM;
 mgc = MGM;
 s4c = S4M;	 
#endif
 /* get Ksp ratio = Ksp(past)/Ksp(modern) */
#ifdef KCORR18
 kspfun(&kspv[1],&kspv[0],25.,35.,0.,CAM,MGM,S4M);		 
 kspfun(&kspv[2],&kspv[0],25.,35.,0.,cac,mgc,s4c);		 
 rksp = kspv[2]/kspv[1];
#endif	 
	 
 /*================ mixing parameters ===========*/	 

 /* mixing (modern)      A     I     P          
 double  mv0[KOC+1]={0,5.5e6,4.5e6,6.5e6}; 
 double mhd0[NOC+1]={0,3.0e6,2.0e6,8.0e6};  */	 
	 
 vsclr(mxv,(3.8*YTOSEC),mv0 ,KOC); /* Sv -> (m3/y) 3.8: tuning */
 vsclr(mhd,(1.3*YTOSEC),mhd0,NOC); /* Sv -> (m3/y) 1.3: tuning */


 /* init fdapi. change initial dic, alk, po4, see initstart(). */ 
 /* only control file modifies fdapi. init default: 00.00%     */
 for(k=1;k<=3;k++)
     fdapi[k] = 00.00; /* percent */

	 
 /*================== Biological Pump ===========*/
 eph   = 1.8*ab[10];/* (mol/y) 1.8 H Export, mol C/m2/y*A = mol/y  */
#ifdef FTYS /* TETHYS */
 rrain = 6.7;       /* 6.7 export rain ratio (Corg/CaCO3)          */
#else
 rrain = 6.1;       /* 6.1 export rain ratio (Corg/CaCO3)          */
#endif
#ifdef FSED
 nuwd  = 0.31;      /* 0.31 water column dissolution               */
#else
 nuwd  = 0.0;    
#endif
	 
 fepl  = 0.80;      /* 0.80 LL utilization   */
 frei  = 0.78;      /* 0.78 fraction EPL, remineralized in I boxes */


 /* fraction EPH, remineralized in deep A,I,P boxes */
 /*                    A   I   P                        
 double gp0[NOC+1]={0,0.3,0.3,0.4};                     */
 fsetv(gp,0.0,NB);
 for(k=1;k<=3;k++)
     gp[k+6] = gp0[k];

 /*========== silicate weathering: volc degass ===*/
 nsi   = 0.60;       /* CaSiO3 weath exponent 0.20        */
#ifdef FTYS /* TETHYS */
 pcsi  = 1000.;      /* uatm, std-stt atm pCO2 1000       */
 fvc0  = 05.e12/AOC; /* mol C, degassing /m2/y @1000 uatm */
 fvc0 *= pow(2.,nsi);
#else	 
 pcsi  = 280.;       /* uatm, std-stt atm pCO2 280        */
 fvc0  = 05.e12/AOC; /* mol C, degassing /m2/y @280 uatm  */
#endif	 

 /* Note: fvc0 is steady-state degassing @280 uatm, 
    balanced by steady-state Si weathering @280 uatm.
	Change fvc(t) in loscaderivs(), not here.

 */
	 	 
 /*================ CaCO3 in-flux ================*/
 ncc   = 0.40;           /* CaCO3 weath exponent 0.40    */
#ifdef FTYS /* TETHYS */ /* mol C    /m2/y riverine flux */
 finc0 = 15.83409493e12/AOC;
 /*finc0 = 1.0*12.e12/AOC; 
   finc0 *= pow(2.,ncc); */
#else
 finc0 = 1.0*12.e12/AOC; /* mol C    /m2/y riverine flux */
#endif

 /* Some of the values for the C-Cycle parameters 
	below are only loosely constrained, some have 
	been be tuned within reasonable ranges. A few 
    references that provide ranges are given in 
	brackets:
	  
    WK92 = WalkerKasting92, Berner = GEOCARB 123
    KA99 = KumpArthur99, Veizer = VeizerEtAl99
	Hayes = HayesEtAl99.	  
 */
	 
 /*============== kerogen oxidation ==============*/ 
#ifdef FTYS /* TETHYS */
 fkrg  = 07.e12/AOC;     /* mol C /m2/y 07 [WK92,Berner] */
#else
 fkrg  = 10.e12/AOC;     /* mol C /m2/y 10 [WK92,Berner] */
#endif
	 
 /*=========== Long C-Cycle fluxes 13C ===========*/
#ifdef FTYS /* TETHYS */
 d13cvc  =  -4.0; /* -4.0  degas [Berner,WK92,KA99]           */    
 d13cin  =   2.0; /*  2.0  CaCO3 [Berner,Veizer]              */
 d13ckrg = -22.2; /* -23.2 kerogen [WK92/tuned]               */
 epscorg = -33.0; /* -33.0 eps(Corg-DIC) [tuned+Berner,Hayes] */
#else
 d13cvc  =  -4.0; /* -4.0  degas [Berner,WK92,KA99]           */   
 d13cin  =   1.5; /*  1.5  CaCO3 [Berner,Veizer]              */
 d13ckrg = -21.0; /* -21.0 kerogen [WK92/tuned]               */	
 epscorg = -27.7; /* -27.7 eps(Corg-DIC) [tuned]              */
#endif
 rvccc   =  (d13cvc/1.e3+1.)*RST;
 rincc   =  (d13cin/1.e3+1.)*RST;
 rkrgcc  = (d13ckrg/1.e3+1.)*RST;


 /* Carbon input default values */
 cinp   = 1000.;	 
 dccinp = -55.;	 
 rccinp = (dccinp/1.e3+1.)*RST;
 tcin0  = 0.;
 tcspan = 6.e3;
	 
 /* csys: first guess for H+  */
 for(k=1;k<=NB;k++)	 
    hgssv[k] = HGUESS;
	 

#ifdef FSED

 /* allocate globals */	 
 dsv     = dvector(1,NSD);
 asva    = dvector(1,NSD);     
 asvi    = dvector(1,NSD);     
 asvp    = dvector(1,NSD);     
 klid    = ivector(1,NSD); /*  sed box index */
 nlid    = ivector(1,3);   /* #sed boxes LID */
 fc0a    = dvector(1,NSD);
 fc0i    = dvector(1,NSD);
 fc0p    = dvector(1,NSD);
 phiia   = dvector(1,NSD);     
 phiii   = dvector(1,NSD);     
 phiip   = dvector(1,NSD);
 #ifdef FTYS
 asvt    = dvector(1,NSD);     
 klidt   = ivector(1,NSD); /*  sed box index Tethys */
 fc0t    = dvector(1,NSD);
 phiit   = dvector(1,NSD);     
 #endif	 
 #ifdef FSEDCC
 fcc0a   = dvector(1,NSD);
 fcc0i   = dvector(1,NSD);
 fcc0p   = dvector(1,NSD);
  #ifdef FTYS
  fcc0t   = dvector(1,NSD);
  #endif	 	 
 #endif	/* FSEDCC */ 

 /* rain of 'remainder' */
 frrf  = 0.35;             /*  g/cm2/ky remainder */
 frrf *= 1.e4/1.e3/1.e3;   /* -> kg/ m2/ y        */
	 
	 
 /*================ sediments ====================*/
 /* dissolution parameter                         */
 /* ~fc^0.5*(cs-c)     Kd defnd in DEQ            */
 ncsd = 2.40;      /* 2.40 calc dissolution order */
 kssd = 20.36e10;  /* mol/m2/y                    */
      	 
 /*====== shelf/deep rain */
#ifdef FTYS /* TETHYS */
 fsh  = 4.5;       /* 4.5 change shelf rain       */
 nsht = 0.4;       /* 0.4 Tethys exponent         */
#else
 fsh  = 1.0;       /* 1.0 change shelf rain       */
#endif	 

 /*====== sed box depths  */	 
 vsclr(dsv,1.e3,dsv0,NSD);           /* km => m   */

 /* depth variable (0 <= z <= dsv[NSD])           */
 nz = 6500+1; /* #grid points = #intervals + 1    */
 zv = dvector(1,nz);
 dz = (dsv[NSD])/((double)(nz-1));
 for(k=1;k<=nz;k++)
	 zv[k] = (double)(k-1)*dz;
	 
 /* area fraction A I P	  */
 vsclr(asva,0.01,asva0,NSD);         /* % => frac */
 vsclr(asvi,0.01,asvi0,NSD);         /* % => frac */
 vsclr(asvp,0.01,asvp0,NSD);         /* % => frac */
#ifdef FTYS /* TETHYS */
 vsclr(asvt,0.01,asvt0,NSD);         /* % => frac */
#endif	 
 /* klid: assign sediment to ocean boxes   */
 /*       Low, Interm, or Deep             */
 fsetiv(nlid,0,3); /* # of sed boxes LID   */
 indxv   = ivector(1,NSD);
 /* interm first (default) */	 
 for(k=1;k<=NSD;k++)  klid[k]  = 4; /* 4 = Intm (Atl) */
 /* low */	 
 dfind(indxv,dsv,NSD,SMALLEQ,hlid[1],MZBH); 
 for(k=1;k<=NSD;k++){ 
    if(indxv[k] != 0){
		              klid[k]  = 1; /* 1 = Surf (Atl) */
                      nlid[1] += 1;
	}
 }	 
 /* deep */	 
 dfind(indxv,dsv,NSD,LARGER,hlid[1]+hlid[2],MZBH); 
 for(k=1;k<=NSD;k++){ 
    if(indxv[k] != 0){
		              klid[k]  = 7; /* 7 = Deep (Atl) */
                      nlid[3] += 1;
	}
 }	 
 nlid[2] = NSD - nlid[1] - nlid[3];

#ifdef FTYS /* TETHYS */
 /* klidt: assign sediment to ocean boxes   */
 fsetiv(klidt,0,NSD);
 for(k=1;k<=nlid[1];k++)
		 klidt[k] = 11;
 for(k=nlid[1]+1;k<=nlid[1]+nlid[2];k++)
		 klidt[k] = 12;
 for(k=nlid[1]+nlid[2]+1;k<=NSD;k++)
		 klidt[k] = 13;
#endif
	 
 /*=================== Porosity ==================*/
 phi0    = 0.85;      /* porosity max   0.85      */
 gam     = 0.23;      /* porosity Delta 0.23      */
 phi1    = phi0-gam;  /* porosity min             */
 
 hsl  = 0.08;         /* (m ) bioturbated layer 0.08 */

 free_ivector(indxv,1,NSD);	 
#endif /* FSED */	 

 } /* END flag == 1 */


 /*============================================================*/
 /*============================================================*/
 /* more allocation */	
 if(flag == 2){
 /* check variables and allocate variables depending on */ 
 /* parameters that change if control file is read      */

 /* first: check a few variables */	 
 if((tfinal - t0) <= 1.e-8){
    fprintf(stderr,"\n[tstart tfinal]=[%e %e]",t0,tfinal); 
    ferrx("intifree(): Set tfinal > tstart + 1.e-8");
 }	

 if((epslvr - 1.e-4) > DEPS){
    fprintf(stderr,"\nepslvr = %e",epslvr); 
    fwarn("intifree(): Your accuracy (EPSLV) is sloppy!");
 }	

 if((fcsml - 0.01) < DEPS){
    fprintf(stderr,"\nfcsml = %e",fcsml); 
    fwarn("intifree(): Your FCSML value is quite small!");
 }	

 if(kmax < 2 || kmax > (MAXSTP+1)){
    fprintf(stderr,"\nkmax = %d",kmax); 
    sprintf(mssg,"intifree(): Reset kmax. range: 2 to %d",MAXSTP+1);
    ferrx(mssg);
 }	

 if((pcsi - 200. + DEPS) < 0.){
    fprintf(stderr,"\npcsi = %e",pcsi); 
    fwarn("intifree(): Your PCO2SI value is quite small!");
 }	
	 
 if(fepl < (0.70-DEPS) || fepl > (0.95+DEPS)){
    fprintf(stderr,"\nfepl = %e",fepl); 
    fwarn("intifree(): Your FBIOL value is out of range.");
 }	
	 
 if(fsh < (1.0-DEPS) || fsh > (6.0+DEPS)){
    fprintf(stderr,"\nfsh = %e",fsh); 
    fwarn("intifree(): Your FSHLF value is out of range.");
 }	

 if(cac <= 1.e-12 || mgc <= 1.e-12 || s4c <= 1.e-12){
    fprintf(stderr,"\ncac = %e",cac); 
    ferrx("intifree(): CALC and/or MAGN and/or SULF too small. Increase.");
 }	

 tmv    = dvector(1,kmax);       /* time vector       */
 yy     = dmatrix(1,NEQ,1,kmax); /* solution matrix   */

 /*=========== recalc changed/derived variables ==============*/

 if(cntrfflag == 1){	 
    /* store initial temperature */	 
    for(k=1;k<=NB;k++)
	    tcb0[k] = tcb0c[k];
 }
 /* store SP in matrix */
 fvec2arr(spm,salv,prsv,NB);	 
 /* approx output interval */	 
 dxsav  = (tfinal-t0)/kmax;  
 
 } /* END flag == 2 */
	
 /*============================================================*/
 /*============================================================*/
 /* print */	
 if(flag == 3){
 /* check errors and write initial vars/parameters to file */	 

 /* check ALK/DIC ratio */
 for(k=1;k<=NB;k++){
   if(ystart[k+1*NB] > ystart[k]*HALK2DIC){
	  fprintf(stderr,"\nBox No. %d: ALK/DIC > %.2f",k,HALK2DIC);
	  fwarn("intifree(): Your ALK/DIC ratio is high. csys() may abort.");
      break;
   }	 
 } 
 /* check PO4 high-lat */
 if(ystart[2*NB+10] < 0.5){
      fprintf(stderr,"\nHighLat PO4 = %e mmol/m3",ystart[2*NB+10]);
      fwarn("intifree(): Your HighLat PO4 is low. derivs() may abort.");
 }
	 
 /*============= parms.out =============*/	 
 fpout = fopen("parms.out","w");
 if(fpout == NULL) 
	ferrx("initfree(): Can't open file 'parms.out'");
	 
 fprintf(fpout,"%d ffflag  \n",ffflag);
 fprintf(fpout,"%d tsnsflag\n",tsnsflag);
 fprintf(fpout,"%e sclim   \n",sclim);
 for(k=1;k<=3;k++)
 fprintf(fpout,"%f fdapi[%d]\n",fdapi[k],k);
 fprintf(fpout,"%e fcsml   \n",fcsml);
 fprintf(fpout,"%e epslvr  \n",epslvr);
 fprintf(fpout,"%d kmax    \n",kmax);
 fprintf(fpout,"%f dxsav   \n",dxsav);
 for(k=1;k<=NB;k++)
 fprintf(fpout,"%e vb[%d]  \n",vb[k],k);
 for(k=1;k<=NB;k++)
 fprintf(fpout,"%e ab[%d]  \n",ab[k],k);
 for(k=1;k<=NB;k++)
 fprintf(fpout,"%e hb[%d]  \n",hb[k],k);
 for(k=1;k<=NB;k++)
 fprintf(fpout,"%e tcb0[%d]\n",tcb0[k],k);
 for(k=1;k<=NB;k++)
 fprintf(fpout,"%e salv[%d]\n",salv[k],k);
 for(k=1;k<=NB;k++)
 fprintf(fpout,"%e prsv[%d]\n",prsv[k],k);
 fprintf(fpout,"%e t0      \n",t0);
 fprintf(fpout,"%e tfinal  \n",tfinal);
 fprintf(fpout,"%e cinp    \n",cinp);
 fprintf(fpout,"%e dccinp  \n",dccinp);
 fprintf(fpout,"%e rccinp  \n",rccinp);
 fprintf(fpout,"%e tcin0   \n",tcin0);
 fprintf(fpout,"%e tcspan  \n",tcspan);
 fprintf(fpout,"%e pcsi    \n",pcsi);
 fprintf(fpout,"%e thc0    \n",thc0);
 fprintf(fpout,"%e thc     \n",thc);
 fprintf(fpout,"%e thbra   \n",thbra);
 fprintf(fpout,"%e fepl    \n",fepl);
 fprintf(fpout,"%e eph     \n",eph);
 fprintf(fpout,"%e rrain   \n",rrain);
 fprintf(fpout,"%e fsh     \n",fsh);	 
 fprintf(fpout,"%e fvc0    \n",fvc0);	 
 fprintf(fpout,"%e finc0   \n",finc0);	 
 fprintf(fpout,"%e fkrg    \n",fkrg);	 
 fprintf(fpout,"%e cac     \n",cac);	 
 fprintf(fpout,"%e mgc     \n",mgc);
 fprintf(fpout,"%e s4c     \n",s4c);	 
 for(k=1;k<=KOC;k++)
 fprintf(fpout,"%e mxv[%d] \n",mxv[k],k);
 for(k=1;k<=NOC;k++)
 fprintf(fpout,"%e mhd[%d] \n",mhd[k],k);
#ifdef FSED
 fprintf(fpout,"%e ncsd     \n",ncsd);	 
 fprintf(fpout,"%e kssd     \n",kssd);	 
 for(k=1;k<=NSD;k++)
 fprintf(fpout,"%e dsv[%d]  \n",dsv[k],k);
 for(k=1;k<=NSD;k++)
 fprintf(fpout,"%e asva[%d] \n",asva[k],k);
 for(k=1;k<=NSD;k++)
 fprintf(fpout,"%e asvi[%d] \n",asvi[k],k);
 for(k=1;k<=NSD;k++)
 fprintf(fpout,"%e asvp[%d] \n",asvp[k],k);
#ifdef FTYS
 for(k=1;k<=NSD;k++)
 fprintf(fpout,"%e asvt[%d] \n",asvt[k],k);
#endif		 
 for(k=1;k<=NSD;k++)
 fprintf(fpout,"%d klid[%d] \n",klid[k],k);
 for(k=1;k<=3;k++)
 fprintf(fpout,"%d nlid[%d] \n",nlid[k],k);
#ifdef FTYS
 for(k=1;k<=NSD;k++)
 fprintf(fpout,"%d klidt[%d] \n",klidt[k],k);
#endif		 
 for(k=1;k<=NSD;k++)
 fprintf(fpout,"%e phiia[%d] \n",phiia[k],k);
 for(k=1;k<=NSD;k++)
 fprintf(fpout,"%e phiii[%d] \n",phiii[k],k);
 for(k=1;k<=NSD;k++)
 fprintf(fpout,"%e phiip[%d] \n",phiip[k],k);
#endif		 
 fclose(fpout);	

#ifdef FSED
 /*============= dsv.out =============*/	 
 fpout = fopen("dsv.out","w");
 if(fpout == NULL) 
	ferrx("initfree(): Can't open file 'dsv.out'");
 for(k=1;k<=NSD;k++)
    fprintf(fpout,"%e\n",dsv[k]);	 
 fclose(fpout);	

/*=============  zv.out =============*/	 
 fpout = fopen("zv.out","w");
 if(fpout == NULL) 
	ferrx("initfree(): Can't open file 'zv.out'");
 for(k=1;k<=nz;k++)
    fprintf(fpout,"%e\n",zv[k]);	 
 fclose(fpout);	
#endif		 
	 
	 
 /*============= ystart.out =============*/	 
 fpout = fopen("ystart.out","w");
 if(fpout == NULL) 
	ferrx("initfree(): Can't open file 'ystart.out'");
	 
 for(i=1;i<=NEQ;i++)
     fprintf(fpout,"%e ystart[%d] \n",ystart[i],i);
 fclose(fpout);	

 } /* END flag == 3 */
	

 /*============================================================*/
 /*============================================================*/
 /* free globals */	
 if(flag == 4){

 free_dvector(ystart,1,NEQ);
 free_dvector(vb,1,NB);
 free_dvector(ab,1,NB);
 free_dvector(hb,1,NB);
 free_ivector(kkv,1,NTS);	 
 free_ivector(kiv,1,NOC);	 
 free_dvector(mxv,1,KOC);
 free_dvector(mhd,1,NOC);
 free_dvector(gp,1,NB);
 free_dvector(hgssv,1,NB);
 free_dvector(tcb0,1,NB);
 free_dvector(tcb0c,1,NB);
 free_dvector(salv,1,NB);
 free_dvector(prsv,1,NB);
 free_dmatrix(spm,1,NB,1,2);
 free_dvector(fdapi,1,3);
 free_dvector(tmv,1,kmax);
 free_dmatrix(yy,1,NEQ,1,kmax);
	 
 if(ffflag == 1){	 
   free_dvector(tems,1,ltem);
   free_dvector(yems,1,ltem);	 
 }	 

#ifdef FSED
 free_dvector(dsv,1,NSD);
 free_dvector(zv,1,nz);
 free_dvector(asva,1,NSD);
 free_dvector(asvi,1,NSD);
 free_dvector(asvp,1,NSD);
 free_ivector(klid,1,NSD);
 free_ivector(nlid,1,3);
 free_dvector(fc0a,1,NSD);
 free_dvector(fc0i,1,NSD);
 free_dvector(fc0p,1,NSD);	 
 free_dvector(phiia,1,NSD);
 free_dvector(phiii,1,NSD);
 free_dvector(phiip,1,NSD);
 #ifdef FTYS
 free_dvector(asvt,1,NSD);
 free_ivector(klidt,1,NSD);
 free_dvector(fc0t,1,NSD);
 free_dvector(phiit,1,NSD);
 #endif	 
 #ifdef FSEDCC
 free_dvector(fcc0a,1,NSD);
 free_dvector(fcc0i,1,NSD);
 free_dvector(fcc0p,1,NSD);	 
  #ifdef FTYS
  free_dvector(fcc0t,1,NSD);
  #endif /*FTYS    */	
 #endif	/* FSEDCC */
#endif /* FSED   */

 printf("\n");
 printf("   ************************************************\n");	
 printf("   * LOSCAR V %s  Done.                        *\n",VLOSCAR);	
 printf("   ************************************************\n\n");	
 printf("  Nunc est bibendum.\n\n");
 timer=time(NULL);
 printf("%s",asctime(localtime(&timer)));
	 
 } /* END flag == 4 */ 

	
}
/*============================================================*/
/*==================== initfree() END ========================*/
/*============================================================*/
