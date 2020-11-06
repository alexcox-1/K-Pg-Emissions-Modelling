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
#include "defs.h"
#include "utils.h"

#include "common.h"

/*============================================================*/
/*==================== initstart() ===========================*/
/*============================================================*/
/* initialize start values                                    */
/*
   updates: 

   12/06/17 restart: recalculate sed 13C (global). see "Dec 2017"
   10/29/11 fdapi[1] applied to dicc
   09/16/11 include 13C (dicc)
   09/11/11 include dissolved O2 (dox) 
   04/19/11 added fdapi (change initial dic, alk, po4)
   04/15/11 temperature: added tcb as tracer (rm tcv). start priorities:
            1. control file 2. restart file 3. intern default
            NOTE: restart is read after cntrl => remember tcb0c
   04/13/11 tcb as tracer (temp [deg C] boxes)  
   04/10/11 removed ystart as arg to initstart(). ystart is
	        global and initstart() includes common.h
   04/07/11 rdflag: double -> int (bugfix)
   03/20/11 added load restart. NOTE: fcx and phiix are global
            and need to be recalculated during restart
   02/18/11 new file
 
 */

void initstart()
{
 int i,k,kfin,ktmp,lk=0,rdflag=1;
 double *onv,*dic0,*alk0,*po40,catm0,*d13c0,*rccb0,*dicc0,
	     d13catm0,ccatm0,*dtmp;
 char mssg[BUFSIZ];	
	
 FILE *fpld;

#ifdef FSED
 double ffphi,*mc0a,*mc0i,*mc0p;
 #ifdef FTYS	 
 double *mc0t;
 #endif
#endif
	
 /* init dicc                L    I    D  */
 double d13clid[NOC+1] = {0,2.35,0.5,0.67};
	
 /* default: XDBG (debug) */ 
 for(i=1;i<=NEQ;i++)
    ystart[i] = XDBG; /* XDBG, 0.0, (double)(i) */

 /*	
   Sequence of boxes (N=NB):
      1,2,3 = ATL,IND,PAC Surface
      4,5,6 = ATL,IND,PAC Interm
      7,8,9 = ATL,IND,PAC Deep
	     10 = High
  (11,12,13 = Surface,Interm,Deep TETH)		   
	 
   Sequence of equations:
 
   ***OCN tracers (N=NOCT)
   1 = dic
   2 = alk
   3 = po4
   4 = tcb
   5 = dox
   6 = dicc (total OCN equations = NB*NOCT)
   Catm
   C13atm   (NOATM = NOCT*NB+NCATM+NCCATM)
	   
   ***SED tracers (N=NSED per basin) 
   NOATM      +1... = fc ATL  	   
   NOATM+1*NSD+1... = fc IND  	   
   NOATM+2*NSD+1... = fc PAC  	   
  (NOATM+3*NSD+1... = fc TETH) 
  
   NOATM+(3+KTY)*NSD+1... = fcc ATL
   NOATM+(4+KTY)*NSD+1... = fcc IND
   NOATM+(5+KTY)*NSD+1... = fcc PAC
   NOATM+(6+KTY)*NSD+1... = fcc TETH

   where KTY = 0 modern, KTY = 1 paleo  
	   
   NEQ = NOCT*NB+NCATM+NCCATM+NOC*NSD+NOC*NSDCC
 */
	
 onv   = dvector(1,NB);
 dic0  = dvector(1,NB);
 alk0  = dvector(1,NB);
 po40  = dvector(1,NB);
 d13c0 = dvector(1,NB);
 rccb0 = dvector(1,NB);
 dicc0 = dvector(1,NB);

 /* 1-vector */ 	
 fsetv(onv,1.0,NB);

 /* init DIC */	
 vsclr(dic0,(2.30*1.e-3*RHO),onv,NB); /* 2.3 mmol/kg => mol/m3 */

 /* init ALK */	
 vsclr(alk0,(2.40*1.e-3*RHO),onv,NB); /* 2.4 mmol/kg => mol/m3 */
	
 /* init PO4 */	
#ifdef FTYS	
#define PO4PE (0.87*1.0133122)	
 vsclr(po40,(2.50*1.e-3*PO4PE),onv,  NB); /* mmol/kg (PO4 at t=0)  */
#else
 vsclr(po40,(2.50*1.e-3*0.87) ,onv,  NB); /* mmol/kg (PO4 at t=0)  */
#endif	
 vsclr(po40,      (1.   *RHO) ,po40 ,NB); /* mmol/m3 in derivs:    */
 /* x0.87 tunes to modern PO4 data */    /* convert to mol/m3 !!! */

if(NOCT >= 6){
 /* init dicc                L    I    D 
 double d13clid[NOC+1] = {0,2.35,0.5,0.67}; */
 for(k=1;k<=3;k++){
       d13c0[k]   = d13clid[1];	
       d13c0[k+3] = d13clid[2];	
       d13c0[k+6] = d13clid[3];	
 } 
 d13c0[10] = 1.63;	/* HighLat */
#ifdef FTYS
 for(k=1;k<=3;k++)
       d13c0[k+10]   = d13clid[k];	
#endif 	
 for(k=1;k<=NB;k++)
		rccb0[k] = (d13c0[k]/1e3+1.)*RST;
 vvelm(dicc0,rccb0,dic0,NB); /* mol/m3 */
}
	
 /* copy all to ystart */
 for(k=1;k<=NB;k++){
	 ystart[k]      = dic0[k];
	 ystart[k+1*NB] = alk0[k];
	 ystart[k+2*NB] = po40[k];
     kfin = k+2*NB;
 }
	
 
if(NOCT >= 4){
 /* copy tcb0 to ystart. tcb0 was set 
	as intern default, see initfree()
 */
 for(k=1;k<=NB;k++){
	 ystart[k+3*NB] = tcb0[k]/TSCAL; /* scaling !!! */
     kfin = k+3*NB;
 }	
} /* NOCT */


if(NOCT >= 5){
 /* dox */
 for(k=1;k<=NB;k++){
	 ystart[k+4*NB] = 0.2; /* mol/m3 */
     kfin = k+4*NB;
 }	
} /* NOCT */

if(NOCT >= 6){
 /* dicc */	
 for(k=1;k<=NB;k++){                 /* mol/m3 ->         */
	 ystart[k+5*NB] = dicc0[k]/CNTI; /* scaling /CNTI !!! */
     kfin = k+5*NB;
 }	
} /* NOCT */


if(NCATM == 1){
 catm0 = 280.*PPMTOMMSQ; 
 /* ppmv => (mol/m2) atm. CO2 inventory / m2      */
 /* 1 ppmv = 2.2 Pg C	                          */
 ystart[NOCT*NB+1] = catm0*CNTI; /* scaling *CNTI !!!  */
 kfin = NOCT*NB+1;
}

if(NCCATM == 1){
 /* atmospheric 13CO2 */
 d13catm0 = -6.45;
 ccatm0   = catm0*(d13catm0/1e3+1.)*RST;	
 ystart[NOCT*NB+2] = ccatm0;    /* mol/m2, not scaled! */
 kfin = NOCT*NB+2;
}

 free_dvector(onv,1,NB);
 free_dvector(dic0,1,NB);
 free_dvector(alk0,1,NB);
 free_dvector(po40,1,NB);
 free_dvector(d13c0,1,NB);
 free_dvector(rccb0,1,NB);
 free_dvector(dicc0,1,NB);


#ifdef FSED
 mc0a   = dvector(1,NSD);
 mc0i   = dvector(1,NSD);
 mc0p   = dvector(1,NSD);
#ifdef FTYS
 mc0t   = dvector(1,NSD);
#endif

 /* set initial calcite fraction */
 fsetv(fc0a,0.46,NSD);
 fsetv(fc0i,0.46,NSD);
 fsetv(fc0p,0.46,NSD);
#ifdef FTYS	 
 fsetv(fc0t,0.46,NSD);
#endif

 /* copy all to ystart             */
 /* NOATM = (NOCT*NB+NCATM+NCCATM) */
 for(k=1;k<=NSD;k++){
	 ystart[NOATM+k]       = fc0a[k];
	 ystart[NOATM+k+1*NSD] = fc0i[k];
	 ystart[NOATM+k+2*NSD] = fc0p[k];	 
     kfin = NOATM+k+2*NSD;
#ifdef FTYS	 
	 ystart[NOATM+k+3*NSD] = fc0t[k];	 
     kfin = NOATM+k+3*NSD;
#endif
 } 

 /* initial porosities and CaCO3 mass */	 
 ffphi = (phi1-phi0)/(1.-phi1);
 for(k=1;k<=NSD;k++){
     phiia[k] = (phi0+ffphi*fc0a[k])/(1.+ffphi*fc0a[k]);
     phiii[k] = (phi0+ffphi*fc0i[k])/(1.+ffphi*fc0i[k]);
     phiip[k] = (phi0+ffphi*fc0p[k])/(1.+ffphi*fc0p[k]);
      mc0a[k] = fc0a[k]*RHOS*(1.-phiia[k]); /* kg CaCO3/m3 */
      mc0i[k] = fc0i[k]*RHOS*(1.-phiii[k]); /* kg CaCO3/m3 */
      mc0p[k] = fc0p[k]*RHOS*(1.-phiip[k]); /* kg CaCO3/m3 */
	 #ifdef FTYS
     phiit[k] = (phi0+ffphi*fc0t[k])/(1.+ffphi*fc0t[k]);
      mc0t[k] = fc0t[k]*RHOS*(1.-phiit[k]); /* kg CaCO3/m3 */
     #endif
 }

 #ifdef FSEDCC
 vsclr(fcc0a,rincc,fc0a,NSD);
 vsclr(fcc0i,rincc,fc0i,NSD); 
 vsclr(fcc0p,rincc,fc0p,NSD);      
#ifdef FTYS
 vsclr(fcc0t,rincc,fc0t,NSD);      
#endif

 /* copy all to ystart             */
 /* NOATM = (NOCT*NB+NCATM+NCCATM) */
 for(k=1;k<=NSD;k++){
	 ystart[NOATM+k+(3+KTY)*NSD] = fcc0a[k]/CNTI; /* scaling /CNTI !!! */
	 ystart[NOATM+k+(4+KTY)*NSD] = fcc0i[k]/CNTI; /* scaling /CNTI !!! */
	 ystart[NOATM+k+(5+KTY)*NSD] = fcc0p[k]/CNTI; /* scaling /CNTI !!! */
     kfin = NOATM+k+(5+KTY)*NSD;
#ifdef FTYS	 
	 ystart[NOATM+k+(6+KTY)*NSD] = fcc0t[k]/CNTI; /* scaling /CNTI !!! */
     kfin = NOATM+k+(6+KTY)*NSD;
#endif
 } 
 #endif

 free_dvector(mc0a,1,NSD);
 free_dvector(mc0i,1,NSD);
 free_dvector(mc0p,1,NSD);
#ifdef FTYS
 free_dvector(mc0t,1,NSD);
#endif

#endif

 printf("\n@ No. default start values : %3d\n",kfin);
 if(kfin != NEQ)
   ferrx("initstart(): # start values don't match # equations.");


 /*============================================================*/
 /*============= alternatively: load restart ==================*/
 if(ldrestart == 1){

 /* open restart file */
 fpld = fopen(fpldstr,"r");
 if(fpld == NULL){ 
    sprintf(mssg,"initstart(): Can't open restart file '%s'",fpldstr);
	ferrx(mssg);
 } else {
 	printf("\n@ This is a restart: loading '%s'\n",fpldstr);
 }

 /* read restart values */
 dtmp = dvector(1,NEQ);
 k = 1;
 while(rdflag != EOF){
    rdflag = fscanf(fpld,"%d %le",&ktmp,&dtmp[k]);
    if(rdflag == 2){
	  lk = k;
      /*printf("%d %d %e\n",lk,ktmp,dtmp[k]);*/
	  if(lk > NEQ)
         ferrx("initstart(): Too many restart values (>NEQ)");
      k += 1;
	} else {
	break;
	}
 }
	 
 /* close restart file */
 fclose(fpld);

 if(lk < NEQ)
    ferrx("initstart(): Insufficient restart values (<NEQ)");
	 
 /* copy numbers read from file to ystart */	 
 for(k=1;k<=NEQ;k++){
	 ystart[k] = dtmp[k];
     /*printf("%d %e\n",k,ystart[k]);*/
 }

 /*================= handle temperature =============*/
 /*	 
  start priorities for temperature are:
  1. control file 2. restart file 3. intern default 
  NOTE: restart is read after cntrl => remember tcb0c
 */
	 
if(NOCT >= 4){	 
 /* set initial temperature from restart */	 
 for(k=1;k<=NB;k++)
	 tcb0[k] = ystart[k+3*NB];	/* deg C */
 /* BUT: control file has priority!
    if tcb0c was set in control file (readparms()), then
    - set initial temperature from control
	- copy tcb0c to ystart 
 */
 if(cntrfflag == 1){
   for(k=1;k<=NB;k++){
	   tcb0[k]        = tcb0c[k];
	   ystart[k+3*NB] = tcb0c[k];
   }
 }
 /* scale ystart temperature (input is in deg C) */
 for(k=1;k<=NB;k++)
	 ystart[k+3*NB] /= TSCAL;	 
} /* NOCT */


#ifdef FSED
 /* fc0_i and phii_i are global and need to be  */
 /* recalculated during restart.                */

	 /* NOATM = (NOCT*NB+NCATM+NCCATM) */
 for(k=1;k<=NSD;k++){
	 fc0a[k] = ystart[NOATM+k]      ;
	 fc0i[k] = ystart[NOATM+k+1*NSD];
	 fc0p[k] = ystart[NOATM+k+2*NSD];
     #ifdef FTYS /* TETHYS */	 
	 fc0t[k] = ystart[NOATM+k+3*NSD];
     #endif
 } 

 /* initial porosities */	 
 for(k=1;k<=NSD;k++){
     phiia[k] = (phi0+ffphi*fc0a[k])/(1.+ffphi*fc0a[k]);
     phiii[k] = (phi0+ffphi*fc0i[k])/(1.+ffphi*fc0i[k]);
     phiip[k] = (phi0+ffphi*fc0p[k])/(1.+ffphi*fc0p[k]);
     #ifdef FTYS /* TETHYS */	 
     phiit[k] = (phi0+ffphi*fc0t[k])/(1.+ffphi*fc0t[k]);
     #endif
 }
#endif	 
#ifdef FSEDCC /* 13C (Dec 2017) */
 for(k=1;k<=NSD;k++){
	 fcc0a[k] = ystart[NOATM+k+(3+KTY)*NSD]*CNTI; 
	 fcc0i[k] = ystart[NOATM+k+(4+KTY)*NSD]*CNTI;
	 fcc0p[k] = ystart[NOATM+k+(5+KTY)*NSD]*CNTI;
     #ifdef FTYS /* TETHYS */	 
	 fcc0t[k] = ystart[NOATM+k+(6+KTY)*NSD]*CNTI;
     #endif
 }
#endif
	 
 free_dvector(dtmp,1,NEQ);
 
 } /* ldrestart */


 /* apply fdapi (change initial dic, alk, po4)        */ 
 /* only control file modifies fdapi. default: 00.00% */

 for(k=1;k<=NB;k++){
	 ystart[k]      *= (1.+fdapi[1]*CNTI);
	 ystart[k+1*NB] *= (1.+fdapi[2]*CNTI);
	 ystart[k+2*NB] *= (1.+fdapi[3]*CNTI);
 }

if(NOCT >= 6){
 /* dicc */	
 for(k=1;k<=NB;k++){    
	 ystart[k+5*NB] *= (1.+fdapi[1]*CNTI);
 }	
} /* NOCT */

}
/*============================================================*/
/*==================== initstart() END =======================*/
/*============================================================*/
