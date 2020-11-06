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
#include <math.h>
#include <stdlib.h>
#include "defs.h"
#include "utils.h"

#include "common.h"

/* compiler flag TIMESTEP, if undefined (ifndef), then default to 5000 */

/*============================================================*/
/*==================== derivs() ==============================*/
/*============================================================*/
/* User supplied routine: right-hand side of DEQs             */
/*
 
 updates: 

  12/02/11 1./bbf/bbf -> fcsml
  10/25/11 co3satm[aipt] all basins. nsd -> sdn 
  10/23/11 helper indices for boxes: m1,m2,m3. T-sensitivity
  10/22/11 Tethys sediments
  10/21/11 including Tethys. new thmfun() 
  10/18/11 eplv  = dvector(1,3=>NLS) etc.
  10/15/11 include shelf/deep rain
  10/08/11 add Ca,Mg to fcsys() 
  09/23/11 include sediment 13C	   
  09/16/11 include 13C (dicc).
  09/15/11 added fflush(stdout), fflush(stderr).
  09/12/11 MM kinetics dissolved oxygen
           High-Lat PO4 < 0: error -> warning.
  09/11/11 include dissolved O2 (dox). o2 -> o2sat. 
           added vask (gas exchange coeff. O2).
  05/06/11 added comments on sediment equations
  04/20/11 added MM kinetics HighLat PO4
  04/13/11 added tcb as tracer (temp [deg C] boxes).  
  04/10/11 separate: fvc0 and fvc(t)
           renamed vars: ga>gtha, gi>gthi
  04/09/11 TempSens: surface relax: 10->20 y
  04/03/11 replaced femiss() by finterp()
  03/14/11 1st sediment test OK. fcsml as parameter
  03/06/11 set epsj=epslvr
  02/14/11 new file	  

 NOTE:
	  "It is absolutely crucial to scale your variables
	  properly when integrating stiff problems with automatic 
	  stepsize adjustment." (NR)

	  Thus variables have been scaled (if necessary) to be of 
      order 1 (including temperature) before passing to the solver.
      Here (in derivs()) variables are scaled to their usual
      values before calculating dy/dt. dy/dt is scaled back to 
      order 1 at the end. See also: defs.h, initstart.c, solver.c
	  
      solver uses: yscal[i]=DMAX(CSCAL,fabs(y[i])); where CSCAL
      is set in defs.h
	  
 */
void derivs(double t,double *y, double *yp)
{
 int k,kk;	
 double *dic,*alk,*po4,*tcb,*dox,*dicc,*dicp,*alkp,*po4p,
	    *tcbp,*doxp,*diccp;
 double *co2,*pco2,*co3,*hpls,*ph,*kh,*o2sat;
 double *alpdb,*alpdg,*alpcb,alpu,alpcorg;
 double *eplv,*ealv,*enlv,*pplv,*eclv,*exlv,eah,pph,enh,ech,exh,
	    oi,gtha,gthi,sdn,pco2a,catmp,*kasv0,*kasv,fvc,fsi,finc,
        fvccc,fsicc,fincc,fkrgcc,*fkrgccb,d13cems,remscc,
	    *vask0,*vask,*fmmo,fems=0,femscc,tmp;
    
 /*DP: flux of H2SO4 from volcanic degassing*/
 double fso4;
    
 double *eplvcc,*ealvcc,*eclvcc,*exlvcc,ephcc,eahcc,echcc,exhcc,
	    *rccb,pcco2a,ccatmp;
 char mssg[BUFSIZ];	
#ifdef FSED
 int i,j,fc0flag=0,fcc0flag=0,fc1flag=0,*kda,*kdi,*kdp,*lea,*lei,*lep,
	 *jja,*jji,*jjp,ns1,ns2,jsh=2;
 double *fcva,*fcvi,*fcvp,*fcvap,*fcvip,*fcvpp,
	    *frrfa,*frrfi,*frrfp,eaa,eai,eap,fpra,fpri,fprp,
	    *fprva,*fprvi,*fprvp,*tsedv,*ssedv,kspcsed,**co3satm,
        *omva,*omvi,*omvp,ffphi,*phia,*phii,*phip,
        *rscva,*rscvi,*rscvp,*rsrva,*rsrvi,*rsrvp,
	    *rsva,*rsvi,*rsvp,*dka,*dki,*dkp,*dissa,*dissi,*dissp,
	    *rdva,*rdvi,*rdvp,*wva,*wvi,*wvp,*wcva,*wcvi,*wcvp,
	    *fba,*fbi,*fbp,*gsda,*gsdi,*gsdp,**dissm,bbf,rcak,
        *fdr,*fshva,*fshvi,*fshvp;
  #ifdef FTYS
  int *kdt,*lett,*jjt;
  double *fcvt,*fcvtp,*frrft,eatt,fprt,*fprvt,*omvt,*phit,
	     *rscvt,*rsrvt,*rsvt,*dkt,*disst,*rdvt,*wvt,*wcvt,
	     *fbt,*gsdt,fsht,*fshvt;
  #endif
  #ifdef FSEDCC
  double *fccva,*fccvi,*fccvp,*fccvap,*fccvip,*fccvpp,
	     eaacc,eaicc,eapcc,*fpracc,*fpricc,*fprpcc,
	     *rccsa,*rccsi,*rccsp,*rscvacc,*rscvicc,*rscvpcc,
         *rdvacc,*rdvicc,*rdvpcc,*wcvacc,*wcvicc,*wcvpcc,
         **dissmcc;
   #ifdef FTYS
   double *fccvt,*fccvtp,eatcc,*fprtcc,*rccst,*rscvtcc,
	      *rdvtcc,*wcvtcc;
   #endif /* FTYS   */
 #endif	 /* FSEDCC */
#endif	/* FSED   */

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

 /* Tethys:
    Set helper indices for boxes. To maintain modern order,
	Tethys is added at the end (Tethys LID = 11,12,13).
	The following array indices (wich do not work for 
    k=4=tethys)	are hence replaced by:
    [k]      =>  [m1[k]]
    [k+3]    =>  [m2[k]+3]
    [k+6]    =>  [m3[k]+6]
	with m1,m2,m3 set here.
 */ 
 int m1[4+1]={0,1,2,3,11};
 int m2[4+1]={0,1,2,3, 9};
 int m3[4+1]={0,1,2,3, 7};

 if(t < 0.0) ferrwrt("derivs(): time is negative");
	
 /* default: all derivatives to zero */
 for(k=1;k<=NEQ;k++)
    yp[k]    =  0.0; /* 0.0, XDBG, -y[k] */

 if(ffflag == 1){
      fems = finterp(t,tems,yems,ltem,0);
	 fems *= 1.e15/12.; /* Pg C/y => mol C/y   */ 
   d13cems = -5.0;     /* d13C emissions -25. */
    remscc = (d13cems/1.e3+1.)*RST;
    femscc = remscc*fems;
 }

 /* rename variable: y => dic,alk,po4,tcb ... 
	avoids index counting disaster */

 dic   = dvector(1,NB);
 alk   = dvector(1,NB);
 po4   = dvector(1,NB);
 tcb   = dvector(1,NB);
 dox   = dvector(1,NB);
 dicc  = dvector(1,NB);
 dicp  = dvector(1,NB);
 alkp  = dvector(1,NB);
 po4p  = dvector(1,NB);
 tcbp  = dvector(1,NB);
 doxp  = dvector(1,NB);
 diccp = dvector(1,NB);
 rccb  = dvector(1,NB);

 for(k=1;k<=NB;k++){
    dic[k] = y[k];
    alk[k] = y[k+1*NB];
    po4[k] = y[k+2*NB]*MMTOM; /* PO4 convert to mol/m3 */
    /* dydt default: all to zero */
   dicp[k] = 0.0;
   alkp[k] = 0.0;
   po4p[k] = 0.0;
 }

 /* warning if high-lat PO4 < 0 */
 if(po4[10] < 0.0){
   fprintf(stderr,"\nHighLat PO4 = %e mol/m3",po4[10]);
   fwarn("derivs(): HighLat PO4 is negative. Increase initial PO4?");
 }

 /* error if low-lat PO4 < 0 */
 for(k=1;k<=NLS;k++){	
  if(po4[k] < 0.0){
    fprintf(stderr,"\nLowLat PO4[%d] = %e mol/m3",k,po4[k]);
    ferrwrt("derivs(): LowLat PO4 is negative. Reduce FBIOL?");
  }
 }	

	
if(NOCT >= 4){
 for(k=1;k<=NB;k++){
    tcb[k] = y[k+3*NB]*TSCAL; /* Temp convert to degC  */
   tcbp[k] = 0.0;
 }	 
} else {
 for(k=1;k<=NB;k++)
    tcb[k] = tcb0[k];         /* Temp = Temp0 = const. */
 }

if(NOCT >= 5){
 /* dox */
 for(k=1;k<=NB;k++){
    dox[k] = y[k+4*NB];       /* mol/m3 */
   doxp[k] = 0.0;
 }	 
 /* warning if oxygen < 0 */
 for(k=1;k<=NB;k++){	
  if(dox[k] < 0.0){
    fprintf(stderr,"\nDissolved Oxygen: dox[%d] = %e mol/m3",k,dox[k]);
    fwarn("derivs(): Diss. O2 is negative (<anoxic). Reduce CBIOH?");
  }
 }
}

if(NOCT >= 6){
 /* dicc */
 for(k=1;k<=NB;k++){          /* scaling *CNTI !!! */
    dicc[k] = y[k+5*NB]*CNTI; /* -> mol/m3         */
   diccp[k] = 0.0;
    rccb[k] = dicc[k]/dic[k]; /* 13R-DIC of boxes  */
 }
}


if(NCATM == 1){	
 /* atmospheric CO2 */	
 pco2a = y[NOCT*NB+1]/PPMTOMMSQ/CNTI; /* scaling /CNTI */
 catmp = 0.0;             /* dydt default: all to zero */
}	

if(NCCATM == 1){	
 /* atmospheric 13CO2 */	
 pcco2a = y[NOCT*NB+2]/PPMTOMMSQ;    /* not scaled !!! */
 ccatmp = 0.0;            /* dydt default: all to zero */
}	


#ifdef FSED
 fcva  = dvector(1,NSD);
 fcvi  = dvector(1,NSD);
 fcvp  = dvector(1,NSD);
 fcvap = dvector(1,NSD);
 fcvip = dvector(1,NSD);
 fcvpp = dvector(1,NSD);
 #ifdef FTYS
 fcvt  = dvector(1,NSD);
 fcvtp = dvector(1,NSD); 
 #endif

 /* rename */
 for(k=1;k<=NSD;k++){
    fcva[k] = y[NOATM+k];
    fcvi[k] = y[NOATM+k+1*NSD];
    fcvp[k] = y[NOATM+k+2*NSD];
    #ifdef FTYS
    fcvt[k] = y[NOATM+k+3*NSD];
    #endif	 
    /* dydt default: all to zero */
   fcvap[k] = 0.0;
   fcvip[k] = 0.0;
   fcvpp[k] = 0.0;
   #ifdef FTYS
   fcvtp[k] = 0.0;
   #endif	 
 } 

#ifdef FSEDCC
 fccva  = dvector(1,NSD);
 fccvi  = dvector(1,NSD);
 fccvp  = dvector(1,NSD);
 fccvap = dvector(1,NSD);
 fccvip = dvector(1,NSD);
 fccvpp = dvector(1,NSD);
 #ifdef FTYS
 fccvt  = dvector(1,NSD);
 fccvtp = dvector(1,NSD);
 #endif

 /* rename */
 for(k=1;k<=NSD;k++){
    fccva[k] = y[NOATM+k+(3+KTY)*NSD]*CNTI; /* scaling *CNTI */
    fccvi[k] = y[NOATM+k+(4+KTY)*NSD]*CNTI; /* scaling *CNTI */
    fccvp[k] = y[NOATM+k+(5+KTY)*NSD]*CNTI; /* scaling *CNTI */
    #ifdef FTYS
    fccvt[k] = y[NOATM+k+(6+KTY)*NSD]*CNTI; /* scaling *CNTI */
    #endif	 
    /* dydt default: all to zero */
    fccvap[k] = 0.0;
    fccvip[k] = 0.0;
    fccvpp[k] = 0.0;
    #ifdef FTYS
    fccvtp[k] = 0.0;
    #endif	 
 } 
#endif /* FSEDCC */

 /* error if any fc < 0 */
 for(k=1;k<=NSD;k++){
    if(fcva[k] < 0.0 || fcvi[k] < 0.0 || fcvp[k] < 0.0){
		fc0flag = 1;
		break;
	}
#ifdef FSEDCC
    if(fccva[k] < 0.0 || fccvi[k] < 0.0 || fccvp[k] < 0.0){
		fcc0flag = 1;
		break;
	}
#endif
 }
 if(fc0flag == 1){
	printf("\n\n@ at t=%e\n",t);
	ferrwrt("derivs(): fc < 0.0. Reduce EPSLV? Raise FCSML?");
 }
 if(fcc0flag == 1){
	printf("\n\n@ at t=%e\n",t);
	ferrwrt("derivs(): fcc < 0.0. Reduce EPSLV? Raise FCSML?");
 }
	 
 /* issue warning if any fc > 1.0 */
 for(k=1;k<=NSD;k++){
    if(fcva[k] > 1.0 || fcvi[k] > 1.0 || fcvp[k] > 1.0)
		 fc1flag = 1;
 }
 if(fc1flag == 1){
    sprintf(mssg,"derivs(): fc > 1 at t=%.2e. Check final fc values.",t);
    fwarn(mssg);	 
 }

 /* #sed boxes low and low+intm */
 ns1 = nlid[1];
 ns2 = nlid[1]+nlid[2];
#endif /* FSED */

 /* init locals */
 eplv  = dvector(1,NLS);
 ealv  = dvector(1,NLS);
 enlv  = dvector(1,NLS);
 pplv  = dvector(1,NLS);
 eclv  = dvector(1,NLS);
 exlv  = dvector(1,NLS);
 co2   = dvector(1,NB);
 pco2  = dvector(1,NB); /* ocean CO2 */
 co3   = dvector(1,NB);
 hpls  = dvector(1,NB);
 ph    = dvector(1,NB);
 kh    = dvector(1,NB);
 o2sat = dvector(1,NB);
 kasv0 = dvector(1,NB);
 kasv  = dvector(1,NB);
 vask0 = dvector(1,NB);
 vask  = dvector(1,NB);
 fmmo  = dvector(1,NB);
 alpdb = dvector(1,NB);
 alpdg = dvector(1,NB);
 alpcb = dvector(1,NB);
 eplvcc = dvector(1,NLS);
 ealvcc = dvector(1,NLS);
 eclvcc = dvector(1,NLS);
 exlvcc = dvector(1,NLS);
 fkrgccb = dvector(1,NB);
#ifdef FSED
  /* init locals */
 fprva = dvector(1,NSD);
 fprvi = dvector(1,NSD);
 fprvp = dvector(1,NSD);
 frrfa = dvector(1,NSD);
 frrfi = dvector(1,NSD);
 frrfp = dvector(1,NSD);

 tsedv = dvector(1,NSD);
 ssedv = dvector(1,NSD);

 co3satm = dmatrix(1,NOC,1,NSD);

 omva  = dvector(1,NSD);
 omvi  = dvector(1,NSD);
 omvp  = dvector(1,NSD);

 kda   = ivector(1,NSD);
 kdi   = ivector(1,NSD);
 kdp   = ivector(1,NSD);

 phia  = dvector(1,NSD);
 phii  = dvector(1,NSD);
 phip  = dvector(1,NSD);

 rscva = dvector(1,NSD);
 rscvi = dvector(1,NSD);
 rscvp = dvector(1,NSD);

 rsrva = dvector(1,NSD);
 rsrvi = dvector(1,NSD);
 rsrvp = dvector(1,NSD);

 rsva  = dvector(1,NSD);
 rsvi  = dvector(1,NSD);
 rsvp  = dvector(1,NSD);

 dka   = dvector(1,NSD);
 dki   = dvector(1,NSD);
 dkp   = dvector(1,NSD);

 dissa = dvector(1,NSD);
 dissi = dvector(1,NSD);
 dissp = dvector(1,NSD);

 rdva  = dvector(1,NSD);
 rdvi  = dvector(1,NSD);
 rdvp  = dvector(1,NSD);

 wva   = dvector(1,NSD);
 wvi   = dvector(1,NSD);
 wvp   = dvector(1,NSD);

 lea   = ivector(1,NSD);
 lei   = ivector(1,NSD);
 lep   = ivector(1,NSD);

 wcva  = dvector(1,NSD);
 wcvi  = dvector(1,NSD);
 wcvp  = dvector(1,NSD);

 fba   = dvector(1,NSD);
 fbi   = dvector(1,NSD);
 fbp   = dvector(1,NSD);

 gsda  = dvector(1,NSD);
 gsdi  = dvector(1,NSD);
 gsdp  = dvector(1,NSD);

 jja   = ivector(1,NSD);
 jji   = ivector(1,NSD);
 jjp   = ivector(1,NSD);

 fdr   = dvector(1,NOC);
 fshva = dvector(1,NSD);
 fshvi = dvector(1,NSD);
 fshvp = dvector(1,NSD);
 #ifdef FTYS
  fprvt = dvector(1,NSD);
  frrft = dvector(1,NSD);
  omvt  = dvector(1,NSD);
  kdt   = ivector(1,NSD);
  phit  = dvector(1,NSD);
  rscvt = dvector(1,NSD);
  rsrvt = dvector(1,NSD);
  rsvt  = dvector(1,NSD);
  dkt   = dvector(1,NSD);
  disst = dvector(1,NSD);
  rdvt  = dvector(1,NSD);
  wvt   = dvector(1,NSD);
  lett  = ivector(1,NSD);
  wcvt  = dvector(1,NSD);
  fbt   = dvector(1,NSD);
  gsdt  = dvector(1,NSD);
  jjt   = ivector(1,NSD);
  fshvt = dvector(1,NSD);
 #endif

 dissm = dmatrix(1,NOC,1,NSD);
 #ifdef FSEDCC
 rccsa   = dvector(1,NSD);
 rccsi   = dvector(1,NSD);
 rccsp   = dvector(1,NSD);

 fpracc  = dvector(1,NSD);
 fpricc  = dvector(1,NSD);
 fprpcc  = dvector(1,NSD);

 rscvacc = dvector(1,NSD);
 rscvicc = dvector(1,NSD);
 rscvpcc = dvector(1,NSD);

 rdvacc  = dvector(1,NSD);
 rdvicc  = dvector(1,NSD);
 rdvpcc  = dvector(1,NSD);

 wcvacc  = dvector(1,NSD);
 wcvicc  = dvector(1,NSD);
 wcvpcc  = dvector(1,NSD);
 #ifdef FTYS
  rccst   = dvector(1,NSD);
  fprtcc  = dvector(1,NSD);
  rscvtcc = dvector(1,NSD);
  rdvtcc  = dvector(1,NSD);
  wcvtcc  = dvector(1,NSD);
 #endif /* FTYS */
 dissmcc = dmatrix(1,NOC,1,NSD);
 #endif /* FSEDCC */
#endif /* FSED */

/* #include "test/x1.c" * @REZ */

 /* CO2 system and O2 of boxes (surface: 1=LA, 2=LI, 3=LP, 10=H)
   requires mol/kg                                           */

 for(k=1;k<=NB;k++){
   fcsys(&co2[k],&pco2[k],&co3[k],&hpls[k],&ph[k],&kh[k],&o2sat[k],
         dic[k]/RHO,alk[k]/RHO,hgssv[k],tcb[k],salv[k],prsv[k],
         cac,mgc,s4c);
   hgssv[k] = hpls[k];
 }

 /*********
 printf("%e  CO2\n",co2[1]);	 
 printf("%e pCO2\n",pco2[1]);	 
 printf("%e  CO3\n",co3[1]);	 
 printf("%e   pH\n",ph[1]);	 
 printf("%e   kh\n",kh[1]);	 
 printf("%e   O2\n",o2sat[1]);	 
 exit(0);
 **********/

 /*===== air-sea CO2 exchange coeff. =======*/
 fsetv(kasv0,0.0,NB);        /* set all to zero                 */
 for(k=1;k<=NTS;k++)
   kasv0[kkv[k]] = 0.06;     /* (mol/uatm/y/m2) CO2 Broecker    */
 vvelm(kasv,kasv0,ab,NB);    /* (mol/uatm/y)                    */
 	
 /*===== air-sea  O2 exchange coeff. =======*/
 fsetv(vask0,0.0,NB);       /* set all to zero                  */
 for(k=1;k<=NTS;k++)        /* Toggweiler 1999                  */
   vask0[kkv[k]] = 3.*365.; /* (m/day) -> (m/y) piston velocity */
 vvelm(vask,vask0,ab,NB);   /* (m3/y)                           */

 /*===== air-sea  13CO2 alphas =============*/
 falpha(alpdb,alpdg,alpcb,&alpu,tcb,kkv,NTS);

 /*=========== Biological Pump ==================*/
 /* Low Lat Export Corg                          */
   
    
    double ExpReduction = 1.;
    
    if(t>390.e3 && t<2160.e3)
        ExpReduction=0.5+((t-390.e3)*0.5/1770.e3);
    else ExpReduction = 1.;
    
    
    eplv[1]=ExpReduction*8.0956926743874953125e13; /*Atlantic*/
    eplv[2]=ExpReduction*8.2702707475210890625e13; /*Indian*/
    eplv[3]=ExpReduction*1.39160302612395453125e14; /*Pacific*/
    eplv[4]=ExpReduction*6.4594032486555375e13; /*Tethys*/
    
    if(t > 390.e3 && t < 2160.e3) {
        rrain = 5.8-((t-390.e3)*(5.8-6.7)/1770.e3);
    }
    else {
        rrain = 6.7;
    }
    
    
    /* Low Lat Export Other */
    vsclr(ealv,2./rrain,eplv,NLS);  /* ALK mol/y */
    vsclr(pplv,REDPC   ,eplv,NLS);  /* PO4 mol/y */
    vsclr(enlv,REDNC   ,eplv,NLS);  /* NO3 mol/y */
    
    /* total carbon: Corg+CaCO3 */
    for(k=1;k<=NLS;k++)
        eclv[k] = eplv[k]+0.5*ealv[k];
    
    
    /* High Lat Export */
    eah =       0.0;    /* ALK = 2.*eph/rrain ? */
    pph = eph*REDPC;    /* PO4 */
    enh = eph*REDNC;    /* NO3 */
    /* total carbon: Corg+CaCO3 */
    ech = eph+0.5*eah;
    
    if(NOCT >= 6){
        /*=========== Biological Pump 13C ==============*/
        alpcorg = epscorg/1.e3+1.;
        for(k=1;k<=3;k++){
            eplvcc[k] = alpcorg*rccb[k]*eplv[k];
            ealvcc[k] =         rccb[k]*ealv[k];
            eclvcc[k] = eplvcc[k]+0.5*ealvcc[k];
        }
        ephcc = alpcorg*rccb[10]*eph;
        eahcc =         rccb[10]*eah;
        echcc = ephcc+0.5*eahcc;
#ifdef FTYS
        k = 4;
        eplvcc[k] = alpcorg*rccb[11]*eplv[k];
        ealvcc[k] =         rccb[11]*ealv[k];
        eclvcc[k] =  eplvcc[k]+0.5*ealvcc[k];
#endif
    }
    
#ifdef MMPO4H
    /* MM kinetics if HighLat PO4 < PMMK mol/m3 */
    if(po4[10] < PMMK){
        pph *= PMMV*po4[10]/(PMM0+po4[10]);
        if(po4[10] < 0.0)
            pph = 0.0;
    }
#endif
    
    if(NOCT >= 5){
        /* MM kinetics dissolved oxygen */
        for(k=1;k<=NB;k++){
            fmmo[k] = dox[k]/(dox[k]+KMMOX);
            if(dox[k] < 0.0)
                fmmo[k] = 0.0;
        }
    }
    
    if (t>390.e3 && t<2160.e3)
        frei=0.95-((t-390.e3)/1770.e3*(0.95-0.78));
    else
        frei=0.78;
    
    
 /* fraction EPL, remineralized in I boxes */
 oi = 1.-frei;

#ifdef FSED
 /* CaCO3 export AIP (Atl, Ind, Pac) mol/y */
 eaa   = ealv[1]+eah*gp[7]; /* +EAH (no H seds, see below) */
 eai   = ealv[2]+eah*gp[8];
 eap   = ealv[3]+eah*gp[9];

 /* CaCO3 rain to sediments AIP (not equal to export) */
 /* fpr = export - water column diss                  */
 fpra  = (1.-nuwd)*eaa*.5/ab[1]; /* mol C/m2/y */
 fpri  = (1.-nuwd)*eai*.5/ab[2]; /* mol C/m2/y */
 fprp  = (1.-nuwd)*eap*.5/ab[3]; /* mol C/m2/y */

 /* CaCO3 rain to sediments AIP, vector        */
 fsetv(fprva,fpra,NSD);          /* mol C/m2/y */
 fsetv(fprvi,fpri,NSD);          /* mol C/m2/y */
 fsetv(fprvp,fprp,NSD);          /* mol C/m2/y */

 /* clay/remainder rain to sediments AIP */
 fsetv(frrfa,frrf*1.0,NSD);
 fsetv(frrfi,frrf*1.0,NSD);
 fsetv(frrfp,frrf*0.5,NSD); /* remote Pac: less terr. input */

 #ifdef FTYS /* TETHYS */
  eatt  = ealv[4];
  fprt  = (1.-nuwd)*eatt*.5/ab[11];
  fsetv(fprvt,fprt,NSD);     
  fsetv(frrft,frrf*6.0,NSD); /* x6.0 !!! @REZ */
 #endif

 #ifdef FSEDCC
 /* Sediment 13C */
 vdelm(rccsa,fccva,fcva,NSD);    /* 13Rsed */
 vdelm(rccsi,fccvi,fcvi,NSD);    /* 13Rsed */
 vdelm(rccsp,fccvp,fcvp,NSD);    /* 13Rsed */
 eaacc  = ealvcc[1]+eahcc*gp[7]; /* +EAH (no H seds, see below) */
 eaicc  = ealvcc[2]+eahcc*gp[8]; 
 eapcc  = ealvcc[3]+eahcc*gp[9];  
 fsetv(fpracc,(1.-nuwd)*eaacc*.5/ab[1],NSD); /* mol C/m2/y */
 fsetv(fpricc,(1.-nuwd)*eaicc*.5/ab[2],NSD); /* mol C/m2/y */
 fsetv(fprpcc,(1.-nuwd)*eapcc*.5/ab[3],NSD); /* mol C/m2/y */
  #ifdef FTYS /* TETHYS */
   vdelm(rccst,fccvt,fcvt,NSD);
   eatcc  = ealvcc[4]; 
   fsetv(fprtcc,(1.-nuwd)*eatcc*.5/ab[11],NSD);
  #endif /* FTYS   */
 #endif /* FSEDCC */

 /* shelf/deep rain partitioning   */
 /* calculate deep factor from fsh */
 fdr[1] = (1.-fsh*vsumlh(asva,NSD,    1,jsh))/
                  vsumlh(asva,NSD,jsh+1,NSD);
 fdr[2] = (1.-fsh*vsumlh(asvi,NSD,    1,jsh))/
                  vsumlh(asvi,NSD,jsh+1,NSD);
 fdr[3] = (1.-fsh*vsumlh(asvp,NSD,    1,jsh))/
                  vsumlh(asvp,NSD,jsh+1,NSD);
 #ifdef FTYS /* TETHYS */
  fsht = pow(fsh,nsht);
  fdr[4] = (1.-fsht*vsumlh(asvt,NSD,    1,jsh))/
                    vsumlh(asvt,NSD,jsh+1,NSD);
 #endif
 /* set shelf and deep factor */
 for(j=1;j<=jsh;j++){
	 fshva[j] = fsh;
	 fshvi[j] = fsh;
	 fshvp[j] = fsh;
     #ifdef FTYS /* TETHYS */
	 fshvt[j] = fsht;
     #endif	 
 }
 for(j=jsh+1;j<=NSD;j++){
	 fshva[j] = fdr[1];
	 fshvi[j] = fdr[2];
	 fshvp[j] = fdr[3];
     #ifdef FTYS /* TETHYS */
	 fshvt[j] = fdr[4];
     #endif	 
 }
 /* shelf/deep: modify CaCO3 rain */
 vvelm(fprva,fshva,fprva,NSD);
 vvelm(fprvi,fshvi,fprvi,NSD);
 vvelm(fprvp,fshvp,fprvp,NSD);
 /* shelf/deep: modify clay  rain */
 vvelm(frrfa,fshva,frrfa,NSD);
 vvelm(frrfi,fshvi,frrfi,NSD);
 vvelm(frrfp,fshvp,frrfp,NSD);
 #ifdef FTYS /* TETHYS */
 vvelm(fprvt,fshvt,fprvt,NSD);
 vvelm(frrft,fshvt,frrft,NSD);
 #endif
 #ifdef FSEDCC
 /* shelf/deep: modify Ca13CO3 rain */
 vvelm(fpracc,fshva,fpracc,NSD);
 vvelm(fpricc,fshvi,fpricc,NSD);
 vvelm(fprpcc,fshvp,fprpcc,NSD);
 #ifdef FTYS /* TETHYS */
 vvelm(fprtcc,fshvt,fprtcc,NSD);
 #endif
 #endif

 /* calc. saturation at sediment levels    */
 for(i=1;i<=3;i++){ /* AIP */
  for(k=1;k<=NSD;k++){ /* sed levels */
     tsedv[k] =  tcb[klid[k]+i-1];
     ssedv[k] = salv[klid[k]+i-1];
     kspfun(&kspcsed,&tmp,tsedv[k],ssedv[k],0.1*dsv[k],cac,mgc,s4c);
	 co3satm[i][k] = kspcsed/cac;  
  }
 }
#ifdef FTYS /* TETHYS */
  for(k=1;k<=NSD;k++){ /* sed levels */
     tsedv[k] =  tcb[klidt[k]];
     ssedv[k] = salv[klidt[k]];
     kspfun(&kspcsed,&tmp,tsedv[k],ssedv[k],0.1*dsv[k],cac,mgc,s4c);
	 co3satm[4][k] = kspcsed/cac;  
  }
#endif

 /* calc. omega and find supersat indices */
 for(k=1;k<=NSD;k++){
     omva[k] = co3[klid[k]]  /co3satm[1][k]; 
     omvi[k] = co3[klid[k]+1]/co3satm[2][k]; 
     omvp[k] = co3[klid[k]+2]/co3satm[3][k];
     #ifdef FTYS /* TETHYS */
     omvt[k] = co3[klidt[k]] /co3satm[4][k];
     #endif	 
 }
 dfind(kda,omva,NSD,LARGEQ,1.0,MZOM);
 dfind(kdi,omvi,NSD,LARGEQ,1.0,MZOM);
 dfind(kdp,omvp,NSD,LARGEQ,1.0,MZOM);
 #ifdef FTYS /* TETHYS */
 dfind(kdt,omvt,NSD,LARGEQ,1.0,MZOM);
 #endif	 

 /* porosities as functions of fc */
 ffphi = (phi1-phi0)/(1.-phi1);
 for(k=1;k<=NSD;k++){
     phia[k] = (phi0+ffphi*fcva[k])/(1.+ffphi*fcva[k]);
     phii[k] = (phi0+ffphi*fcvi[k])/(1.+ffphi*fcvi[k]);
     phip[k] = (phi0+ffphi*fcvp[k])/(1.+ffphi*fcvp[k]);
     #ifdef FTYS /* TETHYS */
     phit[k] = (phi0+ffphi*fcvt[k])/(1.+ffphi*fcvt[k]);
     #endif	 
 }

 /* sed rate, m/y (kg/m2/y / kg*m3 = m/y) */
 /* CaCO3 */
 vsclr(rscva,MTOKGCACR/RHOS/(1.-phi1),fprva,NSD);
 vsclr(rscvi,MTOKGCACR/RHOS/(1.-phi1),fprvi,NSD);
 vsclr(rscvp,MTOKGCACR/RHOS/(1.-phi1),fprvp,NSD);
 /* remainder/clay */
 vsclr(rsrva,1./RHOS/(1.-phi0),frrfa,NSD);
 vsclr(rsrvi,1./RHOS/(1.-phi0),frrfi,NSD);
 vsclr(rsrvp,1./RHOS/(1.-phi0),frrfp,NSD);
 /* CaCO3 + clay */
 vvsum(rsva,rscva,rsrva,NSD);
 vvsum(rsvi,rscvi,rsrvi,NSD);
 vvsum(rsvp,rscvp,rsrvp,NSD);

 #ifdef FTYS /* TETHYS */
 vsclr(rscvt,MTOKGCACR/RHOS/(1.-phi1),fprvt,NSD);
 vsclr(rsrvt,1./RHOS/(1.-phi0),frrft,NSD);
 vvsum(rsvt,rscvt,rsrvt,NSD);
 #endif	 

 #ifdef FSEDCC
 /* Sediment 13C */
 /* sed rate, m/y (kg/m2/y / kg*m3 = m/y) */
 /* Ca13CO3 */
 vsclr(rscvacc,MTOKGCACR/RHOS/(1.-phi1),fpracc,NSD); 
 vsclr(rscvicc,MTOKGCACR/RHOS/(1.-phi1),fpricc,NSD); 
 vsclr(rscvpcc,MTOKGCACR/RHOS/(1.-phi1),fprpcc,NSD); 
 #ifdef FTYS /* TETHYS */
 vsclr(rscvtcc,MTOKGCACR/RHOS/(1.-phi1),fprtcc,NSD); 
 #endif	 
 #endif


 /* dissolution, all to zero */
 fsetv(dissa,0.0,NSD);
 fsetv(dissi,0.0,NSD);
 fsetv(dissp,0.0,NSD);
 #ifdef FTYS /* TETHYS */
 fsetv(disst,0.0,NSD);
 #endif	 

 /* correction for Ca,Mg effect on kinetics */
 rcak = (cac/CAM)*(1./(1.-ALPKC*(MGM/CAM-mgc/cac)));

 /* dissolution, mol/m2/y */
 for(k=1;k<=NSD;k++){
   if(kda[k] == 0){
      dka[k] = kssd*pow((co3satm[1][k]-co3[klid[k]]  )*rcak,ncsd);
    dissa[k] = sqrt(fcva[k])*dka[k];
   }
   if(kdi[k] == 0){
      dki[k] = kssd*pow((co3satm[2][k]-co3[klid[k]+1])*rcak,ncsd);
	dissi[k] = sqrt(fcvi[k])*dki[k];
   }
   if(kdp[k] == 0){
      dkp[k] = kssd*pow((co3satm[3][k]-co3[klid[k]+2])*rcak,ncsd);
	dissp[k] = sqrt(fcvp[k])*dkp[k];
   }
   #ifdef FTYS /* TETHYS */
   if(kdt[k] == 0){
      dkt[k] = kssd*pow((co3satm[4][k]-co3[klidt[k] ])*rcak,ncsd);
	disst[k] = sqrt(fcvt[k])*dkt[k];
   }
   #endif	 
 }

 /* numerics (avoid negative fc and NaN): 
	linear drop in fc, as fc -> 0       */
 dfind(jja,fcva,NSD,SMALLER,fcsml,fcsml);
 dfind(jji,fcvi,NSD,SMALLER,fcsml,fcsml);
 dfind(jjp,fcvp,NSD,SMALLER,fcsml,fcsml);
 #ifdef FTYS /* TETHYS */
 dfind(jjt,fcvt,NSD,SMALLER,fcsml,fcsml);
 #endif

 bbf = sqrt(1./fcsml);          
 for(k=1;k<=NSD;k++){
   if(kda[k] == 0 && jja[k] != 0)
     dissa[k] = fcva[k]*dka[k]*bbf; 
   if(kdi[k] == 0 && jji[k] != 0)
     dissi[k] = fcvi[k]*dki[k]*bbf; 
   if(kdp[k] == 0 && jjp[k] != 0)
     dissp[k] = fcvp[k]*dkp[k]*bbf; 
   #ifdef FTYS /* TETHYS */
   if(kdt[k] == 0 && jjt[k] != 0)
     disst[k] = fcvt[k]*dkt[k]*bbf; 
   #endif
 }

 /* diss rate, m/y [(mol/m2/y)*kg/mol / kg*m3 = m/y] */
 /* pure calcite/(1-phi1) = Delta h1                 */
 vsclr(rdva,MTOKGCACR/RHOS/(1.-phi1),dissa,NSD); /* m/y */
 vsclr(rdvi,MTOKGCACR/RHOS/(1.-phi1),dissi,NSD); /* m/y */
 vsclr(rdvp,MTOKGCACR/RHOS/(1.-phi1),dissp,NSD); /* m/y */
 
 /* total burial rate, m/y (sed - diss) */
 vvsub(wva,rsva,rdva,NSD); 
 vvsub(wvi,rsvi,rdvi,NSD); 
 vvsub(wvp,rsvp,rdvp,NSD); 

 /* find erosion indices */ 
 dfind(lea,wva,NSD,SMALLER,0.0,MZER);
 dfind(lei,wvi,NSD,SMALLER,0.0,MZER);
 dfind(lep,wvp,NSD,SMALLER,0.0,MZER);

 #ifdef FTYS /* TETHYS */
 vsclr(rdvt,MTOKGCACR/RHOS/(1.-phi1),disst,NSD);
 vvsub(wvt,rsvt,rdvt,NSD); 
 dfind(lett,wvt,NSD,SMALLER,0.0,MZER);
 #endif

 /* calcite burial (w>0) */
 for(k=1;k<=NSD;k++){
     wcva[k] = fcva[k]*wva[k]*(1.-phia[k])/(1.-phi1);
     wcvi[k] = fcvi[k]*wvi[k]*(1.-phii[k])/(1.-phi1);
     wcvp[k] = fcvp[k]*wvp[k]*(1.-phip[k])/(1.-phi1);
	 #ifdef FTYS /* TETHYS */
     wcvt[k] = fcvt[k]*wvt[k]*(1.-phit[k])/(1.-phi1);
     #endif
 }
 /* flags default: 1, erosion: 0 */
 fsetv(fba,1.0,NSD);
 fsetv(fbi,1.0,NSD);
 fsetv(fbp,1.0,NSD);
 #ifdef FTYS /* TETHYS */
 fsetv(fbt,1.0,NSD);
 #endif

 /* calcite erosion (w<0) */
 /* Note on signs: 
	dh1/dt = -B*(-w)*.. -rrs (eqs. in model paper)
	wc     = -B*(+w)*.. +rrs (model code) because
	dh1/dt = ... - wc       
  */
 for(k=1;k<=NSD;k++){
   if(lea[k] != 0){
     wcva[k]  = -(1.-fc0a[k])*wva[k]*(1.-phiia[k])/(1.-phi0);
     wcva[k] +=  rsrva[k];
      fba[k]  = 0.0;	   
   }
   if(lei[k] != 0){
     wcvi[k]  = -(1.-fc0i[k])*wvi[k]*(1.-phiii[k])/(1.-phi0);
     wcvi[k] +=  rsrvi[k];
      fbi[k]  = 0.0;	   
   }
   if(lep[k] != 0){
     wcvp[k]  = -(1.-fc0p[k])*wvp[k]*(1.-phiip[k])/(1.-phi0);
     wcvp[k] +=  rsrvp[k];
      fbp[k]  = 0.0;	   
   }
   #ifdef FTYS /* TETHYS */
   if(lett[k] != 0){
     wcvt[k]  = -(1.-fc0t[k])*wvt[k]*(1.-phiit[k])/(1.-phi0);
     wcvt[k] +=  rsrvt[k];
      fbt[k]  = 0.0;	   
   }
   #endif 
 } 

 /* dphi/dfc */
 for(k=1;k<=NSD;k++){
      tmp = ffphi*(1.-phi0)/(1.+ffphi*fcva[k])/(1.+ffphi*fcva[k]);
  gsda[k] =   hsl*(1.-phia[k]-fcva[k]*tmp)/(1.-phi1);
      tmp = ffphi*(1.-phi0)/(1.+ffphi*fcvi[k])/(1.+ffphi*fcvi[k]);
  gsdi[k] =   hsl*(1.-phii[k]-fcvi[k]*tmp)/(1.-phi1);
      tmp = ffphi*(1.-phi0)/(1.+ffphi*fcvp[k])/(1.+ffphi*fcvp[k]);
  gsdp[k] =   hsl*(1.-phip[k]-fcvp[k]*tmp)/(1.-phi1);
  #ifdef FTYS /* TETHYS */
      tmp = ffphi*(1.-phi0)/(1.+ffphi*fcvt[k])/(1.+ffphi*fcvt[k]);
  gsdt[k] =   hsl*(1.-phit[k]-fcvt[k]*tmp)/(1.-phi1);
  #endif 
 }

 /* dissolution (w>0) in mol/m2/y, see above */

 /* dissolution (w<0) in mol/m2/y            */
 for(k=1;k<=NSD;k++){
   if(lea[k] != 0)
     dissa[k] = (rsva[k]-wva[k])*(1.-phi1)*RHOS/MTOKGCACR;
   if(lei[k] != 0)
     dissi[k] = (rsvi[k]-wvi[k])*(1.-phi1)*RHOS/MTOKGCACR;
   if(lep[k] != 0)
     dissp[k] = (rsvp[k]-wvp[k])*(1.-phi1)*RHOS/MTOKGCACR;
   #ifdef FTYS /* TETHYS */
   if(lett[k] != 0)
     disst[k] = (rsvt[k]-wvt[k])*(1.-phi1)*RHOS/MTOKGCACR;
   #endif 
 }

 /* store in dissolution matrix, mol C/y */
 for(k=1;k<=NSD;k++){
     dissm[1][k] = dissa[k]*asva[k]*ab[1];
     dissm[2][k] = dissi[k]*asvi[k]*ab[2];
     dissm[3][k] = dissp[k]*asvp[k]*ab[3];
     #ifdef FTYS /* TETHYS */
     dissm[4][k] = disst[k]*asvt[k]*ab[11];
     #endif 
 }

 #ifdef FSEDCC
 /* Sediment 13C */
 /* dissolution rate */
 vvelm(rdvacc,rccsa,rdva,NSD);
 vvelm(rdvicc,rccsi,rdvi,NSD);
 vvelm(rdvpcc,rccsp,rdvp,NSD);
 /* burial rate */
 vvelm(wcvacc,rccsa,wcva,NSD);
 vvelm(wcvicc,rccsi,wcvi,NSD);
 vvelm(wcvpcc,rccsp,wcvp,NSD);
 #ifdef FTYS /* TETHYS */
 vvelm(rdvtcc,rccst,rdvt,NSD);
 vvelm(wcvtcc,rccst,wcvt,NSD);
 #endif 

 /* dissolution matrix, mol C/y */
 for(k=1;k<=NSD;k++){
     dissmcc[1][k] = rccsa[k]*dissa[k]*asva[k]*ab[1];
     dissmcc[2][k] = rccsi[k]*dissi[k]*asvi[k]*ab[2];
     dissmcc[3][k] = rccsp[k]*dissp[k]*asvp[k]*ab[3];
     #ifdef FTYS /* TETHYS */
     dissmcc[4][k] = rccst[k]*disst[k]*asvt[k]*ab[11];
     #endif 
 }
 #endif

#endif

 /* conveyor */	
 gtha = (1-thbra);       /* export into Deep Ind */
 gthi = (1-thbra-thbri); /* export into Deep Pac */

 /*=========== Long C-Cycle fluxes ==============*/
 /* set volcanic degassing */
 fvc = 1.0*fvc0;
 /* CaSiO3 weathering */
 fsi  =  fvc0*pow(pco2a/pcsi,nsi); /* CaSiO3 */
 /* CaCO3  weathering */
 finc = finc0*pow(pco2a/pcsi,ncc); /* CaCO3  */
 
    
    /*DP set sulfate degassing*/
    /* think about continuously variable pulses */
    /* multiply Gt by 3.125e13 */
 /* index = (int) floor((t-t0)/TIMESTEP);
 fso4 = so4array[index]; */
 
     if(t > 390.e3 && t < 390001.) {
        fso4 = 1.015625e16/AOC; /*Boundary emissions*/
    }
    else if(t > 32.e3 && t < 172.e3) {
        fso4 = 1.650670e12/AOC; /* 1st pulse, 7395GtS = 2.31094E+17molS over 140kyr in molS/m2/y*/
    }
    
    else if(t > 390001. && t < 745.e3) {
        fso4 = 9.727113e10/AOC; /*2nd pulse, 1105GtS = 3.45313E+16molS over 355kyr*/
    }
    else {
        fso4 = 0;
    }


    
 /*=========== Long C-Cycle fluxes 13C ===========*/
 fvccc  = rvccc*fvc;
 fsicc  = rincc*fsi;
 fincc  = rincc*finc;
 fkrgcc = rkrgcc*fkrg;               /* oxidation */
 vsclr(fkrgccb,alpcorg*fkrg,rccb,NB); /* burial   */

 /* sediment flags */	
#ifdef FSED
 vsclr(exlv,1.0,eplv,NLS);
 exh   = eph;	
if(NOCT >= 6){	
 vsclr(exlvcc,1.0,eplvcc,NLS);
 exhcc = ephcc;	
}
 sdn   = 0.;
#else
 vsclr(exlv,1.0,eclv,NLS);
 exh   = ech;
if(NOCT >= 6){	
 vsclr(exlvcc,1.0,eclvcc,NLS);
 exhcc = echcc;	
}
 sdn   = 1.;
#endif		 


/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Here is the right-hand side of the DEQs
%
%     Units ocean tracer: [mol/m3/y]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

 /*=================== DIC =======================*/
 /* TH and mixing */
 thmfun(dic,dicp,fconv,thc,tso,tto,mxv,mhd,vb,gtha,thbra,gthi,thbri);
 /* air-sea */
if(NCATM == 1){	
 for(k=1;k<=NTS;k++){
	kk = kkv[k];
 	dicp[kk] += kasv[kk]*(pco2a-pco2[kk])/vb[kk];
 }
}	
/* bio pump Corg */
 for(k=1;k<=NOC;k++){
     dicp[m1[k]]   +=         -eclv[k]/vb[m1[k]]  ; /* L                 */
     dicp[m2[k]+3] +=     frei*exlv[k]/vb[m2[k]+3]; /* I #!  EC or EP    */
     dicp[m3[k]+6] +=       oi*exlv[k]/vb[m3[k]+6]; /* D #!(tot or Corg) */
     dicp[m3[k]+6] += 0.5*nuwd*ealv[k]/vb[m3[k]+6]; /* D ClmnDiss        */
 }	
 dicp[10] -=  ech/vb[10];                   /* H            */
 for(k=7;k<=9;k++){                         /* DA,DI,DP     */
     dicp[k]   +=          exh/vb[k]*gp[k]; /* #!  EC or EP */
     dicp[k]   += 0.5*nuwd*eah/vb[k]*gp[k]; 
 }   
#ifdef FSED 
 /* riverine and sediment fluxes */
 for(k=1;k<=NOC;k++){
     dicp[m1[k]] += 2.*finc*AOC/vb[m1[k]]/NOC;         /* Fin  L wthr:2 */
     dicp[m1[k]] += 2.*fsi *AOC/vb[m1[k]]/NOC;         /* ... #! Si     */
     dicp[m1[k]] -= 1.*fkrg*AOC/vb[m1[k]]/NOC;         /*     #! krgn   */
   dicp[m1[k]]   += msumlh(dissm,NOC,k,NSD,1,    ns1)/vb[m1[k]];   /* diss L  */
   dicp[m2[k]+3] += msumlh(dissm,NOC,k,NSD,1+ns1,ns2)/vb[m2[k]+3]; /* diss I  */
   dicp[m3[k]+6] += msumlh(dissm,NOC,k,NSD,1+ns2,NSD)/vb[m3[k]+6]; /* diss D  */
 }
#endif
	
 /*=================== ALK   ======================*/
 /* TH and mixing   */
 thmfun(alk,alkp,fconv,thc,tso,tto,mxv,mhd,vb,gtha,thbra,gthi,thbri);
 /* bio pump CaCO3, Aorg */
 for(k=1;k<=NOC;k++){
   alkp[m1[k]]   +=          -ealv[k]/vb[m1[k]]  +enlv[k]/vb[m1[k]]   ; 
   alkp[m2[k]+3] += frei*(sdn*ealv[k]/vb[m2[k]+3]-enlv[k]/vb[m2[k]+3]); 
   alkp[m3[k]+6] +=   oi*(sdn*ealv[k]/vb[m3[k]+6]-enlv[k]/vb[m3[k]+6]);
   alkp[m3[k]+6] +=	     nuwd*ealv[k]/vb[m3[k]+6]  ; /* D ClmnDiss */
 }
 alkp[10] += -eah/vb[10] + enh/vb[10];
 for(k=7;k<=9;k++){                           /* DA,DI,DP   */
     alkp[k] +=  (sdn*eah/vb[k]-enh/vb[k])*gp[k];  
     alkp[k] +=  nuwd*eah/vb[k]*gp[k];
 }
#ifdef FSED 
 /* riverine and sediment fluxes */
 for(k=1;k<=NOC;k++){
    alkp[m1[k]] += 2.*finc*AOC/vb[m1[k]]/NOC;         /* Fin  L wthr:2 */
    alkp[m1[k]] += 2.*fsi *AOC/vb[m1[k]]/NOC;         /* ... #! Si     */
     alkp[m1[k]] -= 2.*fso4 *AOC/vb[m1[k]]/NOC;  /*sulfate:TA = 1:-2*/
  alkp[m1[k]]   += 2.*msumlh(dissm,NOC,k,NSD,1,    ns1)/vb[m1[k]];   /* diss L */
  alkp[m2[k]+3] += 2.*msumlh(dissm,NOC,k,NSD,1+ns1,ns2)/vb[m2[k]+3]; /* diss I */
  alkp[m3[k]+6] += 2.*msumlh(dissm,NOC,k,NSD,1+ns2,NSD)/vb[m3[k]+6]; /* diss D */
 }
#endif
	
 /*=================== PO4   ======================*/
 /* TH and mixing   */
 thmfun(po4,po4p,fconv,thc,tso,tto,mxv,mhd,vb,gtha,thbra,gthi,thbri);
 /* bio pump Porg */
 for(k=1;k<=NOC;k++){
     po4p[m1[k]]   -=      pplv[k]/vb[m1[k]]  ;
     po4p[m2[k]+3] += frei*pplv[k]/vb[m2[k]+3];
     po4p[m3[k]+6] +=   oi*pplv[k]/vb[m3[k]+6];
 }
 po4p[10] -= pph/vb[10];
 for(k=7;k<=9;k++)
     po4p[k]   +=      pph/vb[k]*gp[k]; /* DA,DI,DP */
 
 /*=================== temperature ================*/
 /* Ocean temperature change, using simple         */
 /* climate sensitivity (Archer, 2005)             */
 /* TAU: Relax times (y) Surf, Intrm, Deep         */

if(NOCT >= 4){	
  if(tsnsflag == 1 && pco2a >= 150.){	
    tmp = sclim*log(pco2a/pcsi)/LOG2;
    for(k=1;k<=NOC;k++){
      tcbp[m1[k]]   = (tcb0[m1[k]]  +tmp-tcb[m1[k]]  )/TAUS;
      tcbp[m2[k]+3] = (tcb0[m2[k]+3]+tmp-tcb[m2[k]+3])/TAUI;
      tcbp[m3[k]+6] = (tcb0[m3[k]+6]+tmp-tcb[m3[k]+6])/TAUD;
    }
    tcbp[10] = (tcb0[10]+tmp-tcb[10])/TAUS;
  }
}


 /*=================== Oxygen =====================*/
if(NOCT >= 5){	
 /* TH and mixing   */
 thmfun(dox,doxp,fconv,thc,tso,tto,mxv,mhd,vb,gtha,thbra,gthi,thbri);
  /* air-sea, o2sat = O2 at saturation (atm) */
 for(k=1;k<=NTS;k++){
	 kk = kkv[k];
     doxp[kk] += vask[kk]*(o2sat[kk]-dox[kk])/vb[kk];
 }	
 /* bio pump O2 */
 for(k=1;k<=NOC;k++){
   doxp[m1[k]]   +=                    eplv[k]*REDO2C/vb[m1[k]]  ; /* L */
   doxp[m2[k]+3] -= frei*fmmo[m2[k]+3]*eplv[k]*REDO2C/vb[m2[k]+3]; /* I */
   doxp[m3[k]+6] -=   oi*fmmo[m3[k]+6]*eplv[k]*REDO2C/vb[m3[k]+6]; /* D */
 }	
 doxp[10] +=  eph*REDO2C/vb[10];               /* H */
 for(k=7;k<=9;k++){                            /* DA,DI,DP */
     doxp[k]   -= gp[k]*fmmo[k]*eph*REDO2C/vb[k]; 
 }   	
}


 /*=================== DICC: 13C ==================*/
if(NOCT >= 6){	
 /* TH and mixing */
 thmfun(dicc,diccp,fconv,thc,tso,tto,mxv,mhd,vb,gtha,thbra,gthi,thbri);
 /* air-sea */
if(NCCATM == 1){	
 for(k=1;k<=NTS;k++){
	 kk = kkv[k];
	tmp = (alpdg[kk]*pcco2a-alpdb[kk]*rccb[kk]*pco2[kk]);
 	diccp[kk] += kasv[kk]*alpu*tmp/vb[kk];
 }
}	
 /* bio pump Corg */
 for(k=1;k<=NOC;k++){
     diccp[m1[k]]   +=         -eclvcc[k]/vb[m1[k]]  ; /* L                 */
     diccp[m2[k]+3] +=     frei*exlvcc[k]/vb[m2[k]+3]; /* I #!  EC or EP    */
     diccp[m3[k]+6] +=       oi*exlvcc[k]/vb[m3[k]+6]; /* D #!(tot or Corg) */
     diccp[m3[k]+6] += 0.5*nuwd*ealvcc[k]/vb[m3[k]+6]; /* D ClmnDiss        */
 }	
 diccp[10] -=  echcc/vb[10];                   /* H            */
 for(k=7;k<=9;k++){                            /* DA,DI,DP     */
     diccp[k]   +=          exhcc/vb[k]*gp[k]; /* #!  EC or EP */
     diccp[k]   += 0.5*nuwd*eahcc/vb[k]*gp[k]; 
 }
#ifdef FSEDCC 
 /* riverine and sediment fluxes */
 for(k=1;k<=NOC;k++){
    diccp[m1[k]] += 2.*fincc         *AOC/vb[m1[k]]/NOC;   /* Fin  L wthr:2 */
    diccp[m1[k]] += 2.*fsicc         *AOC/vb[m1[k]]/NOC;   /* ... #! Si     */
    diccp[m1[k]] -= 1.*fkrgccb[m1[k]]*AOC/vb[m1[k]]/NOC;   /*     #! krgn   */
  diccp[m1[k]]   += msumlh(dissmcc,NOC,k,NSD,1,    ns1)/vb[m1[k]];   /* diss L */
  diccp[m2[k]+3] += msumlh(dissmcc,NOC,k,NSD,1+ns1,ns2)/vb[m2[k]+3]; /* diss I */
  diccp[m3[k]+6] += msumlh(dissmcc,NOC,k,NSD,1+ns2,NSD)/vb[m3[k]+6]; /* diss D */
 }
#endif

} /* NOCT */

 /*=================== C  atm =====================*/
if(NCATM == 1){	
 for(k=1;k<=NTS;k++)
 	catmp -= kasv[kkv[k]]*(pco2a-pco2[kkv[k]])/AOC; /* mol/m2/y */

#ifdef FSED
    catmp +=  fvc - finc - 2.*fsi + fkrg;           /* wthr #!  */
#endif
	
 /* fossil fuel */
 if(ffflag == 1)
    catmp += fems/AOC; /* input(t) (mol/y)-> mol/m2/y */   
	
} /* NCATM */

 /*=================== 13C  atm ===================*/
if(NCCATM == 1){	
 for(k=1;k<=NTS;k++){
         kk = kkv[k];
	 	tmp = (alpdg[kk]*pcco2a-alpdb[kk]*rccb[kk]*pco2[kk]);
 	ccatmp -= kasv[kk]*alpu*tmp/AOC; /* mol/m2/y */
 }

#ifdef FSED
    ccatmp +=  fvccc - fincc - 2.*fsicc + fkrgcc;           /* wthr #!  */
#endif	
	
 /* fossil fuel */
 if(ffflag == 1)
    ccatmp += femscc/AOC; /* input(t) (mol/y)-> mol/m2/y */   
	
} /* NCCATM */	


/*========= carbon input scenarios ============*/	
if(cinpflag){	
 if(t > tcin0 && t < tcin0+tcspan){
     catmp +=        cinp*1.e15/12./AOC/tcspan;
    ccatmp += rccinp*cinp*1.e15/12./AOC/tcspan;
 }
}

/*
if(t > 1700. && t < 1700.+5.)
    catmp += 3000.e15/12./AOC/5.;
*/
/*
if(t > 2050.)
     catmp -= 1*2.e15/12./AOC;
*/

 /********* all into yp *********/	
 /* ocean tracers */	
 for(k=1;k<=NB;k++){
    yp[k]      = dicp[k];
    yp[k+1*NB] = alkp[k];
    yp[k+2*NB] = po4p[k]/MMTOM; /* convert back to mmol/m3 */
 }

if(NOCT >= 4){	
 for(k=1;k<=NB;k++)
    yp[k+3*NB] = tcbp[k]/TSCAL; /* Temp: scale to order 1  */
}

if(NOCT >= 5){	
 for(k=1;k<=NB;k++)
    yp[k+4*NB] = doxp[k];       /* mol/m3 */
}

if(NOCT >= 6){	
 for(k=1;k<=NB;k++)
    yp[k+5*NB] = diccp[k]/CNTI; /* scaling /CNTI !!!       */
}


if(NCATM == 1){	
 /* C atm */	
 yp[NOCT*NB+1] = catmp*CNTI;    /* scaling *CNTI !!!       */
}

if(NCCATM == 1){	
 /* 13C atm */	
 yp[NOCT*NB+2] = ccatmp;        /* not scaled    !!!       */
}

#ifdef FSED
 /*============== Sediment Boxes =================*/

 for(k=1;k<=NSD;k++){
     fcvap[k]  = fba[k]*rscva[k] - fba[k]*rdva[k] - wcva[k];
     fcvap[k] /= gsda[k];
     fcvip[k]  = fbi[k]*rscvi[k] - fbi[k]*rdvi[k] - wcvi[k];
     fcvip[k] /= gsdi[k];
     fcvpp[k]  = fbp[k]*rscvp[k] - fbp[k]*rdvp[k] - wcvp[k];
     fcvpp[k] /= gsdp[k];
     #ifdef FTYS /* TETHYS */ 
     fcvtp[k]  = fbt[k]*rscvt[k] - fbt[k]*rdvt[k] - wcvt[k];
     fcvtp[k] /= gsdt[k];
     #endif	 
 }

 /************ all into yp *********/
 /* NOATM = (NOCT*NB+NCATM+NCCATM) */
 for(k=1;k<=NSD;k++){
    yp[NOATM+k]       = fcvap[k];
    yp[NOATM+k+1*NSD] = fcvip[k];
    yp[NOATM+k+2*NSD] = fcvpp[k];
    #ifdef FTYS /* TETHYS */ 
    yp[NOATM+k+3*NSD] = fcvtp[k];
    #endif
 }

 #ifdef FSEDCC
 /*============== Sediment Boxes 13C ==============*/

 for(k=1;k<=NSD;k++){
     fccvap[k]  = fba[k]*rscvacc[k] - fba[k]*rdvacc[k] - wcvacc[k];
     fccvap[k] /= gsda[k];
     fccvip[k]  = fbi[k]*rscvicc[k] - fbi[k]*rdvicc[k] - wcvicc[k];
     fccvip[k] /= gsdi[k];
     fccvpp[k]  = fbp[k]*rscvpcc[k] - fbp[k]*rdvpcc[k] - wcvpcc[k];
     fccvpp[k] /= gsdp[k];
	 #ifdef FTYS /* TETHYS */ 
     fccvtp[k]  = fbt[k]*rscvtcc[k] - fbt[k]*rdvtcc[k] - wcvtcc[k];
     fccvtp[k] /= gsdt[k];
     #endif	 
 }


 /************ all into yp *********/
 /* NOATM = (NOCT*NB+NCATM+NCCATM) */
 for(k=1;k<=NSD;k++){
    yp[NOATM+k+(3+KTY)*NSD] = fccvap[k]/CNTI; /* scaling /CNTI */
    yp[NOATM+k+(4+KTY)*NSD] = fccvip[k]/CNTI; /* scaling /CNTI */
    yp[NOATM+k+(5+KTY)*NSD] = fccvpp[k]/CNTI; /* scaling /CNTI */
    #ifdef FTYS /* TETHYS */ 
    yp[NOATM+k+(6+KTY)*NSD] = fccvtp[k]/CNTI; /* scaling /CNTI */
    #endif
 }

 #endif /* FSEDCC */
#endif /* FSED */


 /* free locals */
 free_dvector(eplv,1,NLS);
 free_dvector(ealv,1,NLS);
 free_dvector(enlv,1,NLS);
 free_dvector(pplv,1,NLS);
 free_dvector(eclv,1,NLS);
 free_dvector(exlv,1,NLS);
 free_dvector(co2,1,NB);
 free_dvector(pco2,1,NB);
 free_dvector(co3,1,NB);
 free_dvector(hpls,1,NB);
 free_dvector(ph,1,NB);
 free_dvector(kh,1,NB);
 free_dvector(o2sat,1,NB);
 free_dvector(kasv0,1,NB);
 free_dvector(kasv,1,NB);
 free_dvector(vask0,1,NB);
 free_dvector(vask,1,NB);
 free_dvector(fmmo,1,NB);
 free_dvector(rccb,1,NB);
 free_dvector(alpdb,1,NB);
 free_dvector(alpdg,1,NB);
 free_dvector(alpcb,1,NB);
 free_dvector(eplvcc,1,NLS);
 free_dvector(ealvcc,1,NLS);
 free_dvector(eclvcc,1,NLS);
 free_dvector(exlvcc,1,NLS);
 free_dvector(fkrgccb,1,NB);


 /* free renamed */
 free_dvector(dic,1,NB);
 free_dvector(alk,1,NB);
 free_dvector(po4,1,NB);
 free_dvector(tcb,1,NB);
 free_dvector(dox,1,NB);
 free_dvector(dicc,1,NB);
 free_dvector(dicp,1,NB);
 free_dvector(alkp,1,NB);
 free_dvector(po4p,1,NB);
 free_dvector(tcbp,1,NB);
 free_dvector(doxp,1,NB);
 free_dvector(diccp,1,NB);
#ifdef FSED
 free_dvector(fcva,1,NSD);
 free_dvector(fcvi,1,NSD);
 free_dvector(fcvp,1,NSD);
 free_dvector(fcvap,1,NSD);
 free_dvector(fcvip,1,NSD);
 free_dvector(fcvpp,1,NSD);
 #ifdef FTYS
 free_dvector(fcvt,1,NSD);
 free_dvector(fcvtp,1,NSD);
 #endif
 #ifdef FSEDCC
 free_dvector(fccva,1,NSD);
 free_dvector(fccvi,1,NSD);
 free_dvector(fccvp,1,NSD);
 free_dvector(fccvap,1,NSD);
 free_dvector(fccvip,1,NSD);
 free_dvector(fccvpp,1,NSD);
 #ifdef FTYS
 free_dvector(fccvt,1,NSD);
 free_dvector(fccvtp,1,NSD);
 #endif
 #endif
 free_dvector(fprva,1,NSD);
 free_dvector(fprvi,1,NSD);
 free_dvector(fprvp,1,NSD);
 free_dvector(frrfa,1,NSD);
 free_dvector(frrfi,1,NSD);
 free_dvector(frrfp,1,NSD);

 free_dvector(tsedv,1,NSD);
 free_dvector(ssedv,1,NSD);

 free_dmatrix(co3satm,1,NOC,1,NSD);

 free_dvector(omva,1,NSD);
 free_dvector(omvi,1,NSD);
 free_dvector(omvp,1,NSD);

 free_ivector(kda,1,NSD);
 free_ivector(kdi,1,NSD);
 free_ivector(kdp,1,NSD);

 free_dvector(phia,1,NSD);
 free_dvector(phii,1,NSD);
 free_dvector(phip,1,NSD);

 free_dvector(rscva,1,NSD);
 free_dvector(rscvi,1,NSD);
 free_dvector(rscvp,1,NSD);

 free_dvector(rsrva,1,NSD);
 free_dvector(rsrvi,1,NSD);
 free_dvector(rsrvp,1,NSD);

 free_dvector(rsva,1,NSD);
 free_dvector(rsvi,1,NSD);
 free_dvector(rsvp,1,NSD);

 free_dvector(dka,1,NSD);
 free_dvector(dki,1,NSD);
 free_dvector(dkp,1,NSD);

 free_dvector(dissa,1,NSD);
 free_dvector(dissi,1,NSD);
 free_dvector(dissp,1,NSD);

 free_dvector(rdva,1,NSD);
 free_dvector(rdvi,1,NSD);
 free_dvector(rdvp,1,NSD);

 free_dvector(wva,1,NSD);
 free_dvector(wvi,1,NSD);
 free_dvector(wvp,1,NSD);

 free_ivector(lea,1,NSD);
 free_ivector(lei,1,NSD);
 free_ivector(lep,1,NSD);

 free_dvector(wcva,1,NSD);
 free_dvector(wcvi,1,NSD);
 free_dvector(wcvp,1,NSD);

 free_dvector(fba,1,NSD);
 free_dvector(fbi,1,NSD);
 free_dvector(fbp,1,NSD);

 free_dvector(gsda,1,NSD);
 free_dvector(gsdi,1,NSD);
 free_dvector(gsdp,1,NSD);

 free_ivector(jja,1,NSD);
 free_ivector(jji,1,NSD);
 free_ivector(jjp,1,NSD);

 free_dvector(fdr,1,NOC);
 free_dvector(fshva,1,NSD);
 free_dvector(fshvi,1,NSD);
 free_dvector(fshvp,1,NSD);

 #ifdef FTYS
  free_dvector(fprvt,1,NSD);
  free_dvector(frrft,1,NSD);
  free_dvector(omvt,1,NSD);
  free_ivector(kdt,1,NSD);
  free_dvector(phit,1,NSD);
  free_dvector(rscvt,1,NSD);
  free_dvector(rsrvt,1,NSD);
  free_dvector(rsvt,1,NSD);
  free_dvector(dkt,1,NSD);
  free_dvector(disst,1,NSD);
  free_dvector(rdvt,1,NSD);
  free_dvector(wvt,1,NSD);
  free_ivector(lett,1,NSD);
  free_dvector(wcvt,1,NSD);
  free_dvector(fbt,1,NSD);
  free_dvector(gsdt,1,NSD);
  free_ivector(jjt,1,NSD);
  free_dvector(fshvt,1,NSD);
 #endif

 free_dmatrix(dissm,1,NOC,1,NSD);
 #ifdef FSEDCC
 free_dvector(rccsa,1,NSD);
 free_dvector(rccsi,1,NSD);
 free_dvector(rccsp,1,NSD);

 free_dvector(fpracc,1,NSD);
 free_dvector(fpricc,1,NSD);
 free_dvector(fprpcc,1,NSD);

 free_dvector(rscvacc,1,NSD);
 free_dvector(rscvicc,1,NSD);
 free_dvector(rscvpcc,1,NSD);

 free_dvector(rdvacc,1,NSD);
 free_dvector(rdvicc,1,NSD);
 free_dvector(rdvpcc,1,NSD);

 free_dvector(wcvacc,1,NSD);
 free_dvector(wcvicc,1,NSD);
 free_dvector(wcvpcc,1,NSD);
 #ifdef FTYS
  free_dvector(rccst,1,NSD);
  free_dvector(fprtcc,1,NSD);
  free_dvector(rscvtcc,1,NSD);
  free_dvector(rdvtcc,1,NSD);
  free_dvector(wcvtcc,1,NSD);
 #endif

 free_dmatrix(dissmcc,1,NOC,1,NSD);
 #endif
#endif

}
/*============================================================*/
/*==================== derivs() END ==========================*/
/*============================================================*/


/*============================================================*/
/*==================== jacobn() ==============================*/
/*============================================================*/
/* Numerical Jacobian                                         */
/*
   dydx = [f1 f2 ... fn]  (output from derivs) i=1,..,n

	            k ----> var
           |  df1/dy1 df1/dy2 ... df1/dyn  |   i
           |    .                    .     |   |
      J =  |    .                    .     |   |   = dfdy[i][k]
	       |    .                    .     |   v     
           |  dfn/dy1 dfn/dy2 ... dfn/dyn  |  eqn

   Call derivs twice (y[i],y[i]+epsj) for each i (=> n x n values).
   Then dfdy ~= [dydx(y+epsj)-dydx(y)]/epsj
   Note: derivs call corresponds to one COLUMN of J (e.g. dfi/dy1)
*/

void jacobn(double x,double *y,double *dfdx,double **dfdy,int n)
{
 int i,k;
 double *ye,*dydx,*dydxe,epsj;

 ye=dvector(1,n);     /* y+epsj       */
 dydx=dvector(1,n);   /* dydx(y)      */
 dydxe=dvector(1,n);  /* dydx(y+epsj) */

 /* If RHS is NOT a function of x, then dfdx=0. Check! */
 for (i=1;i<=n;i++) dfdx[i]=0.0;
	
 /* get dydx(y) */
 derivs(x,y,dydx);
 /*for(i=1;i<=n;i++)
    printf("x=%10.7f y =%10.7f dydx =%10.7f\n",x,y[i],dydx[i]); */

 /* ye: make a copy of y */
 for(i=1;i<=n;i++)
    ye[i] = y[i];

 /* ye = y + epsj. note dfdy[i=eqn][k=var]  */
 for(k=1;k<=n;k++){     /* var */
	/* epsj   = DMAX(ye[k]*epslvr,epslvr); */
	epsj   = epslvr; 
    ye[k] += epsj;
	derivs(x,ye,dydxe); /* get dydx(y+epsj) */
    for(i=1;i<=n;i++){  /* eqn */
      dfdy[i][k] = (dydxe[i]-dydx[i])/epsj; /* given k => 1 column of J */
      /*printf("x=%10.7f ye=%10.7f dydxe=%10.7f\n",x,ye[i],dydxe[i]);*/
	}
	ye[k] = y[k];  /* reset ye[k] ! next call of derivs  */
                   /* needs ye[i != k] = y[i]            */
 }
 free_dvector(ye,1,n);
 free_dvector(dydx,1,n);
 free_dvector(dydxe,1,n);
}
/*============================================================*/
/*==================== jacobn() END ==========================*/
/*============================================================*/

 /* error output */
 void ferrwrt(char error_text[])
 {
 fflush(stdout);
 fflush(stderr);
 fprintf(stderr,"\n");	 
 fprintf(stderr,"@ ======================== Error.\n");
 fprintf(stderr,"%s\n",error_text);
 if(kount > 0){
   /* write current data to files */
   fprintf(stderr,"Writing output (.dat) files for error check\n");
   writedat(tmv,yy,kount,tcb0,spm,cac,mgc,s4c,dsv,zv,nz,svrestart,fpsvstr);
 } else {
   fprintf(stderr,"Nothing to write yet.\n");
  }
 fprintf(stderr,"@ ====================== Exiting.\n");
 fprintf(stderr,"  Errare humanum est.\n");
 exit(1);
 }

