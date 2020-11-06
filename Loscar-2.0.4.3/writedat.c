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
#include <string.h>
#include <math.h>
#include "defs.h"
#include "utils.h"

/*============================================================*/
/*==================== writedat() ============================*/
/*============================================================*/
/* writedat(tmv,yy,lt,...)                                     */
/* 

 updates: 

  09/03/18 new K* corrections (ZeebeTyrrell18)
  10/22/11 calc Tethys CCD
  10/17/11 include Tethys
  10/08/11 add Ca,Mg to fcsys(). now writedat(..,Ca,Mg,..)
  09/16/11 include 13C (dicc).
  09/15/11 issue warning if PO4, O2 was low during run.
  09/11/11 include dissolved O2 (dox).
  04/21/11 reversed search order in fminv(..,1->2) 
  04/15/11 deleted write tcm.dat. tcb now tracer (NEQ >= 4).
           fixed CCD jumps (fccderr=2, see finterp() in utils.c)
  04/10/11 hgss -> hgssvl[k]. file name: pco2.dat -> pco2a.dat
  02/22/11 new file
 
   tmv.dat (time)
	  
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

   Indices:   n   i    k
             var time box

   Example 
             rows   columns    (boxes = columns)
   dic.dat [1...lt][1...NB]
   ======================
        box1  box2 ... boxNB
     t1 c11   c12      c1NB
     t2 c21   c22      c2NB    
     ...
     lt clt1  clt2     cltNB    

	   
   alk.dat [1...lt][1...NB]
   ======================
        box1  box2 ... boxNB
     t1 a11   a12      a1NB
     t2 a21   a22      a2NB    
     ...
     lt alt1  alt2     altNB    

   However, solution matrix is

   yy[1...NOCEQ][1...lt], NOCEQ = NOCT*NB

   dic.dat
   write lines for all boxes and times
   k=1...NB, i=1...lt
   yy[1...NB][1...lt]
   yy[k][i]                 ny = k

   alk.dat
   write lines for all boxes and times
   k=1...NB, i=1...lt
   yy[NB+1...2NB][1...lt],  ny = k + NB
	   
   po4.dat
   yy[2NB+1...3NB][1...lt], ny = k + 2NB

	...
                      =>    ny = k + (n-1)*NB
*/

void writedat(double *tmv,double **yy,int lt,double *tcb0,
              double **spm,double cacl,double mgcl,double s4cl,
              double *dsv,double *zv,int nz,
              int svrestart,char *fpsvstr)
{ 
 int i,n,k,ny;
 double pco2ocn,co3,hpls,ph,*hgssvl,kspc,kspa,***yyy,**tmprm,tmp;
#ifdef FSED
 int *jccda,*jccdi,*jccdp,**jccdm; 
 double fccd,fccderr,*fca,*fci,*fcp,**fcam,**fcim,**fcpm,
	    *ya,*yi,*yp; 
 #ifdef FTYS
 int *jccdt;
 double *fct,**fctm,*yt;
 #endif
#endif
	
 char fpstr[BUFSIZ]; /* max length of file name*/
 char mssg[BUFSIZ];		
 FILE *fpout;
	
 tmprm  = dmatrix(1,NB,1,lt);	
 yyy    = d3tensor(1,NB,1,lt,1,NCSWRT);	
 hgssvl = dvector(1,NB);	
 
#ifdef FSED
 fca   = dvector(1,NSD);	
 fci   = dvector(1,NSD);	
 fcp   = dvector(1,NSD);	
 fcam  = dmatrix(1,NSD,1,lt);	
 fcim  = dmatrix(1,NSD,1,lt);	
 fcpm  = dmatrix(1,NSD,1,lt);	
 ya    = dvector(1,nz);	
 yi    = dvector(1,nz);	
 yp    = dvector(1,nz);	
 jccda = ivector(1,lt);	
 jccdi = ivector(1,lt);	
 jccdp = ivector(1,lt);	
 jccdm = imatrix(1,lt,1,NOC);
#ifdef FTYS		
 fct   = dvector(1,NSD);	
 fctm  = dmatrix(1,NSD,1,lt);		
 yt    = dvector(1,nz);	
 jccdt = ivector(1,lt);	
#endif
#endif
	

if(NOCT >= 4){		
 /* 
  convert temperature units (write restart values in degC)
	
        solver      output
  ==================================
  TMP: degC/TSCAL => degC     *TSCAL
	  
 */	
 for(i=1;i<=lt;i++){ /* time */
	for(k=1;k<=NB;k++){ /* box */
       yy[k+3*NB][i] *= TSCAL; 
	}
 }
} /* NOCT */
	
 /*=============== save restart (yfinal) ==============*/
 if(svrestart == 1){

 fpout = fopen(fpsvstr,"w");
 if(fpout == NULL){
    sprintf(mssg,"writedat(): Can't open file to write: '%s'",fpsvstr);
	ferrx(mssg);
 } else {
	printf("\n@ Saving restart values in '%s'\n",fpsvstr);
 }
		 
 for(k=1;k<=NEQ;k++)
	    fprintf(fpout,"%4d %18.15f \n",k,yy[k][lt]);	 
 fclose(fpout);
	 
 } /* svrestart */
	
 /*============ write tmv.dat (time vector) ===========*/
 fpout = fopen("tmv.dat","w");
 if(fpout == NULL) 
	ferrx("writedat(): Can't open file to write: tmv.dat");
 for(i=1;i<=lt;i++)
    fprintf(fpout,"%18.15f \n",tmv[i]);
 fclose(fpout);	

 /*============ check HL-PO4 and dox ==================*/
 for(i=1;i<=lt;i++){ 
   if(yy[3*10][i] < PMMK/MMTOM){
    sprintf(mssg,"writedat(): High-Lat PO4 was low during this run.");
    sprintf(mssg,"%s\nActual HL-Corg export may have been less than CBIOH.",mssg);
    sprintf(mssg,"%s\nPO4 Uptake reduced using Michaelis-Menten kinetics.",mssg);
    fwarn(mssg);
	break;
   }
 }	
if(NOCT >= 5){		
 for(i=1;i<=lt;i++){ /* time */
  for(k=1;k<=NB;k++){ /* box */
   if(yy[k+4*NB][i] < 10.*KMMOX){
    sprintf(mssg,"writedat(): Dissolved O2 was low during this run.");
    sprintf(mssg,"%s\nCheck output. Oxygen demand exceeded supply.",mssg);
    sprintf(mssg,"%s\nO2 Uptake reduced using Michaelis-Menten kinetics.",mssg);
    fwarn(mssg);
	i = lt;
	break;
   }
  }
 }
}	

 /*=============== write OCN tracer.dat ===============*/
 /* convert ocean tracer units 
	(for start units, see initstart)
	
        solver        output
 1  DIC:  mol/m3     => mmol/kg: *(1.e3/RHO)
 2  ALK:  mol/m3     => mmol/kg: *(1.e3/RHO)
 3  PO4: mmol/m3     => umol/kg: *(1.e3/RHO) in derivs: => mol/m3 !!!
 4  TMP: degC/TSCAL  =>    degC: *TSCAL (see above)
 5  DOX:  mol/m3     =>  mol/m3  *1
 6 DICC: mol/m3/CNTI => mmol/kg: *(1.e3/RHO)*CNTI
	  
 */


 for(i=1;i<=lt;i++){ /* time */
	for(k=1;k<=NB;k++){ /* box */
       yy[k][i]      *= (1.e3/RHO); 
       yy[k+1*NB][i] *= (1.e3/RHO); 
       yy[k+2*NB][i] *= (1.e3/RHO); 
	}
 }	

if(NOCT >= 6){	
 /* dicc */
 for(i=1;i<=lt;i++){ /* time */
	for(k=1;k<=NB;k++){ /* box */
       yy[k+5*NB][i] *= (1.e3/RHO)*CNTI; 
	}
 }	
}

 /* write to file */
 for(n=1;n<=NOCT;n++){ /* output var */
   strcpy(fpstr,getocnstr(n));
   fpout = fopen(fpstr,"w");
   if(fpout == NULL){
     sprintf(mssg,"writedat(): Can't open file to write: %s",fpstr);
     ferrx(mssg);
   }
   for(i=1;i<=lt;i++){ /* time */
     for(k=1;k<=NB;k++){ /* box */
		ny = k+(n-1)*NB;
	    fprintf(fpout,"%18.15f ",yy[ny][i]);
	    /*printf("n i k    ny %d %d %d     %d\n",n,i,k,ny);*/
	 }
	 fprintf(fpout,"\n");
   }
   fclose(fpout);	
 }

 /*============== write OCN d13c.dat ============*/
if(NOCT >= 6){	
 /* d13c */
 fpout = fopen("d13c.dat","w");
 if(fpout == NULL) 
	ferrx("writedat(): Can't open file to write: d13c.dat");
 for(i=1;i<=lt;i++){ /* time */
	for(k=1;k<=NB;k++){ /* box */
	   tmp = (yy[k+5*NB][i]/yy[k][i])/RST;
	   fprintf(fpout,"%18.15f ",(tmp-1.)*1.e3);
	}
	fprintf(fpout,"\n");
 }
 fclose(fpout);			
}

 /*============== write OCN CO2 system.dat ============*/

 /* store temperature in matrix */
 for(i=1;i<=lt;i++){ /* time */
   for(k=1;k<=NB;k++){ /* box */
	   if(NOCT >= 4){
	     tmprm[k][i] = yy[k+3*NB][i];  /* temp as var */
	   } else {
	     tmprm[k][i] = tcb0[k];        /* temp static */
	   }
   }
 }
	
 for(k=1;k<=NB;k++)
     hgssvl[k] = HGUESS;
	   
 /* calculate CO2 paramater */
 for(i=1;i<=lt;i++){ /* time */
   for(k=1;k<=NB;k++){ /* box */
     fcsys(&tmp,&pco2ocn,&co3,&hpls,&ph,&tmp,&tmp,
            yy[k][i]*1.e-3,yy[k+1*NB][i]*1.e-3,hgssvl[k],
            tmprm[k][i],spm[k][1],spm[k][2],cacl,mgcl,s4cl);
	 kspfun(&kspc,&kspa,tmprm[k][i],spm[k][1],spm[k][2],cacl,mgcl,s4cl);
	 yyy[k][i][1] = pco2ocn;  /* uatm    */
	 yyy[k][i][2] = co3*1.e6; /* umol/kg */
	 yyy[k][i][3] = ph;
	 yyy[k][i][4] = co3*cacl/kspc;
	 yyy[k][i][5] = co3*cacl/kspa;
	 hgssvl[k] = hpls;
   }
 }

 /* write to file */
 for(n=1;n<=NCSWRT;n++){ /* output var */
   strcpy(fpstr,getco2str(n));
   fpout = fopen(fpstr,"w");
   if(fpout == NULL){
     sprintf(mssg,"writedat(): Can't open file to write: %s",fpstr);
     ferrx(mssg);
   }
   for(i=1;i<=lt;i++){ /* time */
     for(k=1;k<=NB;k++){ /* box */
       fprintf(fpout,"%18.15f ",yyy[k][i][n]);
     }
     fprintf(fpout,"\n");
   }
   fclose(fpout);	
 }

if(NCATM == 1){
 /* write atm pco2a.dat NOTE: Catm is scaled! */
 fpout = fopen("pco2a.dat","w");
 if(fpout == NULL) 
	ferrx("writedat(): Can't open file to write: pco2a.dat");
 for(i=1;i<=lt;i++)
    fprintf(fpout,"%18.15f \n",yy[NOCT*NB+1][i]/PPMTOMMSQ/CNTI);
 fclose(fpout);	
} /* NCATM */


if(NCCATM == 1){
 /* write atm d13ca.dat NOTE: 13Catm is NOT scaled! */
 fpout = fopen("d13ca.dat","w");
 if(fpout == NULL) 
	ferrx("writedat(): Can't open file to write: d13ca.dat");
 for(i=1;i<=lt;i++){
	tmp = (yy[NOCT*NB+2][i]/(yy[NOCT*NB+1][i]/CNTI))/RST;
    fprintf(fpout,"%18.15f \n",(tmp-1.)*1.e3);
 }
 fclose(fpout);	
} /* NCCATM */

#ifdef FSED


 /*=============== write SED fc[a,i,p,(t)].dat ============*/
 for(n=1;n<=NOC;n++){ /* output var */
   strcpy(fpstr,getfcstr(n));
   fpout = fopen(fpstr,"w");
   if(fpout == NULL){
     sprintf(mssg,"writedat(): Can't open file to write: %s",fpstr);
     ferrx(mssg);
   }
   for(i=1;i<=lt;i++){ /* time */
     for(k=1;k<=NSD;k++){ /* sed box */
		ny = NOATM + k + (n-1)*NSD;
	    fprintf(fpout,"%.15e ",yy[ny][i]);
	    /*printf("n i k    ny %d %d %d     %d\n",n,i,k,ny);*/
	 }
	 fprintf(fpout,"\n");
   }
   fclose(fpout);	
 }

 /*=============== write ccd[a,i,p,(t)].dat ================*/
 /* calculate CCD */
 fccd = FCCD; /* CaCO3 fraction at CCD */
 fccderr = 2; /* error flag finterp()  */

 /* extract fc[a,i,p,(t)][time] from yy */
 for(i=1;i<=lt;i++){ /* time */
     for(k=1;k<=NSD;k++){ /* sed box */
		ny = NOATM + k;
        fcam[k][i] = yy[ny][i];
		ny = NOATM + k + 1*NSD;
        fcim[k][i] = yy[ny][i];
		ny = NOATM + k + 2*NSD;
        fcpm[k][i] = yy[ny][i];
#ifdef FTYS
		ny = NOATM + k + 3*NSD;
        fctm[k][i] = yy[ny][i];
#endif		 
	 }
 }

 /*** find jccd's and zv[jccd] over time ***/
 for(i=1;i<=lt;i++){ /* time */
    for(k=1;k<=NSD;k++){ /* sed box */
		fca[k] = fcam[k][i]; /* fc at time step i */
		fci[k] = fcim[k][i]; 
		fcp[k] = fcpm[k][i]; 
#ifdef FTYS	
		fct[k] = fctm[k][i]; 
#endif		 
	}
	/* interpolate fc(dsv) on depth var zv and  */
    /* find fc value closest to fccd:           */
	/* 1. subtract fccd from fc and take abs    */
	/* 2. then just need to find minimum        */
	for(k=1;k<=nz;k++){
        ya[k] = finterp(zv[k],dsv,fca,NSD,fccderr);
        yi[k] = finterp(zv[k],dsv,fci,NSD,fccderr);
        yp[k] = finterp(zv[k],dsv,fcp,NSD,fccderr);
#ifdef FTYS	
        yt[k] = finterp(zv[k],dsv,fct,NSD,fccderr);
#endif		 
		/* 1. subtract fccd and take abs */
		ya[k] = fabs(ya[k]-fccd);
		yi[k] = fabs(yi[k]-fccd);
		yp[k] = fabs(yp[k]-fccd);
#ifdef FTYS	
		yt[k] = fabs(yt[k]-fccd);
#endif		 
	}
	/* 2. find minimum */
    jccda[i] = fminv(ya,nz,2);
    jccdi[i] = fminv(yi,nz,2); 	 
    jccdp[i] = fminv(yp,nz,2);
#ifdef FTYS	
    jccdt[i] = fminv(yt,nz,2);
#endif		 
 } /* end time */

 /* store jccd's in matrix */
 for(i=1;i<=lt;i++){
    jccdm[i][1] = jccda[i];
    jccdm[i][2] = jccdi[i];
    jccdm[i][3] = jccdp[i];
#ifdef FTYS	
    jccdm[i][4] = jccdt[i];
#endif		 
 }

 /* write to file */
 for(n=1;n<=NOC;n++){ /* output var */
   strcpy(fpstr,getfcstr(NOC+n));
   fpout = fopen(fpstr,"w");
   if(fpout == NULL){
     sprintf(mssg,"writedat(): Can't open file to write: %s",fpstr);
     ferrx(mssg);
   }
   for(i=1;i<=lt;i++) /* time */
	   fprintf(fpout,"%7.2f\n",zv[jccdm[i][n]]);	   
   fclose(fpout);
 }


 #ifdef FSEDCC 
 #ifdef WRTSEDCC

 /*=============== write SED f13c[a,i,p,(t)].dat ============*/
 for(n=1;n<=NOC;n++){ /* output var */
   strcpy(fpstr,getfcstr(2*NOC+n));
   fpout = fopen(fpstr,"w");
   if(fpout == NULL){
     sprintf(mssg,"writedat(): Can't open file to write: %s",fpstr);
     ferrx(mssg);
   }
   for(i=1;i<=lt;i++){ /* time */
     for(k=1;k<=NSD;k++){ /* sed box */
		ny = NOATM + (3+KTY)*NSD + k + (n-1)*NSD;
	    fprintf(fpout,"%.15e ",yy[ny][i]*CNTI); /* scaling *CNTI */
	    /* printf("n i k    ny %d %d %d     %d\n",n,i,k,ny);*/
	 }
	 fprintf(fpout,"\n");
   }
   fclose(fpout);	
 }
 #endif
 #endif

 free_dvector(fca,1,NSD);	
 free_dvector(fci,1,NSD);	
 free_dvector(fcp,1,NSD);	
 free_dmatrix(fcam,1,NSD,1,lt);	
 free_dmatrix(fcim,1,NSD,1,lt);	
 free_dmatrix(fcpm,1,NSD,1,lt);	
 free_dvector(ya,1,nz);	
 free_dvector(yi,1,nz);	
 free_dvector(yp,1,nz);	
 free_ivector(jccda,1,lt);	
 free_ivector(jccdi,1,lt);	
 free_ivector(jccdp,1,lt);	
 free_imatrix(jccdm,1,lt,1,NOC);
#ifdef FTYS	
 free_dvector(fct,1,NSD);	
 free_dmatrix(fctm,1,NSD,1,lt);	
 free_dvector(yt,1,nz);	
 free_ivector(jccdt,1,lt);	
#endif

#endif /* FSED */

 free_dmatrix(tmprm,1,NB,1,lt);	
 free_d3tensor(yyy,1,NB,1,lt,1,NCSWRT);	
 free_dvector(hgssvl,1,NB);	

}
/*============================================================*/
/*==================== writedat() END ========================*/
/*============================================================*/

/*============================================================*/
/*==================== getocnstr() ===========================*/
/*============================================================*/
char *getocnstr(int n)
{
 char *str;

switch(n){
	case 1:
		str = "dic.dat";
		break;
	case 2:
		str = "alk.dat";
		break;
	case 3:
		str = "po4.dat";
		break;
	case 4:
		str = "tcb.dat";
		break;
	case 5:
		str = "dox.dat";
		break;
	case 6:
		str = "dicc.dat";
		break;
	default:	
		str = "tmp.dat";
		break;
}
	
 return(str);
}
/*============================================================*/
/*==================== getocnstr() END =======================*/
/*============================================================*/


/*============================================================*/
/*==================== getco2str() ===========================*/
/*============================================================*/
char *getco2str(int n)
{
 char *str;

switch(n){
	case 1:
		str = "pco2ocn.dat";
		break;
	case 2:
		str = "co3.dat";
		break;
	case 3:
		str = "ph.dat";
		break;
	case 4:
		str = "omegaclc.dat";
		break;
	case 5:
		str = "omegaarg.dat";
		break;
	default:	
		str = "tmp.dat";
		break;
}
	
 return(str);
}
/*============================================================*/
/*==================== getco2str() END =======================*/
/*============================================================*/


/*============================================================*/
/*==================== getfcstr() ============================*/
/*============================================================*/
char *getfcstr(int n)
{
 char *str;

#ifdef FTYS	
switch(n){
	case 1:
		str = "fca.dat";
		break;
	case 2:
		str = "fci.dat";
		break;
	case 3:
		str = "fcp.dat";
		break;
	case 4:
		str = "fct.dat";
		break;
	case 5:
		str = "ccda.dat";
		break;
	case 6:
		str = "ccdi.dat";
		break;
	case 7:
		str = "ccdp.dat";
		break;
	case 8:
		str = "ccdt.dat";
		break;
	case 9:
		str = "fcca.dat";
		break;
	case 10:
		str = "fcci.dat";
		break;
	case 11:
		str = "fccp.dat";
		break;		
	case 12:
		str = "fcct.dat";
		break;		
	default:	
		str = "tmp.dat";
		break;
}
#else
switch(n){
	case 1:
		str = "fca.dat";
		break;
	case 2:
		str = "fci.dat";
		break;
	case 3:
		str = "fcp.dat";
		break;
	case 4:
		str = "ccda.dat";
		break;
	case 5:
		str = "ccdi.dat";
		break;
	case 6:
		str = "ccdp.dat";
		break;
	case 7:
		str = "fcca.dat";
		break;
	case 8:
		str = "fcci.dat";
		break;
	case 9:
		str = "fccp.dat";
		break;		
	default:	
		str = "tmp.dat";
		break;
}
#endif
	
 return(str);
}
/*============================================================*/
/*==================== getfcstr() END ========================*/
/*============================================================*/

