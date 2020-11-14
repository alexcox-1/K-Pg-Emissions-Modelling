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
#include <math.h>
#include "defs.h"
#include "utils.h"

/* Order of equilibrium constants */
#define NK1  1 /* CO2 1 */
#define NK2  2 /* CO2 2 */
#define NKB  3 /* Boron */
#define NKW  4 /* H2O   */
#define NKC  5 /* calc  */
#define NKA  6 /* arag  */
#define NTK (6) /* total*/

/*============================================================*/
/*==================== fcsys() ===============================*/
/*============================================================*/
/* input:   dic,alk,hguess,T,S,P,Ca,Mg,SO4. 
   returns: co2,pco2,co3,h,ph,kh,o2 

   Equilibrium constants: 
   Mehrbach et al. (1973) refit by Lueker et al. (2000).
   pH-scale: 'total'. mol/kg-soln
   For more info, see "CO2 in seawater: Equilibrium,
   Kinetics, Isotopes" (Zeebe and Wolf-Gladrow, 2001)

   updates:

   09/03/18 new K* corrections (ZeebeTyrrell18) KCORR18	   
   10/28/11 O2 solubility: ml/kg => ml/l	   
   10/19/11 falpha(...,int ntsl). ntsl=length kksv is now arg.
   10/08/11 add Ca,Mg to fcsys()
   09/18/11 Total boron: 416.0 DOE94 => 432.5 Lee10
   09/17/11 added falpha().
   09/11/11 include dissolved O2 (dox). o2 -> o2sat.
   04/19/11 added error handling if [H+] < 0. for high ALK/DIC
            Follow's CO2 routine gives negative [H+]  
   03/26/11 a2 for kb set to -2.608, not +2.608
   03/12/11 added kspfun and removed kspc and kspa from 
	        fcsys()
   02/28/11 new file
				
 */
void fcsys(double *co2, double *pco2, double *co3,double *h,
           double *ph,double *kh,double *o2sat,double dic,double alk,
           double hgss,double tc,double sal,double prs,double cacl,
           double mgcl,double s4cl)
{
 int i,nmab=5;
 double khx,k1,k2,kb,kw,o2x,bor,tk,**mab;
 /*double kspc,kspa;*/
 double hx,boh4g,fg,calkg,gam,tmp;

 mab = dmatrix(1,nmab,1,NTK);	
	 
 tk    = tc + TKLV;	
 bor   = (BORT*(sal/35.))*1.e-6; /* 416 DOE94, 432.5 Lee10 */ 

 hx    = hgss;

 /* get equilibrium constants and O2 */
 khx  = fkh(tk,sal);
 k1   = fk1(tk,sal);
 k2   = fk2(tk,sal);
 kb   = fkb(tk,sal);
 kw   = fkw(tk,sal);
 /*kspc = fkspc(tk,sal);
   kspa = fkspa(tk,sal);*/
 o2x  = fo2(tk,sal);	 

 /* pressure correction */
 if(prs > 0.0){
   getpcoeff(mab);
   k1 *= fpcorrk(tk,prs,mab,NK1);
   k2 *= fpcorrk(tk,prs,mab,NK2);
   kb *= fpcorrk(tk,prs,mab,NKB);
   kw *= fpcorrk(tk,prs,mab,NKW);
   /*kspc *= fpcorrk(tk,prs,mab,NKC);
     kspa *= fpcorrk(tk,prs,mab,NKA);*/
 }

#define UPRINTCS

#ifdef PRINTCS
 double kspc,kspa;
 kspfun(&kspc,&kspa,tc,sal,prs,CAM,MGM,S4M);	
 printf("\n%f %f %f rc sal prs\n",tc,sal,prs);	 
 printf("\n%20.15e  k1\n",k1);	 
 printf(  "%20.15e  k2\n",k2);	 
 printf("%20.15e  kspc\n",kspc);	 
 printf("%20.15e  kspa\n",kspa);
#endif	
	
 /* Ca,Mg,SO4 corrections K1, K2 */
#define KEPS (1.e-12)
 if(fabs((cacl-CAM)/CAM) > KEPS || fabs((mgcl-MGM)/MGM) > KEPS ||
    fabs((s4cl-S4M)/S4M) > KEPS){
   k1 += fdelk1(k1,cacl,mgcl,s4cl);
   k2 += fdelk2(k2,cacl,mgcl,s4cl);
 }
		
 /* iterative solution for H+, Follows et al. 2006 */ 	 
 for(i=1;i<=HIMAX;i++){
   hgss  = hx;
   boh4g = bor*kb/(hgss + kb);
   fg    = -boh4g-(kw/hgss)+hgss;
   calkg = alk + fg;
   gam   = dic/calkg;
   tmp   = (1.-gam)*(1.-gam)*k1*k1-4.*k1*k2*(1.-2.*gam);
   hx    = 0.5*((gam-1.)*k1 + sqrt(tmp));	
   if((fabs(hx-hgss)) <= hgss*HCONV){
#ifdef PRINTCS	 
     printf("\n %d i\n",i);
#endif	   
     break;
   } 
   if(i == HIMAX){ 
	 fprintf(stderr,"\n[H+] = %e",hx);  
     fprintf(stderr,"\n[H+] iteration did not converge after %d steps",HIMAX);
     if(hx < 0.0)
	   fprintf(stderr,"\n[H+] < 0. Your ALK/DIC ratio may be too high");  
	 ferrwrt("csys(): check CO2 system input");  
   }
 }

 *co2   =  dic/(1.+k1/hx+k1*k2/hx/hx);
 *pco2  = *co2*1.e6/khx;
 *co3   =  dic/(1+hx/k2+hx*hx/k1/k2);
 *h     =  hx;
 *ph    = -log10(hx);
 *kh    =  khx;
 *o2sat =  o2x;

#ifdef PRINTCS
 /*
 int n;
 for(i=1;i<=nmab;i++){
   for(n=1;n<=NTK;n++)
	 printf("%15.5e  mab[%d][%d]\n",mab[i][n],i,n);	 
 } 
 */
 printf("%20.15e  kh\n",khx);	 
 printf("%20.15e  k1\n",k1);	 
 printf("%20.15e  k2\n",k2);	 
 printf("%20.15e  kb\n",kb);	 
 printf("%20.15e  kw\n",kw);
 kspfun(&kspc,&kspa,tc,sal,prs,cacl,mgcl,s4cl);	
 printf("%20.15e  kspc\n",kspc);	 
 printf("%20.15e  kspa\n",kspa);
 printf("%20.15e  O2  \n",o2x);	 
	 exit(0);
#endif

 free_dmatrix(mab,1,nmab,1,NTK);	
	 
}
/*============================================================*/
/*==================== fcsys() END ===========================*/
/*============================================================*/

/*============================================================*/
/*==================== kspfun() ==============================*/
/*============================================================*/
void kspfun(double *kspc,double *kspa,double tc,double sal,
            double prs,double cacl,double mgcl,double s4cl)
{
 int nmab=5;
 double tk,**mab;	 

 mab = dmatrix(1,nmab,1,NTK);	
	 
 tk    = tc + TKLV;	

 /* get equilibrium constants */
 *kspc = fkspc(tk,sal);
 *kspa = fkspa(tk,sal);

 /* pressure correction */
 if(prs > 0.0){
   getpcoeff(mab);
   *kspc *= fpcorrk(tk,prs,mab,NKC);
   *kspa *= fpcorrk(tk,prs,mab,NKA);
 }

 /* Ca,Mg,SO4 corrections Ksp */
 if(fabs((cacl-CAM)/CAM) > KEPS || fabs((mgcl-MGM)/MGM) > KEPS ||
    fabs((s4cl-S4M)/S4M) > KEPS){
#ifdef KCORR18 
   *kspc += fdelksp(*kspc,cacl,mgcl,s4cl);
   *kspa += fdelksp(*kspa,cacl,mgcl,s4cl);
#else	 
   *kspc *= (1.-ALPKC*(MGM/CAM-mgcl/cacl));
#endif	 
 }
   
 free_dmatrix(mab,1,nmab,1,NTK);
	
}	 
/*============================================================*/
/*==================== kspfun() END ==========================*/
/*============================================================*/


/*============================================================*/
/*==================== fkh() =================================*/
/*============================================================*/
double fkh(double tk, double sal)
{
/* ===========================================================

  kh (K Henry)
  CO2(g) <-> CO2(aq.)
  kh 	= [CO2]/ p CO2  
  Weiss (1974)   [mol/kg/atm]
  =========================================================== */
 double kh,tmp;
	
 tmp  = 9345.17/tk - 60.2409 + 23.3585*log(tk/100.);
 tmp += sal*(0.023517-0.00023656*tk+0.0047036e-4*tk*tk);

 kh = exp(tmp);
 return(kh);
	
}   
/*============================================================*/
/*==================== fkh() END =============================*/
/*============================================================*/


/*============================================================*/
/*==================== fk1() =================================*/
/*============================================================*/
double fk1(double tk, double sal)
{
/* ===========================================================
    first acidity constant:
    [H^+] [HCO_3^-] / [H_2CO_3] = K_1
 
    Mehrbach et al. (1973) refit by Lueker et al. (2000).
 
    pH-scale: 'total'. mol/kg-soln     
   =========================================================== */
 double k1,tmp;
	
 tmp  = 3633.86/tk - 61.2172 + 9.6777*log(tk) - 0.011555*sal;
 tmp += 0.0001152*sal*sal;
 k1   = pow(10.0,-tmp);
       
 return(k1);

}  
/*============================================================*/
/*==================== fk1() END =============================*/
/*============================================================*/

/*============================================================*/
/*==================== fk2() =================================*/
/*============================================================*/
double fk2(double tk, double sal)
{    
/* ===========================================================
 
    second acidity constant:
    [H^+] [CO_3^--] / [HCO_3^-] = K_2
 
    Mehrbach et al. (1973) refit by Lueker et al. (2000).
 
    pH-scale: 'total'. mol/kg-soln
   =========================================================== */
 double k2,tmp;

 tmp  = 471.78/tk + 25.9290 - 3.16967*log(tk) - 0.01781*sal;
 tmp += 0.0001122*sal*sal;
 k2   = pow(10.0,-tmp);
	
 return(k2);

}
/*============================================================*/
/*==================== fk2() END =============================*/
/*============================================================*/

/*============================================================*/
/*==================== fkb() =================================*/
/*============================================================*/
double fkb(double tk, double sal)
{
/* ===========================================================
   Boric acid constant: 

   Kbor = [H+][B(OH)4-]/[B(OH)3] = kp7 / km7

   (Dickson, 1990 in Dickson and Goyet, 1994, Chapter 5, p. 14)
   pH-scale: 'total'
   =========================================================== */
 double kb,tmp1,tmp2,tmp3,tmp4;
	
   tmp1  = -8966.90-2890.53*sqrt(sal)-77.942*sal;
   tmp1 +=  1.728*pow(sal,1.5)-0.0996*sal*sal;
   tmp2  =  148.0248+137.1942*sqrt(sal)+1.62142*sal;
   tmp3  = -24.4344-25.085*sqrt(sal)-0.2474*sal;
   tmp3 *= log(tk);

   tmp4  = tmp1/tk + tmp2 + tmp3 + 0.053105*sqrt(sal)*tk; 

   kb = exp(tmp4);

   return(kb);	

}   
/*============================================================*/
/*==================== fkb() END =============================*/
/*============================================================*/


/*============================================================*/
/*==================== fkw() =================================*/
/*============================================================*/
double fkw(double tk, double sal)
{
/* ==============================================================

   Kwater = [H+] [OH-] : units: (mol/kg)**2
   ion product of water as a function of salinity (s [psu]) and 
   absolute temperature (t [K])
  
	Millero (1995)(in Dickson and Goyet (1994, Chapter 5, p.18))
	$K_w$ in mol/kg-soln.
	pH-scale: pH$_{Hansson}$ ('total` scale).

   ============================================================== */
 double kw,tmp1,tmp2;

   tmp1  = -13847.26/tk + 148.96502 - 23.6521*log(tk);
   tmp2  =  118.67/tk - 5.977 + 1.0495*log(tk);
   tmp2 *=  sqrt(sal);
   tmp1 +=  tmp2 - 0.01615*sal;
	
   kw  = exp(tmp1); 
               
   return(kw);

} 
/*============================================================*/
/*==================== fkw() END =============================*/
/*============================================================*/


/*============================================================*/
/*==================== fkspc() ===============================*/
/*============================================================*/
double fkspc(double tk, double sal)
{
/* ===========================================================
   apparent solubility product of calcite
 
   Kspc = [Ca2+]T [CO32-]T
 
   where $[]_T$ refers to the equilibrium total 
  (free + complexed) ion concentration.
 
   Mucci 1983 (mol/kg-soln)^2
   =========================================================== */
 double kspc,tmp;

   tmp  = -171.9065-0.077993*tk+2839.319/tk+71.595*log10(tk);
   tmp += (-0.77712+0.0028426*tk+178.34/tk)*sqrt(sal);
   tmp += -0.07711*sal+0.0041249*pow(sal,1.5);

   kspc = pow(10.,tmp);

   return(kspc);  
       
}   
/*============================================================*/
/*==================== fkspc() END ===========================*/
/*============================================================*/

/*============================================================*/
/*==================== fkspa() ===============================*/
/*============================================================*/
double fkspa(double tk, double sal)
{
/* ===========================================================
   apparent solubility product of aragonite
 
   Kspa = [Ca2+]T [CO32-]T
 
   where $[]_T$ refers to the equilibrium total 
  (free + complexed) ion concentration.
 
   Mucci 1983 (mol/kg-soln)^2
   =========================================================== */
 double kspa,tmp;

   tmp  = -171.945-0.077993*tk+2903.293/tk+71.595*log10(tk);
   tmp += (-0.068393+0.0017276*tk+88.135/tk)*sqrt(sal);
   tmp += -0.10018*sal+0.0059415*pow(sal,1.5);

   kspa = pow(10.,tmp);

   return(kspa);  
       
}      
/*============================================================*/
/*==================== fkspa() END ===========================*/
/*============================================================*/

/*============================================================*/
/*==================== fo2() =================================*/
/*============================================================*/
double fo2(double tk, double sal)
{
/* ===========================================================

   solubility of O2 
   Weiss (1970) DSR, 17, p. 721.

  =========================================================== */
 double o2,tmp1,tmp2,tk1;

 tk1   = (tk/100.);

 /* ml (STP) / kg	
 tmp1  = -177.7888+255.5907*100./tk+146.4813*log(tk1);
 tmp1 += -22.2040*(tk/100.);
 tmp2  = sal*(-0.037362+0.016504*tk1-0.0020564*tk1*tk1);
 */

 /* ml (STP) / l = (l/m3) */
 tmp1  = -173.4292+249.6339*100./tk+143.3483*log(tk1);
 tmp1 += -21.8492*(tk/100.);
 tmp2  = sal*(-0.033096+0.014259*tk1-0.0017000*tk1*tk1);
	
 o2  = exp(tmp1+tmp2)/22.4; /* -> mol/m3 */

 return(o2);
	
}   
/*============================================================*/
/*==================== fo2() END =============================*/
/*============================================================*/


/*============================================================*/
/*==================== getpcoeff() ===========================*/
/*============================================================*/
void getpcoeff(double **mab)
{
 int n;
               /*  k1      k2      kb      kw       kspc   kspa */
 double	
  a0[NTK+1] = {0,-25.50 ,-15.82 ,-29.48 ,-25.60  ,-48.76,-46.00  },
  a1[NTK+1] = {0, 0.1271,-0.0219, 0.1622, 0.2324 , 0.5304, 0.5304},
  a2[NTK+1] = {0, 0.0   , 0.0   ,-2.608 ,-3.6246 , 0.0   , 0.0   },
  b0[NTK+1] = {0,-3.08  , 1.13  ,-2.84  ,-5.13   ,-11.76 ,-11.76 },                    
  b1[NTK+1] = {0, 0.0877,-0.1475, 0.0   , 0.0794 , 0.3692, 0.3692};
	  
 /* NOTE: a2 for kb should be -2.608, not +2.608 (thanks to
    James Rae). use plus when checking old results             */

 
 for(n=1;n<=NTK;n++){
  mab[1][n] = a0[n];
  mab[2][n] = a1[n];
  mab[3][n] = a2[n]*1.e-3;
  mab[4][n] = b0[n]*1.e-3;
  mab[5][n] = b1[n]*1.e-3;
 }
	
}
/*============================================================*/
/*==================== getpcoeff() END =======================*/
/*============================================================*/


/*============================================================*/
/*==================== fpcorrk() =============================*/
/*============================================================*/
/* pressure correction for equilibrium constants              */
/* index: 1=k1, 2=k2, 3=kb, 4=kw, 5=kspc, 6=kspa              */

#define RGP (83.131)

double fpcorrk(double tk,double prs,double **mab, int n)
{
 double pcorrk,tc,deltav,deltak,lnkpok0;

/*	order of coefficients
  mab[1][n] = a0[n];
  mab[2][n] = a1[n];
  mab[3][n] = a2[n];
  mab[4][n] = b0[n];
  mab[5][n] = b1[n];	
*/
	
 tc = tk - TKLV;
 deltav   = mab[1][n] + mab[2][n]*tc + mab[3][n]*tc*tc;
 deltak   = mab[4][n] + mab[5][n]*tc;  
 lnkpok0  = -(deltav/(RGP*tk))*prs; 
 lnkpok0 += (0.5*deltak/(RGP*tk))*prs*prs;
 
 pcorrk = exp(lnkpok0); 

 return(pcorrk); 
	
} 
/*============================================================*/
/*==================== fpcorrk END ===========================*/
/*============================================================*/



/*============================================================*/
/*==================== fdelk1() ==============================*/
/*============================================================*/
double fdelk1(double k1,double cacl,double mgcl,double s4cl)
{
 double delk1,sk1ca,sk1mg,sk1s4,delk1ca,delk1mg,delk1s4;

 /* sensitivity parameters for Ca,Mg,SO4 effect on K* */
 /* new K* corrections (ZeebeTyrrell18)               */
#ifdef KCORR18	
 sk1ca =   5.e-3;
 sk1mg =  17.e-3;
 sk1s4 = 208.e-3;
#else	
 sk1ca =  33.73e-3;
 sk1mg = 155.00e-3;
 sk1s4 =   0.0;
#endif	
	
 /* add Ca,Mg,SO4 correction K* (Ben-Yaakov & Goldhaber, 1973) */
 delk1ca = sk1ca*k1*(cacl/CAM-1.);
 delk1mg = sk1mg*k1*(mgcl/MGM-1.);
 delk1s4 = sk1s4*k1*(s4cl/S4M-1.);
	
 delk1   = delk1ca+delk1mg+delk1s4;
		
 return(delk1);
}
/*============================================================*/
/*==================== fdelk1() END ==========================*/
/*============================================================*/



/*============================================================*/
/*==================== fdelk2() ==============================*/
/*============================================================*/
double fdelk2(double k2,double cacl,double mgcl,double s4cl)
{
 double delk2,sk2ca,sk2mg,sk2s4,delk2ca,delk2mg,delk2s4;

 /* sensitivity parameters for Ca,Mg,SO4 effect on K* */
 /* new K* corrections (ZeebeTyrrell18)               */
#ifdef KCORR18		
 sk2ca = 157.e-3;
 sk2mg = 420.e-3;
 sk2s4 = 176.e-3;
#else	
 sk2ca =  38.85e-3;
 sk2mg = 442.00e-3;
 sk2s4 =   0.00;
#endif
	
 /* add Ca,Mg,SO4 correction K* (Ben-Yaakov & Goldhaber, 1973) */
 delk2ca = sk2ca*k2*(cacl/CAM-1.);
 delk2mg = sk2mg*k2*(mgcl/MGM-1.);
 delk2s4 = sk2s4*k2*(s4cl/S4M-1.);
	
 delk2   = delk2ca+delk2mg+delk2s4;
		
 return(delk2);
}
/*============================================================*/
/*==================== fdelk2() END ==========================*/
/*============================================================*/



/*============================================================*/
/*==================== fdelksp() =============================*/
/*============================================================*/
double fdelksp(double ksp,double cacl,double mgcl,double s4cl)
{
 double delksp,skspca,skspmg,sksps4,delkspca,delkspmg,delksps4;

 /* sensitivity parameters for Ca,Mg,SO4 effect on K* */
 /* new K* corrections (ZeebeTyrrell18)               */
 skspca = 185.e-3;
 skspmg = 518.e-3;
 sksps4 = 106.e-3;
	
 /* add Ca,Mg,SO4 correction K* (Ben-Yaakov & Goldhaber, 1973) */
 delkspca = skspca*ksp*(cacl/CAM-1.);
 delkspmg = skspmg*ksp*(mgcl/MGM-1.);
 delksps4 = sksps4*ksp*(s4cl/S4M-1.);
	
 delksp   = delkspca+delkspmg+delksps4;
		
 return(delksp);
}
/*============================================================*/
/*==================== fdelksp() END =========================*/
/*============================================================*/



/*============================================================*/
/*==================== falpha() ==============================*/
/*============================================================*/
void falpha(double *alpdb,double *alpdg,double *alpcb, 
            double *alpu,double *tcb, int *kksv, int ntsl)
{
/* ===========================================================
  13C alphas for CO2 gas exchange.
  Mook 1986 or Zhang et al. 1995
  =========================================================== */
 int k,kk;
 double tkb,epsdb,epsdg,epscb;
#ifndef MOOK 	
 double epsbg,epscg; /* only used for Zhang */
#endif		
	
 for(k=1;k<=ntsl;k++){
   kk    = kksv[k];
   tkb   = tcb[kk] + TKLV;
#ifdef MOOK	
   epsdb = -9866./tkb + 24.12;
   epsdg =  -373./tkb +  0.19;
   epscb =  -867./tkb +  2.52;
#else /* Zhang */
   epsbg =  -0.114*tcb[kk] + 10.78;
   epsdg =  0.0049*tcb[kk] -  1.31;
   epscg =  -0.052*tcb[kk] +  7.22;
   epsdb = (epsdg-epsbg)/(1.+epsbg/1.e3);
   epscb = (epscg-epsbg)/(1.+epsbg/1.e3);	 
#endif	 
   alpdb[kk] = epsdb/1.e3+1.;
   alpdg[kk] = epsdg/1.e3+1.;
   alpcb[kk] = epscb/1.e3+1.;
 }

#ifdef MOOK	
 *alpu = 0.9995;
#else /* Zhang */
 *alpu = 0.9991;
#endif		
		
}   
/*============================================================*/
/*==================== falpha() END ==========================*/
/*============================================================*/




