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

/*============================================================*/
/*==================== thmfun() ==============================*/
/*============================================================*/
/* returns yp for given conveyor configuration
   and mixing. same for all ocean tracers                     

 updates: 

 10/21/11 include Tethys. thl>thcl
 03/10/11 thc>thl, mv>mvl, mhd>mhdl (local)
 02/18/11 new file

 Boxes: Atlantic, Indic, Pacific, High, Tethys

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

void thmfun(double *y,double *yp,int fconvl,double thcl,
            double tsol,double ttol,double *mvl,
            double *mhdl,double *vb,double ga,double ta,
            double gi,double ti)
{
 int k;
 double fa=0.4,fi=0.3,fp=0.3;
#ifndef FTYS
 ttol = tsol; /* not used */
 tsol = ttol; /* not used */
#endif
	
 for(k=1;k<=NB;k++)
    yp[k] = 0.0; /* 0.0 -y[k] */

 if(fconvl == 1){
  /* circulation: NADW */
  yp[4] = ( ga*thcl*y[5]+ta*thcl*y[7]-   thcl*y[4] )/vb[4]; /* IA */ 
  yp[5] = ( gi*thcl*y[6]+ti*thcl*y[8]-ga*thcl*y[5] )/vb[5]; /* II */ 
  yp[6] =                      gi*thcl*(y[9] -y[6] )/vb[6]; /* IP */
  yp[7] =                         thcl*(y[10]-y[7] )/vb[7]; /* DA */
  yp[8] =                      ga*thcl*(y[7] -y[8] )/vb[8]; /* DI */
  yp[9] =                      gi*thcl*(y[8] -y[9] )/vb[9]; /* DP */
  yp[10]=                         thcl*(y[4]-y[10])/vb[10]; /* H  */
 } 
 if(fconvl == 2){
  /* circulation: NPDW */
  yp[6] =     ( ga*thcl*  y[5]+
                ta*thcl*  y[9]-
                   thcl*  y[6]      ) /vb[6];       /* IP */
  yp[5] =     ( gi*thcl*  y[4]+
                ti*thcl*  y[8]-
                ga*thcl*  y[5]      ) /vb[5];       /* II */ 
  yp[4] =       gi*thcl*( y[7] -y[4]) /vb[4];       /* IA */
  yp[9] =          thcl*(y[10] -y[9]) /vb[9];       /* DP */
  yp[8] =       ga*thcl*( y[9] -y[8]) /vb[8];       /* DI */
  yp[7] =       gi*thcl*( y[8] -y[7]) /vb[7];       /* DA */
  yp[10] =         thcl*( y[6]-y[10])/vb[10];       /* H  */
  /* add: SO contribution */
  yp[4] +=       fa*tsol*( y[7] -y[4]) /vb[4];      /* IA */
  yp[5] +=       fi*tsol*( y[8] -y[5]) /vb[5];      /* II */
  yp[6] +=       fp*tsol*( y[9] -y[6]) /vb[6];      /* IP */
  yp[7] +=       fa*tsol*(y[10] -y[7]) /vb[7];      /* DA */
  yp[8] +=       fi*tsol*(y[10] -y[8]) /vb[8];      /* DI */
  yp[9] +=       fp*tsol*(y[10] -y[9]) /vb[9];      /* DP */
  yp[10]+=      (fa*tsol*( y[4]-y[10])
                +fi*tsol*( y[5]-y[10])
                +fp*tsol*( y[6]-y[10]))/vb[10];     /* H  */
 }
 if(fconvl == 3){
  /* circulation: SO */
  yp[4]  =       fa*thcl*( y[7] -y[4]) /vb[4];      /* IA */
  yp[5]  =       fi*thcl*( y[8] -y[5]) /vb[5];      /* II */
  yp[6]  =       fp*thcl*( y[9] -y[6]) /vb[6];      /* IP */
  yp[7]  =       fa*thcl*(y[10] -y[7]) /vb[7];      /* DA */
  yp[8]  =       fi*thcl*(y[10] -y[8]) /vb[8];      /* DI */
  yp[9]  =       fp*thcl*(y[10] -y[9]) /vb[9];      /* DP */
  yp[10] =      (fa*thcl*( y[4]-y[10])
                +fi*thcl*( y[5]-y[10])
                +fp*thcl*( y[6]-y[10]))/vb[10];     /* H  */
 } /* fconvl */

#ifdef FTYS /* TETHYS */
 /* circulation: Teth LI > LT > DT > DI > II > LI ... etc */
 yp[11] =         ttol*( y[2]-y[11])/vb[11];        /* LT */
 yp[13] =         ttol*(y[11]-y[13])/vb[13];        /* DT */
 yp[8]  = yp[8]+  ttol*(y[13] -y[8]) /vb[8];        /* DI */
 yp[5]  = yp[5]+  ttol*( y[8] -y[5]) /vb[5];        /* II */
 yp[2]  = yp[2]+  ttol*( y[5] -y[2]) /vb[2];        /* LI */
#endif  
	
 /* mixing AIP H */ 
 for(k=1;k<=3;k++){
    yp[k]   +=  mvl[k]*(y[k+3]-y[k]  )/vb[k]  ; /* L-I */
    yp[k+3] +=  mvl[k]*(y[k]  -y[k+3])/vb[k+3]; /* I-L */
    yp[k+6] += mhdl[k]*(y[10] -y[k+6])/vb[k+6]; /* D-H */
    yp[10]  += mhdl[k]*(y[k+6]-y[10] )/vb[10] ; /* H-D */
 }

#ifdef FTYS /* TETHYS */
 /* mixing Tethys */ 
 yp[12] =          mvl[4]*(y[11]-y[12])/vb[12]; /* IT-LT */
 yp[11] = yp[11]+  mvl[4]*(y[12]-y[11])/vb[11]; /* LT-IT */
 yp[12] = yp[12]+  mvl[5]*( y[5]-y[12])/vb[12]; /* IT-II */
  yp[5] =  yp[5]+  mvl[5]*(y[12] -y[5]) /vb[5]; /* II-IT */
 yp[11] = yp[11]+ mhdl[4]*(y[13]-y[11])/vb[11]; /* LT-DT */
 yp[13] = yp[13]+ mhdl[4]*(y[11]-y[13])/vb[13]; /* DT-LT */
#endif 
}
/*============================================================*/
/*==================== thmfun() END ==========================*/
/*============================================================*/
