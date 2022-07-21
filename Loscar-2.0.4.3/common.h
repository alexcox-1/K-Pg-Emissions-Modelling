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

/*============================================================*/
/*==================== common.h ==============================*/
/*============================================================*/


/*============================================================*/
/*==================== global variables ======================*/
/*============================================================*/
/* Agreed, this is bad practice. My excuse: (1) The C program 
   was developed based on a Matlab script. (2) Use of function 
   derivs() is impractical without globals.                   

   common.h is included in 
	   
   emiss.c
   initfree.c
   initstart.c
   loscarderivs.c
   readparms.c
   loscar.c

   updates:

   10/21/11 include Tethys variables
   04/10/11 checked and cleaned up global vars; a few renamed:
            tv>tmv, mv>mxv, ta>thbra, ti>thbri, ei>frei,
	        hgss>hgssv
   04/07/11 always declare *dsv etc. for writedat()
	        ldrestart,svrestart -> int.
   02/19/11 new file    
 */

/* solver control */
int kmax,kount;         /* max # output values,        see solver */
double epslvr;          /* solver accuracy                        */
double *tmv,**yy,dxsav, /* output vars, save interval, see solver */
       h1slv,hminslv;   /* 1st guess and min step size            */

/* model parameters/vars */	
int *kkv,*kiv,ldrestart,svrestart,ffflag,tsnsflag,ltem,cntrfflag,
    fconv,cinpflag,sflag,ltSem,ltExp,expflag,reminflag,ltRemin;
double *vb,*ab,*hb,thc0,thc,t0,tfinal,*ystart,*mxv,*mhd,thbra,thbri,
       fepl,rrain,nuwd,frei,eph,*gp,*hgssv,*salv,*prsv,**spm,cac,mgc,s4c,
       fvc0,finc0,fkrg,pcsi,ncc,nsi,fcsml,*tems,*yems,*tcb0,*tcb0c,
       *fdapi,rincc,rvccc,rkrgcc,tso,tto,epscorg,cinp,dccinp,rccinp,
       tcin0,tcspan,sclim,rksp,*tSems,*ySems,*tExp,*yExp,*tRemin,*yRemin;
/* file names */
char fpldstr[BUFSIZ],fpsvstr[BUFSIZ],ffldstr[BUFSIZ],sldstr[BUFSIZ],expldstr[BUFSIZ],reminldstr[BUFSIZ];

#ifdef FSED
int *klid,*nlid;
double ncsd,kssd,*asva,*asvi,*asvp,phi0,phi1,hsl,fsh,
       *fc0a,*fc0i,*fc0p,frrf,*phiia,*phiii,*phiip;
 #ifdef FTYS
 int *klidt;
 double nsht,*asvt,*fc0t,*phiit;
 #endif
 #ifdef FSEDCC 
 double *fcc0a,*fcc0i,*fcc0p;
  #ifdef FTYS
   double *fcc0t;
  #endif /* FTYS   */
 #endif /* FSEDCC */
#endif /* FSED   */
/* always needed (for writedat()) */
int nz;
double *dsv,*zv;

