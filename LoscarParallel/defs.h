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

/*============================================================*/
/*========================= defs.h ===========================*/
/*============================================================*/
/*
 updates: 

  02/16/11 new file

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
#define VLOSCAR "2.0.4.3" /* Version. Check also: loscar.c */
#define VMONTH     (9)    /* V month */
#define VYEAR   (2018)    /* V year  */

#define XDBG (88.88)
#define CNTI (0.01)
#define LOG2 (0.693147180559945)
#define DEPS (1.e-9)
/* DBL_EPSILON= 2.2204460492503130808472633361816406250000e-16 
   (3.0.0-13-generic-pae #22-Ubuntu, gcc 4.6.1)  
*/

/* solver control etc.*/
#define CSCAL (1.0)   /* error scaling, set !?! */
#define MAXSTP 100000 /* max # of steps         */

/* oceans, vars, boxes, equations */
#ifdef FTYS /* Paleo  */
 /* NOTE: for Paleo version, define FTYS 
    during compile: make loscar PALEO=1 
 */
 #define NOC    (4)    /* # oceans 4, AIPT          */
 #define NB    (13)    /* # of ocean boxes 13       */
 #define NLS    (4)    /* # low-lat surf boxes      */
#define KTY     (1)    /* shift 13C sed index       */
#else       /* Modern */
 #define NOC    (3)    /* # oceans 3, AIP           */
 #define NB    (10)    /* # of ocean boxes 10       */
 #define NLS    (3)    /* # low-lat surf boxes      */
 #define KTY    (0)    /* shift 13C sed index       */
#endif
#define KOC    (NOC+KTY)
#define NOCT   (6)    /* # ocn tracers: dic,alk,.. */
#define NCATM  (1)    /*    C atmosphere 1         */
#define NCCATM (1)    /*   CC atmosphere 1         */
#define NHS    (1)    /* # high-lat surf boxes     */
#define NTS    (NLS+NHS)  /* # total surf boxes    */
#define NOATM  (NOCT*NB+NCATM+NCCATM)


/* sediment on: FSED, off: FSEDU    */
#define FSED
#ifdef FSED
 #define NSD (13)   /* # sed boxes per OCN   */
#else
 #define NSD  (0) 
#endif
#define FCCD (0.10) /* CaCO3 fraction at CCD */

/* 13C */
#define FSEDCC /* on: 12C+13C seds, off: 12C seds   */
#ifdef FSEDCC
 #define NSDCC (NSD)
#else
 #define NSDCC  (0) 
#endif
#define WRTSEDCCU /* write 13C sed to file on/off */

/* # of DEQs */
#define NEQ  (NOCT*NB+NCATM+NCCATM+NOC*NSD+NOC*NSDCC)

/* ocean etc. */
#define VOC (1.2918235e18)    /* (m3) volume ocean   */
#define AOC      (3.49e14)    /* (m2) area ocean     */
#define HAV      (VOC/AOC)    /* (m)  average depth  */
#define RHO      (1.025e3)    /* kg/m3 seawater      */
#define RHOS      (2.50e3)    /* kg/m3 sediment      */
#define SALOCN     (34.72)    /* ocean salinity      */

/* conversions */
#define YTOSEC    (3600.*24.*365.) /* years to secs     */
#define MMTOM     (1.e-3)          /* mmol  to mol      */
#define PPMTOMMSQ (2.2e15/12./AOC) /* ppmv  to mol/m2   */
#define MTOKGCACR      (100./1.e3) /* mol C to kg CaCO3 */

/* biological pump */
#define REDPC    (1./130.) /* 130     Redfield  P:C     */
#define REDNC   (15./130.) /* 15/130  Redfield  N:C     */
#define REDO2C (165./130.) /* 165/130 Redfield O2:C 169 */

/* CO2 chemistry */
#define TKLV  (273.15) /* TC <=> TKelvin          */
#define HGUESS (1.e-8) /* H+ first guess   1.e-8  */
#define HIMAX     (50) /* H+ max iterations   50  */
#define HCONV  (1.e-4) /* H+ rel. accuracy 1.e-4  */
#define NCSWRT     (5) /* # calc csys output vars */
#define HALK2DIC (1.6) /* High ALK/DIC ratio      */
#define BORT   (432.5) /* 416. DOE94, 432.5 Lee10 */
/* Ca,Mg,SO4 */
#define KCORR18 /* new K* correct (ZeebeTyrrell18) */ 	
#define CAM (10.3e-3)  /* (mol/kg) modern Ca  10.3 */
#define MGM (53.0e-3)  /* (mol/kg) modern Mg  53.0 */
#define S4M (28.2e-3)  /* (mol/kg) modern SO4 28.2 */
#define ALPKC (0.0833) /* slope Mg/Ca on Kspc      */

/* fossil fuel scenarios:                      */
/* max emmiss values to read from file 100,000 */
#define NEMSMAX (100000)
/* OCN Temp change (x deg C per 2xCO2)    */
/* Temp relax times (y) Surf, Intrm, Deep */
#define TAUS   (20.)
#define TAUI  (200.)
#define TAUD (1000.)
/* scale temperature to oder 1 (x1/TSCAL), see derivs() */
#define TSCAL  (10.)

/* MM kinetics HighLat PO4 uptake */
#define MMPO4H
#define PMMK (0.10e-3) /* mol/m3 (0.1 umol/l) */
#define PMM0 (0.10e-3)
#define PMMV (1.+PMM0/PMMK)

/* MM kinetics dissolved oxygen       */
#define KMMOX (0.001) /* mol/m3 0.001 */

/* Carbon-13 (dicc, ccatm) */
#define RST (0.011) /* 0.011 R standard (value irrelevant)   */
#define MOOK   /* 13alpha choice. default: Mook. else: Zhang */

/* comparisons, see dfind() */
#define SMALLER 1 /* <  */
#define SMALLEQ 2 /* <= */
#define EQUAL   3 /* == */
#define LARGEQ  4 /* >= */
#define LARGER  5 /* >  */
/* order of magnitude for comparison */
#define MZBH  (1.e2)   /* box height */
#define MZOM   (1.0)   /* omega      */
#define MZER (1.e-4)   /* erosion    */

/*============================================================*/
/*================== function declarations ===================*/
/*============================================================*/

/* csys.c */
void fcsys(double *co2, double *pco2, double *co3,double *h,
           double *ph,double *kh,double *o2sat,double dic,double alk,
           double hgss,double tc,double sal,double prs,double cacl,
           double mgcl,double s4cl);
void kspfun(double *kspc,double *kspa,double tc,double sal,
            double prs,double cacl,double mgcl,double s4cl);
double fkh(double tk, double sal);
double fk1(double tk, double sal);
double fk2(double tk, double sal);
double fkb(double tk, double sal);
double fkw(double tk, double sal);
double fkspc(double tk, double sal);
double fkspa(double tk, double sal);
double fo2(double tk, double sal);
void getpcoeff(double **mab);
double fpcorrk(double tk,double prs,double **mab, int n);
double fdelk1(double k1,double cacl,double mgcl,double s4cl);
double fdelk2(double k2,double cacl,double mgcl,double s4cl);
double fdelksp(double ksp,double cacl,double mgcl,double s4cl);
void falpha(double *alpdb,double *alpdg,double *alpcb, 
            double *alpu,double *tcb, int *kksv, int ntsl);
/* emiss.c */
void reademiss();
void readSemiss();
void readExp();
/* initfree.c */
void initfree(int flag);
/* initstart.c */
void initstart();
/* loscarderivs.c */
void derivs(double t,double *y, double *yp);
void jacobn(double x,double *y,double *dfdx,double **dfdy,
            int n);
void ferrwrt(char error_text[]);
/* matrix.c */
void ludcmp(double **a,int n,int *indx,double *d);
void lubksb(double **a,int n,int *indx,double *b); 
/* readparms.c */
void readparms(char *fparmstr);
void getparms(FILE *fparm,char *target,int *ffound,int *np,
              double *pval,char *pstr);
void setpstr(int *actflag,char *fparstr,char *varstr,int ffound,
             int np,int npok,char *pstr,char *ntxt);
void setpvald(double *parmd,char *varstr,int ffound,int np,
              int npok,double *pval,char *punit,char *ntxt);
void setpvali(int *parmi,char *varstr,int ffound,int np,
              int npok,double *pval,char *punit,char *ntxt);
void setpvaldv(double *parmdv,char *varstr,int ffound,int np,
              int npok,double *pval,char *punit,char *ntxt);
void removecr(char c, char *str);
/* solver.c */
void odeint(double *ystart,int nvar,double x1, double x2, 
     double eps,double h1,double hmin,int *nok,
     int *nbad,void (*derivs)(),void (*stiff)());
void stiff(double *y,double *dydx,int n,double *x,double htry,
     double eps,double *yscal,double *hdid,double *hnext,
     void (*derivs)()); 
/* thmfun.c */
void thmfun(double *y,double *yp,int fconvl,double thcl,
            double tsol,double ttol,double *mvl,
            double *mhdl,double *vb,double ga,double ta,
            double gi,double ti);
/* write.c */
void writedat(double *tmv,double **yy,int lt,double *tcb0,
              double **spm,double cacl,double mgcl,double s4cl,
              double *dsv,double *zv,int nz,
              int svrestart,char *fpsvstr);
char *getocnstr(int n);
char *getco2str(int n);
char *getfcstr(int n);
