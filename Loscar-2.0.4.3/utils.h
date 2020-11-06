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

/*============================================================*/
/*==================== utils.h ===============================*/
/*============================================================*/
/*
 updates: 

  04/21/11 deleted #ifndef _NR_UTILS_H_
  04/03/11 deleted femiss()
  03/02/11 deleted #if defined ANSI 
  02/26/11 deleted #define NRALL block
  02/22/11 new file

  Error/Warning handlers:
	  
  NR:  nrerror fatal (exit)
  REZ: ferrx   fatal (exit)
  REZ: fwarn   issue warning 
  REZ: ferrwrt write data & exit, see loscarderivs.c
	  
 */

void nrerror(char error_text[]);
int      *ivector(long nl, long nh);
double   *dvector(long nl, long nh);
double  **dmatrix(long nrl, long nrh, long ncl, long nch);
int     **imatrix(long nrl, long nrh, long ncl, long nch);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, 
                   long ndh);

void free_ivector(int *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m,long nrl,long nrh,long ncl,long nch);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);


/*===================== REZ additions ================*/

/* error output */
void ferrx(char error_text[]);

/* warning output */
void fwarn(char warn_text[]);

/*** vsclr: vector X scalar  , returns vector ******/
void vsclr(double *w,double s,double *u,int n);

/*** vvdot: vector product   , returns scalar ******/
void vvdot(double *s,double *u,double *v,int n);

/*** vsumlh: sum vector elements, returns scalar **/
double vsumlh(double *u,int n,int nl,int nh);

/*** msumlh: sum matrix column elements, returns scalar **/
double msumlh(double **m,int krmax,int k,int n,int nl,int nh);

/*** vvsum: vector sum, returns vector *************/
void vvsum(double *w,double *u,double *v,int n);

/*** vvsub: vector subtraction, returns vector *****/
void vvsub(double *w,double *u,double *v,int n);

/*** vvelm: v-elmnt X v-elmnt, returns vector ******/
void vvelm(double *w,double *u,double *v, int n);

/*** vdelm: v-elmnt / v-elmnt, returns vector ******/
void vdelm(double *w,double *u,double *v, int n);

/*** fsetiv: init int vector elements , returns vector **/
void fsetiv(int *iv,int k,int n);

/*** fsetv: init double vector elements , returns vector **/
void fsetv(double *u,double x,int n);

/*** fvec2arr: store 2 vectors in matrix, returns matrix ******/
void fvec2arr(double **m,double *x,double *y,int n);

/*** dfind: finds vector entries <,<=,==,>=,> z ****/
/*   CAUTION!!! This function can produce unexpected 
     results, see utils.c                          */
void dfind(int *indxv,double *u,int n,int opr, double z, 
           double mz);

/*** fminv: find vector min, returns int ***/
int fminv(double *u,int n,int ordflag);

/*** finterp: linear interpol., returns double ***/
double finterp(double x,double *xv,double *yv,int n,int errflag);
