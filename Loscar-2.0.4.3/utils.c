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
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "utils.h"

/*============================================================*/
/*==================== utils.c ===============================*/
/*============================================================*/
/* based on nrutils.c. additions at the end                   */
/*
 updates: 

  12/03/11 finterp(): changed xv[] == x to (fabs(xv[]-x)) <= DBEPS
           fminv(): u[i] < umin to (double)(u[i]) < (double)(umin)
  12/01/11 revised double comparison	
           changed '==' etc. using DBL_EPSILON
  09/15/11 added fflush(stdout), fflush(stderr).
  04/21/11 added search order option in fminv()
  04/03/11 added fminv() and finterp()
  03/09/11 added dfind() 
  03/02/11 deleted #if defined ANSI 
  02/22/11 new file

  Error/Warning handlers:
	  
  NR:  nrerror fatal (exit)
  REZ: ferrx   fatal (exit)
  REZ: fwarn   issue warning 
  REZ: ferrwrt write data & exit, see loscarderivs.c
	  
 */

#define NR_END 1
#define FREE_ARG char*


void nrerror(char error_text[])
/* NR standard error handler */
{
    fflush(stdout);
    fflush(stderr);	
	fprintf(stderr,"NR run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;
	
    v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}


double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}


double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, 
                   long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;

	/* allocate pointers to pointers to rows */
	t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}



void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
	if(nh < nl) nrerror("free_ivector(): nh < nl");
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
	if(nh < nl) nrerror("free_dvector(): nh < nl");
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
	if(nrh < nrl) nrerror("free_dmatrix(): nrh < nrl");
	if(nch < ncl) nrerror("free_dmatrix(): nch < ncl");	
}

void free_imatrix(int **m,long nrl,long nrh,long ncl,long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
	if(nrh < nrl) nrerror("free_imatrix(): nrh < nrl");
	if(nch < ncl) nrerror("free_imatrix(): nch < ncl");	
}

void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a double f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
	if(nrh < nrl) nrerror("free_d3tensor(): nrh < nrl");
	if(nch < ncl) nrerror("free_d3tensor(): nch < ncl");	
	if(ndh < ndl) nrerror("free_d3tensor(): ndh < ndl");	
}

/*===================== REZ additions ================*/

 /* error output */
 void ferrx(char error_text[])
 {
    fflush(stdout);
    fflush(stderr);	 
	fprintf(stderr,"\n");
	fprintf(stderr,"@ ======================== Error.\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"@ ====================== Exiting.\n");
	fprintf(stderr,"  Morituri te salutant!\n");
	exit(1);
 }

 /* warning */
 void fwarn(char warn_text[])
 {
    fflush(stdout);
    fflush(stderr);
	fprintf(stderr,"\n");
	fprintf(stderr,"@ ======================== Warning.\n");
	fprintf(stderr,"%s\n",warn_text);
	fprintf(stderr,"@ ===================== Continuing.\n");
    fprintf(stderr,"  Cogito ergo sum.\n");	
    fflush(stderr);
 }

 /****** vsclr: vector X scalar, returns vector ******/
 void vsclr(double *w,double s,double *u,int n)
 {
  int i;
 
  for(i=1;i<=n;i++)
    w[i] = s*u[i];
 }

 /****** vvdot: vector product, returns scalar *******/
 void vvdot(double *s,double *u,double *v,int n)
 {
  int i;

  *s = 0.;
  for(i=1;i<=n;i++)
    *s += u[i]*v[i];
 }

 /*** vsumlh: sum vector elements, returns scalar **/
 double vsumlh(double *u,int n,int nl,int nh)
 {
  int i;
  double s;

  if(nh > n)
     ferrx("vsumlh(): nh exceeds vector length n.");
	 
  s = 0.;
  for(i=nl;i<=nh;i++)
    s += u[i];

  return(s);	 
 }

 /*** msumlh: sum matrix column elements, returns scalar **/
 double msumlh(double **m,int krmax,int k,int n,int nl,int nh)
 {
 /* k = row, n: column */	 
  int i;
  double s;

  if(k > krmax)
     ferrx("msumlh(): k exceeds matrix row length krmax.");
  if(nh > n)
     ferrx("msumlh(): nh exceeds matrix column length n.");
	 
  s = 0.;
  for(i=nl;i<=nh;i++)
    s += m[k][i];

  return(s);	 	 
 }


 /****** vvsum: vector sum, returns vector *******/
 void vvsum(double *w,double *u,double *v,int n)
 {
  int i;

  for(i=1;i<=n;i++)
    w[i] = u[i] + v[i];
 }

 /****** vvsub: vector subtraction, returns vector *******/
 void vvsub(double *w,double *u,double *v,int n)
 {
  int i;

  for(i=1;i<=n;i++)
    w[i] = u[i] - v[i];
 }

 /****** vvelm: v-elmnt X v-elmnt, returns vector *******/
 void vvelm(double *w,double *u,double *v, int n)
 {
  int i;

  for(i=1;i<=n;i++)
    w[i] = u[i]*v[i];
 }

 /****** vdelm: v-elmnt / v-elmnt, returns vector *******/
 void vdelm(double *w,double *u,double *v, int n)
 {
  int i;

  for(i=1;i<=n;i++)
    w[i] = u[i]/v[i];
 }

 /*** fsetiv: init int vector elements , returns vector ******/
 void fsetiv(int *iv,int k,int n)
 {
  int i;

  for(i=1;i<=n;i++)
    iv[i] = k;
 }

 /*** fsetv: init double vector elements , returns vector ******/
 void fsetv(double *u,double x,int n)
 {
  int i;

  for(i=1;i<=n;i++)
    u[i] = x;
 }

 /*** fvec2arr: store 2 vectors in matrix, returns matrix ******/
 void fvec2arr(double **m,double *x,double *y,int n)
 {
  int i;

  for(i=1;i<=n;i++){
       m[i][1] = x[i];
       m[i][2] = y[i];
  }		
 } 

#define DBEPS DBL_EPSILON
 /*** dfind: finds vector entries <,<=,==,>=,> z ****/
 /* mz = z's order of magnitude must be supplied by user.
 
 CAUTION!!! This function can produce unexpected results
 because doubles are compared using a simple method. Make
 sure you understand the pitfalls before using. See

 http://gcc.gnu.org/bugzilla/show_bug.cgi?id=323
 
 Bruce Dawson: Comparing floating point numbers
 http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm

 David Goldberg: What Every Computer Scientist Should Know About 
 Floating-Point Arithmetic, ACM, 1991.
 http://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html

 */
 void dfind(int *indxv,double *u,int n,int opr,double z, double mz)
 {
  int i;

  if( (double)(mz) <= (double)(0.0) )
     ferrx("dfind(): set mz > 0.0");
		 
  switch(opr){
     /******** <  *******/
	 case 1:
      for(i=1;i<=n;i++){
          if( ((u[i]-z)/mz) < -DBEPS ) 
			            indxv[i] = i;
          else          indxv[i] = 0;
	  }
      break;		    
     /******** <= *******/
	 case 2:
      for(i=1;i<=n;i++){ 
          if( ((u[i]-z)/mz) < -DBEPS || (fabs((u[i]-z)/mz)) <= DBEPS )
			            indxv[i] = i;
          else          indxv[i] = 0;
	  }
      break;		    
     /******** == *******/
	 case 3:
      for(i=1;i<=n;i++){
          if( (fabs((u[i]-z)/mz)) <= DBEPS )
			            indxv[i] = i;
          else          indxv[i] = 0;
	  }
      break;		    
     /******** >= *******/
	 case 4:
      for(i=1;i<=n;i++){
          if( ((u[i]-z)/mz) > DBEPS || (fabs((u[i]-z)/mz)) <= DBEPS )
                        indxv[i] = i;
          else          indxv[i] = 0;
	  }
      break;		    
     /******** >  *******/
	 case 5:
      for(i=1;i<=n;i++){
          if( ((u[i]-z)/mz) > DBEPS ) 
			            indxv[i] = i;
          else          indxv[i] = 0;
	  }
      break;		    
     /***************/
	 default:	
	  printf("\ndfind(): Operator not defined\n");
      for(i=1;i<=n;i++) indxv[i] = 0;
	  break;		  
  }
 }

 /*** fminv: find vector min, returns int ***/
 int fminv(double *u,int n,int ordflag)
 {
  int i,imin;
  double umin;

  /* search order 1,..,n */	 
  if(ordflag == 1){
    imin =   1;
    umin = u[1];

    for(i=2;i<=n;i++){
        if((double)(u[i]) < (double)(umin)){ 
           imin =   i;
           umin = u[i];
	    }
    }
  }
	  
  /* search order n,..,1 */	 
  if(ordflag == 2){
    imin =   n;
    umin = u[n];

    for(i=n-1;i>=1;i--){
        if((double)(u[i]) < (double)(umin)){ 
		   imin =   i;
		   umin = u[i];
	    }
    }
  }
	 
  return(imin);

 }

 /*** finterp: linear interpol., returns double ***/
 double finterp(double x,double *xv,double *yv,int n,int errflag)
 {
  int i,k=0;
  double dx,b,y;	

  if( (fabs(x-xv[1])) <= DBEPS )
       return(yv[1]);	
  if( (fabs(x-xv[n])) <= DBEPS )
       return(yv[n]);
	 
  if(x > xv[1] && x < xv[n]){
    for(i=1;i<=n;i++){
        if( (fabs(xv[i]-x)) <= DBEPS )
			return(yv[i]); 
        if(xv[i] >  x){
			k = i-1;
            break;
		}
	}
	dx = x - xv[k];
	b  = (yv[k+1]-yv[k])/(xv[k+1]-xv[k]);
    y  = yv[k] + b*dx;
    return(y);	  
  }

  if(x < xv[1] || x > xv[n]){
   if(errflag == 0)
     return(0.0);
   if(errflag == 1)
     ferrx("finterp(): Input x is out of bounds");
   if(errflag == 2){
	   if(x < xv[1])
		   return(yv[1]);
	   if(x > xv[n])
		   return(yv[n]);
   }
  }
   return(0.0);

 }
#undef DBEPS

