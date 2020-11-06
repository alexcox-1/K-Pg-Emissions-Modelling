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
#include <float.h>
#include "defs.h"
#include "utils.h"

/*============================================================*/
/*==================== ludcmp() ==============================*/
/*============================================================*/
/*
 updates:

 12/01/11 changed if: a[j][j] == 0.0 to if: fabs(a[j][j]) < TINY
          changed if: big == 0.0 to if: fabs(big) <= DBL_EPSILON
 04/21/11 deleted #define NRANSI
 02/26/11 initialized imax=0
 02/14/11 new file
 
 */

#define TINY (1.0e-20) /* NR: semicolon ?*/ 
 
void ludcmp(double **a,int n,int *indx,double *d)
{ 
	int i,imax=0,j,k; 
	double big,dum,sum,temp; 
	double *vv; 
 
	vv=dvector(1,n); 
	*d=1.0; 
	for (i=1;i<=n;i++) { 
		big=0.0; 
		for (j=1;j<=n;j++) 
			if ((temp=fabs(a[i][j])) > big) big=temp; 
		if (fabs(big) <= DBL_EPSILON) 
			ferrwrt("Singular matrix in routine ludcmp"); 
		vv[i]=1.0/big; 
	} 
	for (j=1;j<=n;j++) { 
		for (i=1;i<j;i++) { 
			sum=a[i][j]; 
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j]; 
			a[i][j]=sum; 
		} 
		big=0.0; 
		for (i=j;i<=n;i++) { 
			sum=a[i][j]; 
			for (k=1;k<j;k++) 
				sum -= a[i][k]*a[k][j]; 
			a[i][j]=sum; 
			if ( (dum=vv[i]*fabs(sum)) >= big) { 
				big=dum; 
				imax=i; 
			} 
		} 
		if (j != imax) { 
			for (k=1;k<=n;k++) { 
				dum=a[imax][k]; 
				a[imax][k]=a[j][k]; 
				a[j][k]=dum; 
			} 
			*d = -(*d); 
			vv[imax]=vv[j]; 
		} 
		indx[j]=imax; 
		if ((fabs(a[j][j])) < TINY) a[j][j]=TINY; /* REZ */
		if (j != n) { 
			dum=1.0/(a[j][j]); 
			for (i=j+1;i<=n;i++) a[i][j] *= dum; 
		} 
	} 
	free_dvector(vv,1,n); 
} 
#undef TINY 

/*============================================================*/
/*==================== ludcmp() END ==========================*/
/*============================================================*/


/*============================================================*/
/*==================== lubksb() ==============================*/
/*============================================================*/
void lubksb(double **a,int n,int *indx,double *b)
{ 
	int i,ii=0,ip,j; 
	double sum; 
 
	for (i=1;i<=n;i++) { 
		ip=indx[i]; 
		sum=b[ip]; 
		b[ip]=b[i]; 
		if (ii) 
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j]; 
		else if (sum) ii=i; 
		b[i]=sum; 
	} 
	for (i=n;i>=1;i--) { 
		sum=b[i]; 
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j]; 
		b[i]=sum/a[i][i]; 
	} 
} 
/*============================================================*/
/*==================== lubksb() END ==========================*/
/*============================================================*/

