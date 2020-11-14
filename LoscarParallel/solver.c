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
/*==================== odeint() ==============================*/
/*============================================================*/
/*
 updates: 

  11/30/11 if: hdid == h  changed to if: (fabs(hdid - h))  <= DBL_EPSILON
	       same results (either h=hdid or not). still safer with eps.
	       if: *x == xsav changed to if: (fabs(*x - xsav)) <= DBL_EPSILON
  10/30/11 removed macro SOLVER -> stiff
  04/08/11 moved MAXSTP to defs.h
  03/02/11 #define MAXSTP 100000 (10000 before)
  02/26/11 initialized xsav=0.0
  02/14/11 new file
 

 Runge-Kutta driver. calls derivs and SOLVER.
 Runge-Kutta driver with adaptive stepsize control. Integrate starting 
 values ystart[1..nvar] from x1 to x2 with accuracy eps, storing 
 intermediate results in global variables. h1 should be set as a 
 guessed first stepsize, hmin as the minimum allowed stepsize (can be 
 zero). On output nok and nbad are the number of good and bad (but retried 
 and fixed) steps taken, and ystart is replaced by values at the end of 
 the integration interval. derivs is the user-supplied routine for 
 calculating the right-hand side derivative, while SOLVER is the name 
 of the stepper routine to be used.
 */	 

#define TINY 1.0e-30

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

extern int kmax,kount;
extern double *tmv,**yy,dxsav;

/* User storage for intermediate results. Preset kmax and dxsav 
  in the calling program. If kmax > 0, results are stored at approximate 
  intervals dxsav in the arrays tmv[1..kount], yy[1..nvar] [1..kount], 
  where kount is output by odeint. Defining declarations for these 
  variables, with memory allocations tmv[1..kmax] and yy[1..nvar][1..kmax] 
  for the arrays, should be in the calling program. 
*/

void odeint(double *ystart,int nvar,double x1, double x2, 
     double eps,double h1,double hmin,int *nok,
     int *nbad,void (*derivs)(),void (*stiff)())
{
	/* int k; */
    double *dfdx,**dfdy;
		
	int nstp,i;
	double xsav=0.0,x,hnext,hdid,h;
	double *yscal,*y,*dydx;
    double xcount=0.0; /* REZ */
	
    dfdx=dvector(1,nvar);
    dfdy=dmatrix(1,nvar,1,nvar);	
	
	yscal=dvector(1,nvar);
	y=dvector(1,nvar);
	dydx=dvector(1,nvar);
	x=x1;
	h=SIGN(h1,x2-x1);
	*nok = (*nbad) = kount = 0;
	for (i=1;i<=nvar;i++) y[i]=ystart[i];
	if (kmax > 0) xsav=x-dxsav*2.0;
	/* REZ*/
    printf("\n@ Progress (model time): \n\n");
	for (nstp=1;nstp<=MAXSTP;nstp++) {
		/* REZ */
		if((x-x1) >= xcount*(x2-x1)){
		     printf("%.2e ",x);
 			 fflush(stdout);
			 xcount += 0.05;
		}
		(*derivs)(x,y,dydx);
		for (i=1;i<=nvar;i++)
			/* REZ error scaling */
			/* yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY; */
            /* max(C,|y|) set C !!! */			
			yscal[i]=DMAX(CSCAL,fabs(y[i]));
		if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
			tmv[++kount]=x;
			for (i=1;i<=nvar;i++) yy[i][kount]=y[i];
			xsav=x;
		}
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
		/* SOLVER */
		(*stiff)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
		if ((fabs(hdid - h)) <= DBL_EPSILON){ /* REZ */
			++(*nok); 
			/* printf("%d %.40e ",*nok,DBL_EPSILON); */ 
		} else ++(*nbad);
		if ((x-x2)*(x2-x1) >= 0.0) {
			for (i=1;i<=nvar;i++) ystart[i]=y[i];
			if (kmax) {
				tmv[++kount]=x;
				for (i=1;i<=nvar;i++) yy[i][kount]=y[i];
			}
			free_dvector(dydx,1,nvar);
			free_dvector(y,1,nvar);
			free_dvector(yscal,1,nvar);
			return;
		}
		if (fabs(hnext) <= hmin) ferrwrt("Step size too small in odeint");
		h=hnext;
		/* JACOBN */
		jacobn(x,y,dfdx,dfdy,nvar);
/*		for(i=1;i<=nvar;i++)
   		   printf("%7.4f dfdx \n",dfdx[i]);
		for(i=1;i<=nvar;i++){
		   for(k=1;k<=nvar;k++)
		      printf("%7.4f ",dfdy[i][k]);
		   printf("\n");
		} 	*/
	}
	ferrwrt("Too many steps in routine odeint");

	/* JACOBN */
	free_dmatrix(dfdy,1,nvar,1,nvar);
    free_dvector(dfdx,1,nvar);
}

#undef TINY
/*============================================================*/
/*==================== odeint() END ==========================*/
/*============================================================*/



/*============================================================*/
/*==================== stiff() ===============================*/
/*============================================================*/
/* Fourth-order Rosenbrock step.

 void stiff(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs) 

 Fourth-order Rosenbrock step for integrating stiff o.d.e.'s, with 
 monitoring of local truncation error to adjust stepsize. Input are 
 the dependent variable vector y[1..n] and its derivative dydx[1..n] 
 at the starting value of the independent variable x. Also input 
 are the stepsize to be attempted htry, the required accuracy eps, 
 and the vector yscal[1..n] against which the error is scaled. On 
 output, y and x are replaced by their new values, hdid is the stepsize
 that was actually accomplished, and hnext is the estimated next 
 stepsize. derivs is a usersupplied routine that computes the 
 derivatives of the right-hand side with respect to x, while 
 jacobn (a fixed name) is a user-supplied routine that computes the 
 Jacobi matrix of derivatives of the right-hand side with respect 
 to the components of y.
*/	 

#define NRANSI 
#define SAFETY 0.9 
#define GROW 1.5 
#define PGROW -0.25 
#define SHRNK 0.5 
#define PSHRNK (-1.0/3.0) 
#define ERRCON 0.1296 
#define MAXTRY 40 
#define GAM (1.0/2.0) 
#define A21 2.0 
#define A31 (48.0/25.0) 
#define A32 (6.0/25.0) 
#define C21 -8.0 
#define C31 (372.0/25.0) 
#define C32 (12.0/5.0) 
#define C41 (-112.0/125.0) 
#define C42 (-54.0/125.0) 
#define C43 (-2.0/5.0) 
#define B1 (19.0/9.0) 
#define B2 (1.0/2.0) 
#define B3 (25.0/108.0) 
#define B4 (125.0/108.0) 
#define E1 (17.0/54.0) 
#define E2 (7.0/36.0) 
#define E3 0.0 
#define E4 (125.0/108.0) 
#define C1X (1.0/2.0) 
#define C2X (-3.0/2.0) 
#define C3X (121.0/50.0) 
#define C4X (29.0/250.0) 
#define A2X 1.0 
#define A3X (3.0/5.0) 
 
void stiff(double *y,double *dydx,int n,double *x,double htry,
     double eps,double *yscal,double *hdid,double *hnext,
     void (*derivs)())
{ 

	int i,j,jtry,*indx; 
	double d,errmax,h,xsav,**a,*dfdx,**dfdy,*dysav,*err; 
	double *g1,*g2,*g3,*g4,*ysav; 
 
	indx=ivector(1,n); 
	a=dmatrix(1,n,1,n); 
	dfdx=dvector(1,n); 
	dfdy=dmatrix(1,n,1,n); 
	dysav=dvector(1,n); 
	err=dvector(1,n); 
	g1=dvector(1,n); 
	g2=dvector(1,n); 
	g3=dvector(1,n); 
	g4=dvector(1,n); 
	ysav=dvector(1,n); 
	xsav=(*x); 
	for (i=1;i<=n;i++) { 
		ysav[i]=y[i]; 
		dysav[i]=dydx[i]; 
	} 
	jacobn(xsav,ysav,dfdx,dfdy,n); 
	h=htry; 
	for (jtry=1;jtry<=MAXTRY;jtry++) { 
		for (i=1;i<=n;i++) { 
			for (j=1;j<=n;j++) a[i][j] = -dfdy[i][j]; 
			a[i][i] += 1.0/(GAM*h); 
		} 
		ludcmp(a,n,indx,&d); 
		for (i=1;i<=n;i++) 
			g1[i]=dysav[i]+h*C1X*dfdx[i]; 
		lubksb(a,n,indx,g1); 
		for (i=1;i<=n;i++) 
			y[i]=ysav[i]+A21*g1[i]; 
		*x=xsav+A2X*h; 
		(*derivs)(*x,y,dydx); 
		for (i=1;i<=n;i++) 
			g2[i]=dydx[i]+h*C2X*dfdx[i]+C21*g1[i]/h; 
		lubksb(a,n,indx,g2); 
		for (i=1;i<=n;i++) 
			y[i]=ysav[i]+A31*g1[i]+A32*g2[i]; 
		*x=xsav+A3X*h; 
		(*derivs)(*x,y,dydx); 
		for (i=1;i<=n;i++) 
			g3[i]=dydx[i]+h*C3X*dfdx[i]+(C31*g1[i]+C32*g2[i])/h; 
		lubksb(a,n,indx,g3); 
		for (i=1;i<=n;i++) 
			g4[i]=dydx[i]+h*C4X*dfdx[i]+(C41*g1[i]+C42*g2[i]+C43*g3[i])/h; 
		lubksb(a,n,indx,g4); 
		for (i=1;i<=n;i++) { 
			y[i]=ysav[i]+B1*g1[i]+B2*g2[i]+B3*g3[i]+B4*g4[i]; 
			err[i]=E1*g1[i]+E2*g2[i]+E3*g3[i]+E4*g4[i]; 
		} 
		*x=xsav+h;
		if ((fabs(*x - xsav)) <= DBL_EPSILON) /* REZ */
			ferrwrt("stepsize not significant in stiff"); 
		errmax=0.0; 
		for (i=1;i<=n;i++) errmax=DMAX(errmax,fabs(err[i]/yscal[i])); 
		errmax /= eps; 
		if (errmax <= 1.0) { 
			*hdid=h; 
			*hnext=(errmax > ERRCON ? SAFETY*h*pow(errmax,PGROW) : GROW*h); 
			free_dvector(ysav,1,n); 
			free_dvector(g4,1,n); 
			free_dvector(g3,1,n); 
			free_dvector(g2,1,n); 
			free_dvector(g1,1,n); 
			free_dvector(err,1,n); 
			free_dvector(dysav,1,n); 
			free_dmatrix(dfdy,1,n,1,n); 
			free_dvector(dfdx,1,n); 
			free_dmatrix(a,1,n,1,n); 
			free_ivector(indx,1,n); 
			return; 
		} else { 
			*hnext=SAFETY*h*pow(errmax,PSHRNK); 
			h=(h >= 0.0 ? DMAX(*hnext,SHRNK*h) : DMIN(*hnext,SHRNK*h)); 
		} 
	} 
	ferrwrt("exceeded MAXTRY in stiff"); 
} 
#undef SAFETY 
#undef GROW 
#undef PGROW 
#undef SHRNK 
#undef PSHRNK 
#undef ERRCON 
#undef MAXTRY 
#undef GAM 
#undef A21 
#undef A31 
#undef A32 
#undef C21 
#undef C31 
#undef C32 
#undef C41 
#undef C42 
#undef C43 
#undef B1 
#undef B2 
#undef B3 
#undef B4 
#undef E1 
#undef E2 
#undef E3 
#undef E4 
#undef C1X 
#undef C2X 
#undef C3X 
#undef C4X 
#undef A2X 
#undef A3X 
#undef NRANSI 
/*============================================================*/
/*==================== stiff() END ===========================*/
/*============================================================*/

