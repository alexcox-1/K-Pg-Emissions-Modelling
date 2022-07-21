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
#include "defs.h"
#include "utils.h"

#include "common.h"

/*============================================================*/
/*==================== reademiss() ===========================*/
/*============================================================*/
/*
 updates: 
  11/12/20 Add readSemiss()
  04/03/11 deleted femiss() (replaced by finterp())
  03/26/11 new file
	  
 */

void reademiss()
{
 double *tmp1,*tmp2;
 int i=1,rdflag=0;
 char mssg[BUFSIZ];	

 FILE *fpems;

 ltem = 0;
 tmp1 = dvector(1,NEMSMAX);
 tmp2 = dvector(1,NEMSMAX);
	 
 /* open emission file */
 fpems = fopen(ffldstr,"r");
 if(fpems == NULL){ 
    sprintf(mssg,"reademiss(): Can't open emission file '%s'",ffldstr);
	ferrx(mssg);
 } 

 /* read emissions */
 while(rdflag != EOF){
    rdflag = fscanf(fpems,"%le %le",&tmp1[i],&tmp2[i]);
    if(rdflag == 1)
         ferrx("reademiss(): Emission file: #values/line read = 1. need 2.");		 
    if(rdflag  > 2)
         ferrx("reademiss(): Emission file: #values/line read > 2. need 2.");		 
    if(rdflag == 2){
	  ltem = i;
      /*printf("%d %e %e\n",ltem,tmp1[i],tmp2[i]);*/
	  if(ltem > NEMSMAX)
         ferrx("reademiss(): Too many emission values. Increase NEMSMAX?");
      i += 1;
	} else { /* EOF */ 
		 /* printf("%d %d",rdflag,EOF); */
	     break;
	}
 }
	 
 /* close emission file */
 fclose(fpems);

 if(ltem == 0)
    ferrx("reademiss(): No. of emission values = 0?");

 /* allocate time and emission vectors */	 
 tems  = dvector(1,ltem);
 yems  = dvector(1,ltem);
 /* free: see initfree() */

 /* copy numbers read from file to time and emission vectors */	 
 for(i=1;i<=ltem;i++){
	 tems[i] = tmp1[i];
	 yems[i] = tmp2[i];
     /*printf("%e %e\n",tems[i],yems[i]);*/
 }
	
 free_dvector(tmp1,1,NEMSMAX);
 free_dvector(tmp2,1,NEMSMAX);
	
}

void readSemiss()
{
 double *tSmp1,*tSmp2;
 int i=1,rdflag=0;
 char mssg[BUFSIZ];	

 FILE *Spems;

 ltSem = 0;
 tSmp1 = dvector(1,NEMSMAX);
 tSmp2 = dvector(1,NEMSMAX);
	 
 /* open emission file */
 Spems = fopen(sldstr,"r");
 if(Spems == NULL){ 
    sprintf(mssg,"reademiss(): Can't open emission file '%s'",sldstr);
	ferrx(mssg);
 } 

 /* read emissions */
 while(rdflag != EOF){
    rdflag = fscanf(Spems,"%le %le",&tSmp1[i],&tSmp2[i]);
    if(rdflag == 1)
         ferrx("reademiss(): S Emission file: #values/line read = 1. need 2.");		 
    if(rdflag  > 2)
         ferrx("reademiss(): S Emission file: #values/line read > 2. need 2.");		 
    if(rdflag == 2){
	  ltSem = i;
      /*printf("%d %e %e\n",ltem,tmp1[i],tmp2[i]);*/
	  if(ltSem > NEMSMAX)
         ferrx("reademiss(): Too many emission values. Increase NEMSMAX?");
      i += 1;
	} else { /* EOF */ 
		 /* printf("%d %d",rdflag,EOF); */
	     break;
	}
 }
	 
 /* close emission file */
 fclose(Spems);

 if(ltSem == 0)
    ferrx("reademiss(): No. of emission values = 0?");

 /* allocate time and emission vectors */	 
 tSems  = dvector(1,ltSem);
 ySems  = dvector(1,ltSem);
 /* free: see initfree() */

 /* copy numbers read from file to time and emission vectors */	 
 for(i=1;i<=ltSem;i++){
	 tSems[i] = tSmp1[i];
	 ySems[i] = tSmp2[i];
     /*printf("%e %e\n",tems[i],yems[i]);*/
 }
	
 free_dvector(tSmp1,1,NEMSMAX);
 free_dvector(tSmp2,1,NEMSMAX);
	
}

void readExp()
{
 double *tExp1,*tExp2;
 int i=1,rdflag=0;
 char mssg[BUFSIZ];	

 FILE *Exppems;

 ltExp = 0;
 tExp1 = dvector(1,NEMSMAX);
 tExp2 = dvector(1,NEMSMAX);
	 
 /* open emission file */
 Exppems = fopen(expldstr,"r");
 if(Exppems == NULL){ 
    sprintf(mssg,"reademiss(): Can't open emission file '%s'",expldstr);
	ferrx(mssg);
 } 

 /* read emissions */
 while(rdflag != EOF){
    rdflag = fscanf(Exppems,"%le %le",&tExp1[i],&tExp2[i]);
    if(rdflag == 1)
         ferrx("reademiss(): Export Reduction file: #values/line read = 1. need 2.");		 
    if(rdflag  > 2)
         ferrx("reademiss(): Export Reduction file: #values/line read > 2. need 2.");		 
    if(rdflag == 2){
	  ltExp = i;
      /*printf("%d %e %e\n",ltem,tmp1[i],tmp2[i]);*/
	  if(ltExp > NEMSMAX)
         ferrx("reademiss(): Too many emission values. Increase NEMSMAX?");
      i += 1;
	} else { /* EOF */ 
		 /* printf("%d %d",rdflag,EOF); */
	     break;
	}
 }
	 
 /* close emission file */
 fclose(Exppems);

 if(ltExp == 0)
    ferrx("reademiss(): No. of emission values = 0?");

 /* allocate time and emission vectors */	 
 tExp  = dvector(1,ltExp);
 yExp  = dvector(1,ltExp);
 /* free: see initfree() */

 /* copy numbers read from file to time and emission vectors */	 
 for(i=1;i<=ltExp;i++){
	 tExp[i] = tExp1[i];
	 yExp[i] = tExp2[i];
     /*printf("%e %e\n",tems[i],yems[i]);*/
 }
	
 free_dvector(tExp1,1,NEMSMAX);
 free_dvector(tExp2,1,NEMSMAX);
	
}

void readRemin()
{
 double *tRemin1,*tRemin2;
 int i=1,rdflag=0;
 char mssg[BUFSIZ];	

 FILE *Reminpems;

 ltRemin = 0;
 tRemin1 = dvector(1,NEMSMAX);
 tRemin2 = dvector(1,NEMSMAX);
	 
 /* open emission file */
 Reminpems = fopen(reminldstr,"r");
 if(Reminpems == NULL){ 
    sprintf(mssg,"reademiss(): Can't open emission file '%s'",reminldstr);
	ferrx(mssg);
 } 

 /* read emissions */
 while(rdflag != EOF){
    rdflag = fscanf(Reminpems,"%le %le",&tRemin1[i],&tRemin2[i]);
    if(rdflag == 1)
         ferrx("reademiss(): Export Reduction file: #values/line read = 1. need 2.");		 
    if(rdflag  > 2)
         ferrx("reademiss(): Export Reduction file: #values/line read > 2. need 2.");		 
    if(rdflag == 2){
	  ltRemin = i;
      /*printf("%d %e %e\n",ltem,tmp1[i],tmp2[i]);*/
	  if(ltRemin > NEMSMAX)
         ferrx("reademiss(): Too many emission values. Increase NEMSMAX?");
      i += 1;
	} else { /* EOF */ 
		 /* printf("%d %d",rdflag,EOF); */
	     break;
	}
 }
	 
 /* close emission file */
 fclose(Reminpems);

 if(ltRemin == 0)
    ferrx("reademiss(): No. of emission values = 0?");

 /* allocate time and emission vectors */	 
 tRemin  = dvector(1,ltRemin);
 yRemin  = dvector(1,ltRemin);
 /* free: see initfree() */

 /* copy numbers read from file to time and emission vectors */	 
 for(i=1;i<=ltRemin;i++){
	 tRemin[i] = tRemin1[i];
	 yRemin[i] = tRemin2[i];
     /*printf("%e %e\n",tems[i],yems[i]);*/
 }
	
 free_dvector(tRemin1,1,NEMSMAX);
 free_dvector(tRemin2,1,NEMSMAX);
	
}
/*============================================================*/
/*==================== reademiss() END =======================*/
/*============================================================*/
