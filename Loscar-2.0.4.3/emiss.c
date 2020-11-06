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
/*============================================================*/
/*==================== reademiss() END =======================*/
/*============================================================*/
