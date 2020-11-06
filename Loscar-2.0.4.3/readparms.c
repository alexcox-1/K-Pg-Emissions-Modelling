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
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "utils.h"

#include "common.h"

/*============================================================*/
/*==================== readparms() ===========================*/
/*============================================================*/
/* read parameters from input control file                    */
/*

   updates:

   09/03/18 new K* corrections (ZeebeTyrrell18)
   10/31/11 added sclim
   10/28/11 added d13C carbon input, tcin0, tcspan
   10/23/11 added carbon input
   10/21/11 including Tethys. read thc0, set thc = thc0
   04/16/11 added setpvald(), setpvali(), setpvaldv()
   04/15/11 temperature: added tcb as tracer (rm tcv). start priorities:
            1. control file 2. restart file 3. intern default
            NOTE: restart is read after cntrl => remember tcb0c
   04/11/11 added setpstr() etc. cleaned up.
   03/03/11 new file

   parameters added:

   04/21/11 FBIOL/CBIOH: Low/High-Lat biopump
   04/19/11 FDICI/FALKI/FPO4I: change initial dic, alk, po4 (fdapi)
   04/16/11 TSNS: OCN temperature change (CO2-sensitivity)


 */

#define PVALMAX 50

void readparms(char *fparmstr)
{
 double *pval;
 int n,np=0,ffound=0,npok=0;
 char pstr[BUFSIZ],*varstr,*punit,*ntxt;
 char mssg[BUFSIZ];		
 FILE *fparm;

 /* vector returning parameter values */	
 pval = dvector(1,PVALMAX); 
	
 /* open parameter file */
 fparm = fopen(fparmstr,"r");
 if(fparm == NULL){ 
    sprintf(mssg,"readparms(): Can't open control file '%s'",fparmstr);
	ferrx(mssg);
 } 
 else{
    printf("\n@ Setting up run. Reading parameter file: '%s'\n\n",fparmstr);
 }

 /*========================================================*/	
 /*==== get parms from control file and set parms/vars. ===*/
 /*==== remember to set arguments 1 (and 2) of setp()   ===*/
 /*==== when adding new control parameters !!!          ===*/
 /*                                                        */
 /* ntxt: display text if not found                        */
	
 /*====================== restart file ====================*/	
 varstr = "RESTART";    /* name of parameter     */
 npok   = 1;	        /* # par values expected */
 ntxt   = "none. Using default start values"; 
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpstr: set arguments 1 and 2 */	
 setpstr(&ldrestart,fpldstr,varstr,ffound,np,npok,pstr,ntxt);

	
 /*================== save restart file ===================*/	
 varstr = "SVSTART";	/* name of parameter     */
 npok   = 1;	        /* # par values expected */
 ntxt   = "none. Not saving restart values";  
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpstr: set arguments 1 and 2 */	
 setpstr(&svrestart,fpsvstr,varstr,ffound,np,npok,pstr,ntxt);

	
 /*=========== load fossil fuel emission file =============*/	
 varstr = "EMSFILE";	/* name of parameter     */
 npok   = 1;	        /* # par values expected */
 ntxt   = "none.";
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpstr: set arguments 1 and 2 */	
 setpstr(&ffflag,ffldstr,varstr,ffound,np,npok,pstr,ntxt);

	
 /*===================== time start =======================*/	
 varstr = "TSTART";	    /* name of parameter     */
 npok   = 1;	        /* # par values expected */
 punit  = "(y)";        /* parameter unit        */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvald: set argument 1 */	
 setpvald(&t0,varstr,ffound,np,npok,pval,punit,ntxt);	

	
 /*====================== time final ======================*/
 varstr = "TFINAL";	    /* name of parameter     */
 npok   = 1;	        /* # par values expected */
 punit  = "(y)";        /* parameter unit        */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvald: set argument 1 */	
 setpvald(&tfinal,varstr,ffound,np,npok,pval,punit,ntxt);	

	
 /*================= carbon input ==========================*/	
 varstr = "CINP";       /* name of parameter     */
 npok   = 1;	        /* # par values expected */
 punit  = "  (Pg C)      ";/* parameter unit     */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvald: set argument 1 */	
 setpvald(&cinp,varstr,ffound,np,npok,pval,punit,ntxt);
 if(ffound == 1)
    cinpflag = 1;	

	
 /*================= d13C carbon input =====================*/	
 varstr = "D13CIN";     /* name of parameter     */
 npok   = 1;	        /* # par values expected */
 punit  = "(per mil)  "; /* parameter unit     */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvald: set argument 1 */	
 setpvald(&dccinp,varstr,ffound,np,npok,pval,punit,ntxt);
 if(ffound == 1)
    rccinp = (dccinp/1.e3+1.)*RST;

 /*================= time start C input ====================*/	
 varstr = "TCIN0";      /* name of parameter     */
 npok   = 1;	        /* # par values expected */
 punit  = " (y)";       /* parameter unit     */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvald: set argument 1 */	
 setpvald(&tcin0,varstr,ffound,np,npok,pval,punit,ntxt);


 /*================= time span C input ====================*/	
 varstr = "TCSPAN";     /* name of parameter     */
 npok   = 1;	        /* # par values expected */
 punit  = "(y)";        /* parameter unit     */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvald: set argument 1 */	
 setpvald(&tcspan,varstr,ffound,np,npok,pval,punit,ntxt);


 /*================= epslvr, solver accuracy ==============*/
 varstr = "EPSLV";	    /* name of parameter     */
 npok   = 1;	        /* # par values expected */
 punit  = " (-)";       /* parameter unit        */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvald: set argument 1 */	
 setpvald(&epslvr,varstr,ffound,np,npok,pval,punit,ntxt);	


 /*======== fcsml, numerics: linear f_calcite drop ========*/
 /*======== during dissolution if fc < fcsml.      ========*/
 varstr = "FCSML";	    /* name of parameter     */
 npok   = 1;	        /* # par values expected */
 punit  = " (-)";       /* parameter unit        */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvald: set argument 1 */	
 setpvald(&fcsml,varstr,ffound,np,npok,pval,punit,ntxt);	
	

 /*======== kmax, desired # output values during run ======*/
 varstr = "KMAX";	    /* name of parameter     */
 npok   = 1;	        /* # par values expected */
 punit  = "  (-)";      /* parameter unit        */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvali: set argument 1 */	
 setpvali(&kmax,varstr,ffound,np,npok,pval,punit,ntxt);	


 /*=========== temperature sensitivity to pCO2 ============*/
 varstr = "TSNS";	    /* name of parameter     */
 npok   = 1;	        /* # par values expected */
 punit  = "  (-)";      /* parameter unit        */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvali: set argument 1 */	
 setpvali(&tsnsflag,varstr,ffound,np,npok,pval,punit,ntxt);	

	
 /*=========== value T-sensitivity ========================*/
 varstr = "SCLIM";      /* name of parameter     */
 npok   = 1;            /* # par values expected */
 punit  = " (C)";  /* parameter unit     */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvaldv: set argument 1 */	
 setpvald(&sclim,varstr,ffound,np,npok,pval,punit,ntxt);	


 /*================== temperature boxes ===================*/
 varstr = "TEMP";       /* name of parameter     */
 npok   = NB;           /* # par values expected */
 punit  = "  (C)";      /* parameter unit        */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /* set default tcb0c = tcb0. if not found tcb0c = empty! */
 for(n=1;n<=NB;n++)
         tcb0c[n] = tcb0[n];
 /*setpvaldv: set argument 1 */	
 setpvaldv(tcb0c,varstr,ffound,np,npok,pval,punit,ntxt);	
		

 /*==================== salinity boxes ====================*/
 varstr = "SAL";        /* name of parameter     */
 npok   = NB;           /* # par values expected */
 punit  = "   (-)";     /* parameter unit        */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvaldv: set argument 1 */	
 setpvaldv(salv,varstr,ffound,np,npok,pval,punit,ntxt);	


 /*==================== pCO2 steady-state =================*/
 varstr = "PCO2SI";     /* name of parameter     */
 npok   = 1;            /* # par values expected */
 punit  = "(ppmv)      ";  /* parameter unit     */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvaldv: set argument 1 */	
 setpvald(&pcsi,varstr,ffound,np,npok,pval,punit,ntxt);	
	
	
/*=================== factor initial DIC ==================*/	
 varstr = "FDICI";      /* name of parameter     */
 npok   = 1;	        /* # par values expected */
 punit  = " (%)";       /* parameter unit        */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvald: set argument 1 */	
 setpvald(&fdapi[1],varstr,ffound,np,npok,pval,punit,ntxt);
	

/*=================== factor initial ALK ==================*/	
 varstr = "FALKI";      /* name of parameter     */
 npok   = 1;	        /* # par values expected */
 punit  = " (%)";       /* parameter unit        */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvald: set argument 1 */	
 setpvald(&fdapi[2],varstr,ffound,np,npok,pval,punit,ntxt);

	
/*=================== factor initial PO4 ==================*/	
 varstr = "FPO4I";      /* name of parameter     */
 npok   = 1;	        /* # par values expected */
 punit  = " (%)";       /* parameter unit        */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvald: set argument 1 */	
 setpvald(&fdapi[3],varstr,ffound,np,npok,pval,punit,ntxt);
 
	
 /*=================== conveyor transport =================*/	
 varstr = "THC";        /* name of parameter     */
 npok   = 1;	        /* # par values expected */
 punit  = "  (Sv)";     /* parameter unit        */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvald: set argument 1 */	
 setpvald(&thc0,varstr,ffound,np,npok,pval,punit,ntxt);
 if(ffound == 1){
    thc0 *= (1.e6*YTOSEC);
    thc   = thc0;
 }


 /*=============== LL Biopump utilization =================*/
 varstr = "FBIOL";      /* name of parameter     */
 npok   = 1;	        /* # par values expected */
 punit  = " (-)";       /* parameter unit        */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvald: set argument 1 */	
 setpvald(&fepl,varstr,ffound,np,npok,pval,punit,ntxt);	


 /*================= HL Biopump C-Export ==================*/
 varstr = "CBIOH";      /* name of parameter     */
 npok   = 1;	        /* # par values expected */
 punit  = " (mol C/m2/y)";/* parameter unit      */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvald: set argument 1 */	
 setpvald(&eph,varstr,ffound,np,npok,pval,punit,ntxt);	
 if(ffound == 1)
    eph *= (ab[10]);


 /*====================== rain ratio ======================*/
 varstr = "RRAIN";      /* name of parameter     */
 npok   = 1;	        /* # par values expected */
 punit  = " (-)";       /* parameter unit        */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvald: set argument 1 */	
 setpvald(&rrain,varstr,ffound,np,npok,pval,punit,ntxt);	

	
 /*================= shelf/deep rain ======================*/	
 varstr = "FSHLF";      /* name of parameter     */
 npok   = 1;	        /* # par values expected */
 punit  = " (-)";       /* parameter unit        */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvald: set argument 1 */	
 setpvald(&fsh,varstr,ffound,np,npok,pval,punit,ntxt);

	
 /*================= CaCO3 riverine flux ==================*/	
 varstr = "FINC";       /* name of parameter     */
 npok   = 1;	        /* # par values expected */
 punit  = "  (mol C/y)   ";/* parameter unit     */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvald: set argument 1 */	
 setpvald(&finc0,varstr,ffound,np,npok,pval,punit,ntxt);
 if(ffound == 1)
    finc0 /= AOC;


 /*================= seawater [Ca2+] ======================*/	
 varstr = "CALC";      /* name of parameter     */
 npok   = 1;	       /* # par values expected */
 punit  = "  (mol/kg)    "; /* parameter unit   */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvald: set argument 1 */	
 setpvald(&cac,varstr,ffound,np,npok,pval,punit,ntxt);


 /*================= seawater [Mg2+] ======================*/	
 varstr = "MAGN";      /* name of parameter     */
 npok   = 1;	       /* # par values expected */
 punit  = "  (mol/kg)    "; /* parameter unit   */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvald: set argument 1 */	
 setpvald(&mgc,varstr,ffound,np,npok,pval,punit,ntxt);

 /*================= seawater [SO42-] =====================*/	
 varstr = "SULF";      /* name of parameter     */
 npok   = 1;	       /* # par values expected */
 punit  = "  (mol/kg)    "; /* parameter unit   */
 ntxt   = "none. Using default";              
 getparms(fparm,varstr,&ffound,&np,pval,pstr);
 /*setpvald: set argument 1 */	
 setpvald(&s4c,varstr,ffound,np,npok,pval,punit,ntxt);
	
 printf("\n@ Done reading parameter file: '%s'\n",fparmstr);
	
 /* close parameter file */
 fclose(fparm);

 free_dvector(pval,1,PVALMAX);
	
}
/*============================================================*/
/*==================== readparms() END========================*/
/*============================================================*/


/*============================================================*/
/*==================== getparms() ============================*/
/*============================================================*/
void getparms(FILE *fparm,char *target,int *ffound,int *np,
              double *pval,char *pstr)
{
 int fcmmnt,fcr;
 char *ptr, buff[BUFSIZ];
 char delims[]=" ",cmmnt[]="#",cr[]="\n";
 char *result=NULL;	

 *ffound = 0;
 *np = 0;
 rewind(fparm);
	
 /* find target */
 while(fgets(buff,(int)(BUFSIZ),fparm)) {
    ptr    = strstr(buff,target); 	  
	fcmmnt = strncmp(buff,cmmnt,1); 
  /*printf("     %d\n",strncmp(buff,cmmnt,1));*/
    if(ptr && fcmmnt != 0) {
    /*printf("\n target found\n");	
      printf("buff %s\n",buff);*/
      *ffound = 1;
	  break;
    }
 }
	
 if(*ffound == 1){	
  /* tokenize buff and convert results to string and
	 double. atof => double, see ref below            */
  result = strtok(buff,delims); /* ignore parm string */
  result = strtok(NULL,delims); 
  while(result != NULL){ 
  /*printf("result is '%s' \n",result);
    printf("pstr   is '%s' \n",pstr);
	printf("%d \n",strncmp(result,cr,1));*/
	fcr = strncmp(result,cr,1);
	if(fcr != 0){
       strcpy(pstr,result);
       removecr(cr[0],pstr);
       *np += 1;
       pval[*np] = atof(result);
	}
    result = strtok(NULL,delims);
    /*printf("pval %e\n",*pval);*/
  }
 }
	
}
/*============================================================*/
/*==================== getparms() END ========================*/
/*============================================================*/


/*============================================================*/
/*==================== setpstr() =============================*/
/*============================================================*/
/* arguments:

 *actflag action flag if ffound = 1                (out)
 fparstr  file parameter string                    (out)
 varstr   name of parameter                         (in)
 ffound   flag: found in control file 0/1           (in)
 np       # parameters found                        (in)
 npok     # parameters excpected                    (in)
 pstr     file parm string found in control file    (in)
 ntxt     display text if not found in control file (in)

 */
void setpstr(int *actflag,char *fparstr,char *varstr,int ffound,
             int np,int npok,char *pstr,char *ntxt)
{
 char mssg[BUFSIZ];

 if(ffound == 1){
   *actflag = 1;
   if(np != npok){
     sprintf(mssg,"setpstr(): %s requires %d file name(s). Check control file.",
     varstr,npok);
     ferrx(mssg);
   }
   strcpy(fparstr,pstr); /* SET GLOBAL VALUE */
   printf("%s %s\n",varstr,fparstr);
 } else {
     *actflag = 0;
     printf("%s: %s\n",varstr,ntxt);
 }
		
}	
/*============================================================*/
/*==================== setpstr() END =========================*/
/*============================================================*/

/*============================================================*/
/*==================== setpvald() ============================*/
/*============================================================*/
/* arguments:

 *parmd   parm value (double) if ffound = 1        (out)
 varstr   name of parameter                         (in)
 ffound   flag: found in control file 0/1           (in)
 np       # parameters found                        (in)
 npok     # parameters excpected                    (in)
 *pval    parm value found in control file          (in)
 punit    parameter unit                            (in)
 ntxt     display text if not found in control file (in)

 */
void setpvald(double *parmd,char *varstr,int ffound,int np,
              int npok,double *pval,char *punit,char *ntxt)
{
 char mssg[BUFSIZ];
 
 if(ffound == 1){
    if(np != npok){
      sprintf(mssg,"setpvald(): %s requires %d value(s). Check control file.",
              varstr,npok);
      ferrx(mssg);
	}
    *parmd = pval[1]; /* SET GLOBAL VALUE */
	printf("%s  %s %.2e\n",varstr,punit,*parmd);
 } else {
    printf("%s: %s \n",varstr,ntxt);
 }
		
}	
/*============================================================*/
/*==================== setpvald() END ========================*/
/*============================================================*/


/*============================================================*/
/*==================== setpvali() ============================*/
/*============================================================*/
/* arguments:

 *parmi   parm value (int) if ffound = 1           (out)
 varstr   name of parameter                         (in)
 ffound   flag: found in control file 0/1           (in)
 np       # parameters found                        (in)
 npok     # parameters excpected                    (in)
 *pval    parm value found in control file          (in)
 punit    parameter unit                            (in)
 ntxt     display text if not found in control file (in)

 */
void setpvali(int *parmi,char *varstr,int ffound,int np,
              int npok,double *pval,char *punit,char *ntxt)
{
 char mssg[BUFSIZ];

  
 if(ffound == 1){
    if(np != npok){
      sprintf(mssg,"setpvald(): %s requires %d value(s). Check control file.",
               varstr,npok);
      ferrx(mssg);
	}
    *parmi = (int)(pval[1]); /* SET GLOBAL VALUE */
	printf("%s  %s %d\n",varstr,punit,*parmi);
 } else {
    printf("%s: %s \n",varstr,ntxt);
 }
		
}	
/*============================================================*/
/*==================== setpvali() END ========================*/
/*============================================================*/

/*============================================================*/
/*==================== setpvaldv() ===========================*/
/*============================================================*/
/* arguments:

 *parmdv  parm value (double vec) if ffound = 1    (out)
 varstr   name of parameter                         (in)
 ffound   flag: found in control file 0/1           (in)
 np       # parameters found                        (in)
 npok     # parameters excpected                    (in)
 *pval    parm value found in control file          (in)
 punit    parameter unit                            (in)
 ntxt     display text if not found in control file (in)

 */
void setpvaldv(double *parmdv,char *varstr,int ffound,int np,
              int npok,double *pval,char *punit,char *ntxt)
{
 int n;	
 char mssg[BUFSIZ];

 /*
 printf(">> %s <<\n",varstr);
 printf(">> %d <<\n",ffound);
 printf(">> %d <<\n",np);
 printf(">> %d <<\n",npok);
 for(n=1;n<=np;n++)
     printf(">> %e <<",pval[n]);	
 printf(">> %s <<\n",punit);
 printf(">> %s <<\n",ntxt);
 */ 
  
 if(ffound == 1){
    if(np != npok){
      sprintf(mssg,"setpvaldv(): %s requires %d value(s). Check control file.",
              varstr,npok);
      ferrx(mssg);
	}
    printf("%s  %s ",varstr,punit);
    for(n=1;n<=npok;n++){
         parmdv[n] = pval[n]; /* SET GLOBAL VALUE */
         printf("%.2f ",parmdv[n]);
	}
    printf("\n");
 } else {
    printf("%s: %s \n",varstr,ntxt);
 }
		
}	
/*============================================================*/
/*==================== setpvaldv() END =======================*/
/*============================================================*/



/*============================================================*/
/*==================== removecr() ============================*/
/*============================================================*/
/* remove carriage return from string                         */
void removecr(char c, char *str)
{
 int i=0;
 int len = strlen(str)+1;

 for(i=0;i<len;i++)
 {
 if(str[i] == c)
    strncpy(&str[i],&str[i+1],len-i);
 }

}
/*============================================================*/
/*==================== removecr() END ========================*/
/*============================================================*/

/* 
 http://www.acm.uiuc.edu/webmonkeys/book/c_guide/
 http://www.elook.org/programming/c/strtok.html 
 */

 
