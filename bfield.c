
#include "bfield.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>

//GLOBALS SCOPED IN SOURCE FILE

double* xsurf;
double* ysurf;
double* zsurf;
double* nsurfx;
double* nsurfy;
double* nsurfz;

size_t size_surf;
double* currents;
int Nseg;
int Ncoils;

double* sfilx;
double* sfily;
double* sfilz;

int Nradfil;
int Ntorfil;

double* mfilx;
double* mfily;
double* mfilz;

double* Bsfilx;
double* Bsfily;
double* Bsfilz;
double* Bsfiln;

double* Bmfilx;
double* Bmfily;
double* Bmfilz;
double* Bmfiln;


void CalcSingleFilsB(void){
   double muc = 0.0000001;
   Bsfilx = (double*) malloc(size_surf*size_surf*sizeof(double));
   Bsfily = (double*) malloc(size_surf*size_surf*sizeof(double));
   Bsfilz = (double*) malloc(size_surf*size_surf*sizeof(double));
   Bsfiln = (double*) malloc(size_surf*size_surf*sizeof(double));
 


   double ex, ey, ez, norm, ri, rix, riy, riz, rf, rfx, rfy, rfz, l;
   int i,j,k;
  
   for(i=0;i<size_surf*size_surf;i++){
      for(j=0;j<Ncoils;j++){
         //Before looping over all segments, calculate B from first seg to improve speed of loop
         l = sqrt((*(sfilx+j*Nseg+1)-*(sfilx+j*Nseg))*(*(sfilx+j*Nseg+1)-*(sfilx+j*Nseg)) + (*(sfily+j*Nseg+1)-*(sfily+j*Nseg))*(*(sfily+j*Nseg+1)-*(sfily+j*Nseg)) + (*(sfilz+j*Nseg+1)-*(sfilz+j*Nseg))*(*(sfilz+j*Nseg+1)-*(sfilz+j*Nseg)));  
         ex = (*(sfilx+j*Nseg+1) - *(sfilx+j*Nseg))/l;
         ey = (*(sfily+j*Nseg+1) - *(sfily+j*Nseg))/l;
         ez = (*(sfilz+j*Nseg+1) - *(sfilz+j*Nseg))/l;
 
         rix = *(xsurf+i) - *(sfilx+j*Nseg);
         riy = *(ysurf+i) - *(sfily+j*Nseg); 
         riz = *(zsurf+i) - *(sfilz+j*Nseg); 
         ri = sqrt(rix*rix + riy*riy + riz*riz);
         
         rfx = *(xsurf+i) - *(sfilx+j*Nseg+1);
         rfy = *(ysurf+i) - *(sfily+j*Nseg+1); 
         rfz = *(zsurf+i) - *(sfilz+j*Nseg+1); 
         rf = sqrt(rfx*rfx + rfy*rfy + rfz*rfz);
         
         *(Bsfilx+i) = muc*(*(currents+i))*(ey*riz-ez*riy)*(2*l*(ri+rf)/(ri*rf))*(1/((ri+rf)*(ri+rf)-l*l));
         *(Bsfily+i) = muc*(*(currents+i))*(ez*rix-ex*riz)*(2*l*(ri+rf)/(ri*rf))*(1/((ri+rf)*(ri+rf)-l*l));
         *(Bsfilz+i) = muc*(*(currents+i))*(ex*riy-ey*rix)*(2*l*(ri+rf)/(ri*rf))*(1/((ri+rf)*(ri+rf)-l*l));
  
         rix = rfx;
         riy = rfy;
         riz = rfz;
         ri = rf;
     
         for(k=1;k<Nseg-1;k++){
            l = sqrt((*(sfilx+j*Nseg+k+1)-*(sfilx+j*Nseg+k))*(*(sfilx+j*Nseg+k+1)-*(sfilx+j*Nseg+k)) + (*(sfily+j*Nseg+k+1)-*(sfily+j*Nseg+k))*(*(sfily+j*Nseg+k+1)-*(sfily+j*Nseg+k)) + (*(sfilz+j*Nseg+k+1)-*(sfilz+j*Nseg+k))*(*(sfilz+j*Nseg+k+1)-*(sfilz+j*Nseg+k)));  
           
            ex = (*(sfilx+j*Nseg+k+1) - *(sfilx+j*Nseg+k))/l;
            ey = (*(sfily+j*Nseg+k+1) - *(sfily+j*Nseg+k))/l;
            ez = (*(sfilz+j*Nseg+k+1) - *(sfilz+j*Nseg+k))/l;
              
            rfx = *(xsurf+i) - *(sfilx+j*Nseg+k+1);
            rfy = *(ysurf+i) - *(sfily+j*Nseg+k+1); 
            rfz = *(zsurf+i) - *(sfilz+j*Nseg+k+1);           
            rf = sqrt(rfx*rfx + rfy*rfy + rfz*rfz);
      
            *(Bsfilx+i) = *(Bsfilx+i) + muc*(*(currents+i))*(ey*riz-ez*riy)*(2*l*(ri+rf)/(ri*rf))*(1/((ri+rf)*(ri+rf)-l*l));
            *(Bsfily+i) = *(Bsfily+i) + muc*(*(currents+i))*(ez*rix-ex*riz)*(2*l*(ri+rf)/(ri*rf))*(1/((ri+rf)*(ri+rf)-l*l));
            *(Bsfilx+i) = *(Bsfilz+i) + muc*(*(currents+i))*(ex*riy-ey*rix)*(2*l*(ri+rf)/(ri*rf))*(1/((ri+rf)*(ri+rf)-l*l));
      
           rix = rfx;
           riy = rfy;
           riz = rfz;
           ri = rf;
         }
      }

      l = sqrt((*(sfilx+j*Nseg+k+1)-*(sfilx+j*Nseg+k))*(*(sfilx+j*Nseg+k+1)-*(sfilx+j*Nseg+k)) + (*(sfily+j*Nseg+k+1)-*(sfily+j*Nseg+k))*(*(sfily+j*Nseg+k+1)-*(sfily+j*Nseg+k)) + (*(sfilz+j*Nseg+k+1)-*(sfilz+j*Nseg+k))*(*(sfilz+j*Nseg+k+1)-*(sfilz+j*Nseg+k)));  
           
      ex = (*(sfilx+j*Nseg) - *(sfilx+j*Nseg+Nseg-1))/l;
      ey = (*(sfily+j*Nseg) - *(sfily+j*Nseg+Nseg-1))/l;
      ez = (*(sfilz+j*Nseg) - *(sfilz+j*Nseg+Nseg-1))/l;
              
      rfx = *(xsurf+i) - *(sfilx+j*Nseg+Nseg-1);
      rfy = *(ysurf+i) - *(sfily+j*Nseg+Nseg-1); 
      rfz = *(zsurf+i) - *(sfilz+j*Nseg+Nseg-1);           
      rf = sqrt(rfx*rfx + rfy*rfy + rfz*rfz);
      
      *(Bsfilx+i) = *(Bsfilx+i) + muc*(*(currents+i))*(ey*riz-ez*riy)*(2*l*(ri+rf)/(ri*rf))*(1/((ri+rf)*(ri+rf)-l*l));
      *(Bsfily+i) = *(Bsfily+i) + muc*(*(currents+i))*(ez*rix-ex*riz)*(2*l*(ri+rf)/(ri*rf))*(1/((ri+rf)*(ri+rf)-l*l));
      *(Bsfilz+i) = *(Bsfilz+i) + muc*(*(currents+i))*(ex*riy-ey*rix)*(2*l*(ri+rf)/(ri*rf))*(1/((ri+rf)*(ri+rf)-l*l));
      *(Bsfiln+i) = *(Bsfilx+i)* *(nsurfx+i) + *(Bsfily+i)* *(nsurfy+i) + *(Bsfilz+i)* *(nsurfz+i);  
      double B = sqrt(*(Bsfilx+i)* *(Bsfilx+i) + *(Bsfily+i)* *(Bsfily+i) + *(Bsfilz+i)* *(Bsfilz+i)); 
      printf("%.15f %.15f %.15f %.8f\n", *(xsurf+i),*(ysurf+i),*(zsurf+i),B);   
  }

}


//void CalcMultiFilsB(void){}
//void DiagnoseSingleFilsB(void){}
//void DiagnoseMultiFilsB(void){}


