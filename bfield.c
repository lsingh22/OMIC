
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

   double ex, ey, ez, norm, ri, rix, riy, riz, rf, rfx, rfy, rfz, l,B;
   int i,j,k;
 
   //for(i=0;i<size_surf*size_surf;i++){
   for(i=0;i<10;i++){
      for(j=0;j<Ncoils;j++){
         //Before looping over all segments, calculate B from first seg to improve speed of loop
         l = sqrt( pow((*(sfilx+j*Nseg+1)-*(sfilx+j*Nseg)),2) + pow((*(sfily+j*Nseg+1)-*(sfily+j*Nseg)),2)) + pow((*(sfilz+j*Nseg+1)-*(sfilz+j*Nseg)),2);  
   // fprintf(fb, "periods 1\n begin filament\n mirror NIL\n");
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
 
         *(Bsfilx+i) = 0.0;//muc * *(currents+j)*(ey*riz-ez*riy)*(2*l*(ri+rf)/(ri*rf))*(1/((ri+rf)*(ri+rf)-l*l));
         *(Bsfily+i) = 0.0;//muc * *(currents+j)*(ez*rix-ex*riz)*(2*l*(ri+rf)/(ri*rf))*(1/((ri+rf)*(ri+rf)-l*l));
         *(Bsfilz+i) = 0.0;//muc * *(currents+j)*(ex*riy-ey*rix)*(2*l*(ri+rf)/(ri*rf))*(1/((ri+rf)*(ri+rf)-l*l));
         //B = sqrt(*(Bsfilx+i)* *(Bsfilx+i) + *(Bsfily+i)* *(Bsfily+i) + *(Bsfilz+i)* *(Bsfilz+i)); 
      
        //printf("%.15f %.15f %.15f %.15f\n", *(xsurf+i),*(ysurf+i),*(zsurf+i),B);   
        //printf("%.15f %.15f %.15f %.15f\n", *(Bsfilx+i),*(Bsfily+i),*(Bsfilz+i),B);   
        //printf("%.15f %.15f %.15f\n", *(Bsfilx+i),*(currents+j) , muc*(*(currents+j))*(ey*riz-ez*riy)*(2*l*(ri+rf)/(ri*rf))*(1/((ri+rf)*(ri+rf)-l*l)));
 
    
         rix = rfx;
         riy = rfy;
         riz = rfz;
         ri = rf;
     
         for(k=1;k<Nseg-1;k++){

            l = sqrt( pow((*(sfilx+j*Nseg+k+1)-*(sfilx+j*Nseg+k)),2) ) + pow((*(sfily+j*Nseg+k+1)-*(sfily+j*Nseg+k)),2) + pow((*(sfilz+j*Nseg+k+1)-*(sfilz+j*Nseg+k)),2);  
           
            ex = (*(sfilx+j*Nseg+k+1) - *(sfilx+j*Nseg+k))/l;
            ey = (*(sfily+j*Nseg+k+1) - *(sfily+j*Nseg+k))/l;
            ez = (*(sfilz+j*Nseg+k+1) - *(sfilz+j*Nseg+k))/l;
             
            rix = *(xsurf+i) - *(sfilx+j*Nseg+k);
            riy = *(ysurf+i) - *(sfily+j*Nseg+k); 
            riz = *(zsurf+i) - *(sfilz+j*Nseg+k);           
            ri = sqrt(rix*rix + riy*riy + riz*riz); 

            rfx = *(xsurf+i) - *(sfilx+j*Nseg+k+1);
            rfy = *(ysurf+i) - *(sfily+j*Nseg+k+1); 
            rfz = *(zsurf+i) - *(sfilz+j*Nseg+k+1);           
            rf = sqrt(rfx*rfx + rfy*rfy + rfz*rfz);

 
            *(Bsfilx+i) += muc* *(currents+j) *(ey*riz-ez*riy)*(2*l*(ri+rf)/(ri*rf))*(1/((ri+rf)*(ri+rf)-l*l));
            *(Bsfily+i) += muc* *(currents+j) *(ez*rix-ex*riz)*(2*l*(ri+rf)/(ri*rf))*(1/((ri+rf)*(ri+rf)-l*l));
            *(Bsfilz+i) += muc* *(currents+j) *(ex*riy-ey*rix)*(2*l*(ri+rf)/(ri*rf))*(1/((ri+rf)*(ri+rf)-l*l));
            //rix = rfx;
            //riy = rfy;
            //riz = rfz;
            //ri = rf;
            
         }
   
         l = sqrt((*(sfilx+j*Nseg)-*(sfilx+j*Nseg+Nseg-1))*(*(sfilx+j*Nseg)-*(sfilx+j*Nseg+Nseg-1)) + (*(sfily+j*Nseg)-*(sfily+j*Nseg+Nseg-1))*(*(sfily+j*Nseg)-*(sfily+j*Nseg+Nseg-1)) + (*(sfilz+j*Nseg)-*(sfilz+j*Nseg+Nseg-1))*(*(sfilz+j*Nseg)-*(sfilz+j*Nseg+Nseg-1)));  
           
         ex = (*(sfilx+j*Nseg) - *(sfilx+j*Nseg+Nseg-1))/l;
         ey = (*(sfily+j*Nseg) - *(sfily+j*Nseg+Nseg-1))/l;
         ez = (*(sfilz+j*Nseg) - *(sfilz+j*Nseg+Nseg-1))/l;
              
         rfx = *(xsurf+i) - *(sfilx+j*Nseg+Nseg-1);
         rfy = *(ysurf+i) - *(sfily+j*Nseg+Nseg-1); 
         rfz = *(zsurf+i) - *(sfilz+j*Nseg+Nseg-1);           
         rf = sqrt(rfx*rfx + rfy*rfy + rfz*rfz);
      
         //*(Bsfilx+i) = *(Bsfilx+i) + muc*(*(currents+j))*(ey*riz-ez*riy)*(2*l*(ri+rf)/(ri*rf))*(1/((ri+rf)*(ri+rf)-l*l));
         //*(Bsfily+i) = *(Bsfily+i) + muc*(*(currents+j))*(ez*rix-ex*riz)*(2*l*(ri+rf)/(ri*rf))*(1/((ri+rf)*(ri+rf)-l*l));
         //*(Bsfilz+i) = *(Bsfilz+i) + muc*(*(currents+j))*(ex*riy-ey*rix)*(2*l*(ri+rf)/(ri*rf))*(1/((ri+rf)*(ri+rf)-l*l));
         //*(Bsfiln+i) = *(Bsfilx+i)* *(nsurfx+i) + *(Bsfily+i)* *(nsurfy+i) + *(Bsfilz+i)* *(nsurfz+i);  

         //printf("%.15f %.15f %.15f %.15f\n", *(xsurf+i),*(ysurf+i),*(zsurf+i),B);   

 
      }
      B = sqrt(*(Bsfilx+i) * *(Bsfilx+i) + *(Bsfily+i) * *(Bsfily+i) + *(Bsfilz+i) * *(Bsfilz+i)); 
      printf("%.15f %.15f %.15f %.15f\n", *(Bsfilx+i),*(Bsfily+i),*(Bsfilz+i),B);   
      
   }
}


//void CalcMultiFilsB(void){}
//void DiagnoseSingleFilsB(void){}
//void DiagnoseMultiFilsB(void){}



