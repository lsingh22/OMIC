
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

void CalculateSingleField(double x, double y, double z, \
                          double* Bx, double* By, double* Bz){
  
   double muc = 1.0e-7;
   int i,j;
 
   double l, ex, ey, ez, bx, by, bz;
   double ri, rix, riy, riz;
   double rf, rfx, rfy, rfz;
   
   bx = 0.0;
   by = 0.0;
   bz = 0.0;
   
   for(i=0;i<Ncoils;i++){
      for(j=0;j<(Nseg+1);j++){

         l = sqrt( pow( ((*(sfilx+i*(Nseg+1)+j+1)) - (*(sfilx+i*(Nseg+1)+j))) ,2) + \
                   pow( ((*(sfily+i*(Nseg+1)+j+1)) - (*(sfily+i*(Nseg+1)+j))) ,2) + \
                   pow( ((*(sfilz+i*(Nseg+1)+j+1)) - (*(sfilz+i*(Nseg+1)+j))) ,2) ); 
      
         ex = ( (*(sfilx+i*(Nseg+1)+j+1)) - (*(sfilx+i*(Nseg+1)+j)) ) / l;  
         ey = ( (*(sfily+i*(Nseg+1)+j+1)) - (*(sfily+i*(Nseg+1)+j)) ) / l;    
         ez = ( (*(sfilz+i*(Nseg+1)+j+1)) - (*(sfilz+i*(Nseg+1)+j)) ) / l;   


         rix = x - (*(sfilx+i*(Nseg+1)+j));
         riy = y - (*(sfily+i*(Nseg+1)+j));
         riz = z - (*(sfilz+i*(Nseg+1)+j));
          ri = sqrt( pow(rix,2) + pow(riy,2) + pow(riz,2) );


         rfx = x - (*(sfilx+i*(Nseg+1)+j+1)); 
         rfy = y - (*(sfily+i*(Nseg+1)+j+1));  
         rfz = z - (*(sfilz+i*(Nseg+1)+j+1));   
          rf = sqrt( pow(rfx,2) + pow(rfy,2) + pow(rfz,2) );
 
         
          bx += muc * (*(currents+i)) * (ey*riz-ez*riy) * (2*l*(ri+rf)) \
               / ( (ri*rf) * ( pow((ri+rf),2) - pow(l,2)) );

          by += muc * (*(currents+i)) * (ez*rix-ex*riz) * (2*l*(ri+rf)) \
               / ( (ri*rf) * ( pow((ri+rf),2) - pow(l,2)) );

          bz += muc * (*(currents+i)) * (ex*riy-ey*rix) * (2*l*(ri+rf)) \
               / ( (ri*rf) * ( pow((ri+rf),2) - pow(l,2))  );
      }
   }

   *Bx = bx;
   *By = by;
   *Bz = bz;
}


         
void CalcMultiFilsB(double x, double y, double z, \
                          double* Bx, double* By, double* Bz){


  int i,j,k;
  int Nfils; 

   for(i=0;i<Ncoils;i++){
      for(j=0;j<Nfils;j++){
         for(k=0;k<Nseg+1;k++){
            
         }
      }
   }


}
//void DiagnoseSingleFilsB(void){}
//void DiagnoseMultiFilsB(void){}



