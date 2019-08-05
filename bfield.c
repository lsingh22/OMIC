
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

int Nzeta;
int Nteta;
int Nfp;

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

double cosnfp(int ip){

   double val;
   double twopi = M_PI*2;
   val = cos( (ip-1)*twopi / Nfp);

return val;
}

double sinnfp(int ip){

   double val;
   double twopi = M_PI*2;
   val = sin( (ip-1)*twopi / Nfp);

return val;
}

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
      for(j=0;j<Nseg;j++){

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


         
void CalculateMultiField(double x, double y, double z, \
                          double* Bx, double* By, double* Bz){

   double muc = 1.0e-7;
   int i,j,k;
   int Nfils = Nradfil*Ntorfil; 
   double l, ex, ey, ez, bx, by, bz;
   double ri, rix, riy, riz;
   double rf, rfx, rfy, rfz;
   
   bx = 0.0;
   by = 0.0;
   bz = 0.0;
   
   for(i=0;i<Ncoils;i++){
      for(j=0;j<Nfils;j++){
         for(k=0;k<Nseg;k++){
            l = sqrt( pow( ((*(mfilx+i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1)) - (*(mfilx+i*Nfils*(Nseg+1)+j*(Nseg+1)+k)) ) ,2) + \
                      pow( ((*(mfily+i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1)) - (*(mfily+i*Nfils*(Nseg+1)+j*(Nseg+1)+k)) ) ,2) + \
                      pow( ((*(mfilz+i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1)) - (*(mfilz+i*Nfils*(Nseg+1)+j*(Nseg+1)+k)) ) ,2) ); 
               
            ex = ( (*(mfilx+i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1) ) - (*(mfilx+i*Nfils*(Nseg+1)+j*(Nseg+1)+k)) ) / l;  
            ey = ( (*(mfily+i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1) ) - (*(mfily+i*Nfils*(Nseg+1)+j*(Nseg+1)+k)) ) / l;    
            ez = ( (*(mfilz+i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1) ) - (*(mfilz+i*Nfils*(Nseg+1)+j*(Nseg+1)+k)) ) / l;   


            rix = x - ( *(mfilx+i*Nfils*(Nseg+1)+j*(Nseg+1)+k) );
            riy = y - ( *(mfily+i*Nfils*(Nseg+1)+j*(Nseg+1)+k) );
            riz = z - ( *(mfilz+i*Nfils*(Nseg+1)+j*(Nseg+1)+k) );
             ri = sqrt( pow(rix,2) + pow(riy,2) + pow(riz,2) );


            rfx = x - ( *(mfilx+i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1) ); 
            rfy = y - ( *(mfily+i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1) );  
            rfz = z - ( *(mfilz+i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1) );   
             rf = sqrt( pow(rfx,2) + pow(rfy,2) + pow(rfz,2) );
 
         
             bx += muc * ( *(currents+i) / Nfils ) * (ey*riz-ez*riy) * (2*l*(ri+rf)) \
                  / ( (ri*rf) * ( pow((ri+rf),2) - pow(l,2)) );

             by += muc * ( *(currents+i) / Nfils ) * (ez*rix-ex*riz) * (2*l*(ri+rf)) \
                  / ( (ri*rf) * ( pow((ri+rf),2) - pow(l,2)) );

             bz += muc * ( *(currents+i) / Nfils ) * (ex*riy-ey*rix) * (2*l*(ri+rf)) \
                  / ( (ri*rf) * ( pow((ri+rf),2) - pow(l,2))  );
            
            
         }
      }
   }
   
   *Bx = bx;
   *By = by;
   *Bz = bz;

}


void CalculateMultiFieldSym(double x, double y, double z, \
                          double* Bx, double* By, double* Bz){
   
   int iCoils = Ncoils / Nfp;
   
   double muc = 1.0e-7;
   int ip,i,j,k;
   int Nfils = Nradfil*Ntorfil; 
   double l, ex, ey, ez, bx, by, bz, bxx, byy, bzz;
   double xx, yy, zz;
   double ri, rix, riy, riz;
   double rf, rfx, rfy, rfz;
   
   bx = 0.0;
   by = 0.0;
   bz = 0.0;

   bxx = 0.0;
   byy = 0.0;
   bzz = 0.0;

   
   // Symmetry calculation taken straight from FOCUS
   for(ip=1;ip<Nfp+1;ip++){   
      //Find periodic point
      xx =  x*cosnfp(ip) + y*sinnfp(ip);
      yy = -x*sinnfp(ip) + y*cosnfp(ip);
      zz =  z;
      for(i=0;i<iCoils;i++){
         for(j=0;j<Nfils;j++){
            for(k=0;k<Nseg;k++){
               l = sqrt( pow( ((*(mfilx+i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1)) - (*(mfilx+i*Nfils*(Nseg+1)+j*(Nseg+1)+k)) ) ,2) + \
                      pow( ((*(mfily+i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1)) - (*(mfily+i*Nfils*(Nseg+1)+j*(Nseg+1)+k)) ) ,2) + \
                      pow( ((*(mfilz+i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1)) - (*(mfilz+i*Nfils*(Nseg+1)+j*(Nseg+1)+k)) ) ,2) ); 
               
               ex = ( (*(mfilx+i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1) ) - (*(mfilx+i*Nfils*(Nseg+1)+j*(Nseg+1)+k)) ) / l;  
               ey = ( (*(mfily+i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1) ) - (*(mfily+i*Nfils*(Nseg+1)+j*(Nseg+1)+k)) ) / l;    
               ez = ( (*(mfilz+i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1) ) - (*(mfilz+i*Nfils*(Nseg+1)+j*(Nseg+1)+k)) ) / l;   


               rix = xx - ( *(mfilx+i*Nfils*(Nseg+1)+j*(Nseg+1)+k) );
               riy = yy - ( *(mfily+i*Nfils*(Nseg+1)+j*(Nseg+1)+k) );
               riz = zz - ( *(mfilz+i*Nfils*(Nseg+1)+j*(Nseg+1)+k) );
                ri = sqrt( pow(rix,2) + pow(riy,2) + pow(riz,2) );


               rfx = xx - ( *(mfilx+i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1) ); 
               rfy = yy - ( *(mfily+i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1) );  
               rfz = zz - ( *(mfilz+i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1) );   
                rf = sqrt( pow(rfx,2) + pow(rfy,2) + pow(rfz,2) );
 
         
                bx = muc * ( *(currents+i) / Nfils ) * (ey*riz-ez*riy) * (2*l*(ri+rf)) \
                     / ( (ri*rf) * ( pow((ri+rf),2) - pow(l,2)) );

                by = muc * ( *(currents+i) / Nfils ) * (ez*rix-ex*riz) * (2*l*(ri+rf)) \
                     / ( (ri*rf) * ( pow((ri+rf),2) - pow(l,2)) );

                bz = muc * ( *(currents+i) / Nfils ) * (ex*riy-ey*rix) * (2*l*(ri+rf)) \
                     / ( (ri*rf) * ( pow((ri+rf),2) - pow(l,2))  );
            

                bxx += bx*cosnfp(ip) - by*sinnfp(ip);
                byy += bx*sinnfp(ip) + by*cosnfp(ip);
                bzz += bz;
            
            }
         }
      }
   }

   *Bx = bxx;
   *By = byy;
   *Bz = bzz;

}

