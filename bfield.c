#include "bfield.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <omp.h>
//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

//GLOBALS SCOPED IN SOURCE FILE
//See globals.h for description of parameters and arrays

double* xsurf;
double* ysurf;
double* zsurf;

double* nsurfx;
double* nsurfy;
double* nsurfz;

int Ncoil;
int iCoil;
int Nzeta;
int Nteta;
int Nfp;

double* currents;
int Nseg;

double* sfilx;
double* sfily;
double* sfilz;

int Nradfil;
int Ntorfil;
int Nfils;

double* mfilx;
double* mfily;
double* mfilz;

int Ns;
//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

double cosnfp(int ip){ 
//----------------------------------------------------------------------------------------------------
// Returns the cosine of a vector component at the ip-th field period
//----------------------------------------------------------------------------------------------------
   double val;
   double twopi = M_PI*2;
   val = cos( (ip-1)*twopi / Nfp);
   return val;
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

double sinnfp(int ip){ 
//----------------------------------------------------------------------------------------------------
// Returns the sine of a vector component at the ip-th field period
//----------------------------------------------------------------------------------------------------
   double val;
   double twopi = M_PI*2;
   val = sin( (ip-1)*twopi / Nfp);
   return val;
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

void CalculateSingleField(double x, double y, double z, \
                          double* Bx, double* By, double* Bz){
//----------------------------------------------------------------------------------------------------
// Calculates the field due to single filament representation at the point (x,y,z)
// Field components are implicitly stored in (Bx,By,Bz) in function call
// Based on the method of Hanson and Hirschman
// See publication for parameter definitions (TODO:insert url here)
// TODO: this should be exactly the same as MultiFilFieldSym except for indexing
//---------------------------------------------------------------------------------------------------- 
   double muc = 1.0e-7;
   int i,j;
   
   double l, ex, ey, ez, bx, by, bz;
   double ri, rix, riy, riz;
   double rf, rfx, rfy, rfz;
   
   bx = 0.0;
   by = 0.0;
   bz = 0.0;
   
   for(i=0;i<Ncoil;i++){
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

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

void CalculateMultiFieldSym(double x, double y, double z, \
                          double* Bx, double* By, double* Bz){
//----------------------------------------------------------------------------------------------------
// Calculates the field due to multi-filament representation at the point (x,y,z)
// Field components are implicitly stored in (Bx,By,Bz) in function call
// Based on the method of Hanson and Hirschman
// Periodicity is supported, but not stellarator symmetry
//
// Symmetry calculation taken straight from FOCUS:
// The input point is used to find the same point on each period,
// (x,y,z) --> (xxi_1,yy_1,zz_1) ... (xx_nfp-1,yy_nfp-1,zz_nfp-1)
// The field is then calculated at each of the nfp points, and each contribution is superposed upon proper
// rotation back to the initial field period. 
// Rotations are 'swept under the rug' by the cosnfp,sinnfp functions. 
//----------------------------------------------------------------------------------------------------  
   
   double muc = 1.0e-7;
   double one = 1.00000000000000000000; //these should probably be global
   double two = 2.00000000000000000000;
   register int ip,i,j,k;
   double factor = muc * two / Nfils;
   double l, ex, ey, ez, bx, by, bz, bxx, byy, bzz;
   double xx, yy, zz;
   double ri, xi, yi, zi;
   double rf, xf, yf, zf;
   double cur, eps, eta, coef;
   double rot_cos, rot_sin;

   bx = 0.0;
   by = 0.0;
   bz = 0.0;

   bxx = 0.0;
   byy = 0.0;
   bzz = 0.0;
   zz = z;
   
   for(ip=1;ip<Nfp+1;ip++)
   { 
      //Find coefficients of transformation matrix
      rot_cos = cosnfp(ip);
      rot_sin = sinnfp(ip);

      //Find symmetric point on other field periods
      xx =  x * rot_cos - y * rot_sin; 
      yy =  x * rot_sin + y * rot_cos; 
      for(i=0;i<iCoil;i++)
      {
         //Store current of i-th coil 
         cur = *(currents+i);
         for(j=0;j<Nfils;j++)
         {
            //Store first point of current filament before main loop
            xi = xx - *(mfilx+i*Nfils*(Nseg+1)+j*(Nseg+1));
            yi = yy - *(mfily+i*Nfils*(Nseg+1)+j*(Nseg+1));
            zi = zz - *(mfilz+i*Nfils*(Nseg+1)+j*(Nseg+1));
            ri = sqrt( xi*xi + yi*yi + zi*zi );           
            
            for(k=0;k<Nseg;k++) //TODO: check if this works for odd nseg
            {
               xf = xx - *(mfilx+i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1);
               yf = yy - *(mfily+i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1);
               zf = zz - *(mfilz+i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1);
               rf = sqrt( xf*xf + yf*yf + zf*zf);
               
               l = sqrt( (xf-xi)*(xf-xi) + (yf-yi)*(yf-yi) + (zf-zi)*(zf-zi));  

               ex = ( xf - xi ) / l;  
               ey = ( yf - yi ) / l;  
               ez = ( zf - zi ) / l;  
       
               eps = l / (ri + rf);
               eta = ri * rf;
               coef = cur * eps / (eta * (one - eps * eps)); 
   
               bx += coef * (ey*zi-ez*yi);
               by += coef * (ez*xi-ex*zi);
               bz += coef * (ex*yi-ey*xi);
              
               //End of segment k becomes beginning of segment k+1
               xi = xf;
               yi = yf;
               zi = zf;
               ri = rf;
            }
         }
      }    
      //Rotate back to first field period
      bxx +=  bx * rot_cos + by * rot_sin;  
      byy += -bx * rot_sin + by * rot_cos; 
      bzz +=  bz;

      bx = 0;
      by = 0;
      bz = 0;          
   }
   *Bx = bxx * factor;
   *By = byy * factor;
   *Bz = bzz * factor;
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

