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
//---------------------------------------------------------------------------------------------------- 
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

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
         
void CalculateMultiField(double x, double y, double z, \
                          double* Bx, double* By, double* Bz){
//----------------------------------------------------------------------------------------------------
// Calculates the field due to multi-filament representation at the point (x,y,z)
// Field components are implicitly stored in (Bx,By,Bz) in function call
// Based on the method of Hanson and Hirschman
// Doesn't enforce periodicity or stellarator symmetry
//----------------------------------------------------------------------------------------------------
   double muc = 1.0e-7;
   int i,j,k;
   int Nfils = Nradfil*Ntorfil; 
   double l, ex, ey, ez, bx, by, bz;
   double ri, rix, riy, riz;
   double rf, rfx, rfy, rfz;
   
   bx = 0.0;
   by = 0.0;
   bz = 0.0;
  
   omp_set_num_threads(Nthreads);
   #pragma omp parallel for  
   for(i=0;i<Ncoils;i++){
      for(j=0;j<Nfils;j++){ //index over each filament comprising the multi-filament coil
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
   int iCoils = Ncoils / Nfp;
   
   double muc = 1.0e-7;
   double one = 1.00000000000000000000;
   double two = 2.00000000000000000000;
   register int ip,i,j,k;
   int Nfils = Nradfil*Ntorfil; 
   double factor = muc * two / Nfils;
   double l, ex, ey, ez, bx, by, bz, bxx, byy, bzz;
   double xx, yy, zz;
   double ri, xi, yi, zi;
   double rf, xf, yf, zf;
   double cur, eps, eta;

   bx = 0.0;
   by = 0.0;
   bz = 0.0;

   bxx = 0.0;
   byy = 0.0;
   bzz = 0.0;

   for(ip=1;ip<Nfp+1;ip++)
   { 
      //Find symmetric point on other field periods
      xx =  x*cosnfp(ip) - y*sinnfp(ip); 
      yy =  x*sinnfp(ip) + y*cosnfp(ip); 
      zz =  z; 
      for(i=0;i<iCoils;i++)
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

            for(k=0;k<Nseg;k++) //TODO: this assumes that we have even nseg!
            {
               xf = xx - *(mfilx+i*Nfils*(Nseg+1)+j*(Nseg+1)+ k + 1);
               yf = yy - *(mfily+i*Nfils*(Nseg+1)+j*(Nseg+1)+ k + 1);
               zf = zz - *(mfilz+i*Nfils*(Nseg+1)+j*(Nseg+1)+ k + 1);
               rf = sqrt( xf*xf + yf*yf + zf*zf);
               
               l = sqrt( pow((xf-xi),2) + pow((yf-yi),2) + pow((zf-zi),2));  

               ex = ( xf - xi ) / l;  
               ey = ( yf - yi ) / l;  
               ez = ( zf - zi ) / l;  
       
               eps = l / (ri + rf);
               eta = ri*rf;

               bx = cur * (ey*zi-ez*yi) * eps / (eta * (one - eps * eps));
               by = cur * (ez*xi-ex*zi) * eps / (eta * (one - eps * eps));
               bz = cur * (ex*yi-ey*xi) * eps / (eta * (one - eps * eps));

               //Find equivalent contribution to field period 1
               bxx +=  bx*cosnfp(ip) + by*sinnfp(ip);  
               byy += -bx*sinnfp(ip) + by*cosnfp(ip); 
               bzz +=  bz;        
               
               //End of segment k becomes beginning of segment k+1
               xi = xf;
               yi = yf;
               zi = zf;
               ri = rf;
            }
         }
      }
   }
   *Bx = bxx * factor;
   *By = byy * factor;
   *Bz = bzz * factor;
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
/*
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
   int iCoils = Ncoils / Nfp;
   
   double muc = 1.0e-7;
   register int ip,i,j,k;
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

   //double* filx = mfilx;
   //double* fily = mfily;
   //double* filz = mfilz;
   //double* tcurrents = currents;
   //omp_set_num_threads(Nthreads);
   //#pragma omp parallel for  

   for(ip=1;ip<Nfp+1;ip++){ //ip is used in cosnfp and sinnfp for finding the contribution to field at ip-th field period   
      xx =  x*cosnfp(ip) - y*sinnfp(ip); //find the x component of the periodic point at ip+1-th field period
      yy =  x*sinnfp(ip) + y*cosnfp(ip); //find the y component of the periodic point at ip+1-th field period
      zz =  z; //no change in the z component
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
        
                bx = muc * ( *(currents+i) / Nfils ) * (ey*riz-ez*riy) * (l*2*(ri+rf)) \
                     / ( (ri*rf) * ( pow((ri+rf),2) - pow(l,2)) );

                by = muc * ( *(currents+i) / Nfils ) * (ez*rix-ex*riz) * (l*2*(ri+rf)) \
                     / ( (ri*rf) * ( pow((ri+rf),2) - pow(l,2)) );

                bz = muc * ( *(currents+i) / Nfils ) * (ex*riy-ey*rix) * (l*2*(ri+rf)) \
                     / ( (ri*rf) * ( pow((ri+rf),2) - pow(l,2))  );
        
    
                bxx += bx*cosnfp(ip) + by*sinnfp(ip);  //rotate bx back to the 1st field period
                byy += -bx*sinnfp(ip) + by*cosnfp(ip); //rotate by back to the 1st field period
                bzz += bz;        
            }
         }
      }
   }
   *Bx = bxx;
   *By = byy;
   *Bz = bzz;
}
*/

