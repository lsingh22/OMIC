#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <omp.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "bfield.h"
#include "bfield_gpu.cuh"

#define MY_PI 3.14159265358979323846

double cosnfp(int ip) { 

   // Cosine of vector component at the ip-th field period
   return cos((ip - 1) * MY_PI * 2 / Nfp);
}

double sinnfp(int ip) { 

   // Sine of vector component at the ip-th field period
   return sin((ip - 1) * MY_PI * 2 / Nfp);
}

void CalculateSingleField(double x, double y, double z, \
                          double* Bx, double* By, double* Bz){
//----------------------------------------------------------------------------------------------------
// Calculates the field due to single filament representation at the point (x,y,z)
// Based on the method of Hanson and Hirschman
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

void CalculateFieldSerial(int start, int end) {

   register int i;
	
   // Set the number of threads to the max available
   omp_set_num_threads(Nthreads);   

   // Parallel thread region to chunk boundary points
   #pragma omp parallel for
   for(i = start; i < end; i++) {      
		CalculateFieldAtPoint(xsurf[i], ysurf[i], zsurf[i], Bmfilx+i, Bmfily+i, Bmfilz+i);    
	}
}

void CalculateFieldAtPoint(double x, double y, double z, \
                           double* Bx, double* By, double* Bz){
//----------------------------------------------------------------------------------------------------
// Calculates the field due to multi-filament representation at the point (x,y,z)
// Based on the method of Hanson and Hirschman
// Periodicity is supported, but not stellarator symmetry
//----------------------------------------------------------------------------------------------------  
   
	// Numerical prefactors and index variables
   double muc = 1.0e-7;
   register int ip, i, j, k, is;
   double factor = muc * 2 / Nfils;
   
	// Relevant geometric and EM quantities
	double l, ex, ey, ez, bx, by, bz, bxx, byy, bzz;
   double xx, yy, zz, ri, xi, yi, zi, rf, xf, yf, zf;
   double cur, coef, symfac, rot_cos, rot_sin;

	// Initialize field components to zero 
   bx  = 0.0; by  = 0.0; bz  = 0.0;
   bxx = 0.0; byy = 0.0; bzz = 0.0;

   for(ip = 1; ip < Nfp + 1; ip++) { 
      
		//Find coefficients of transformation matrix
      rot_cos = cosnfp(ip);
      rot_sin = sinnfp(ip);
      is = Ns; 
      
		while(is > -1) {  
         
         //Find symmetric point on other field periods
			symfac = pow(-1,is);         
         xx = x * rot_cos + y * rot_sin; 
         yy = (-x * rot_sin + y * rot_cos) * symfac; 
         zz = z * symfac;
         
			for(i = 0; i < iCoil; i++) {
            
				//Load current of i-th coil 
            cur = currents[i];

            for(j = 0; j < Nfils; j++) {
               
					//Calculate first point of current filament before main loop
               xi = xx - mfilx[i*Nfils*(Nseg+1)+j*(Nseg+1)];
               yi = yy - mfily[i*Nfils*(Nseg+1)+j*(Nseg+1)];
               zi = zz - mfilz[i*Nfils*(Nseg+1)+j*(Nseg+1)];
               ri = sqrt(xi * xi + yi * yi + zi * zi);           
                
               for(k = 0; k < Nseg; k++) {
						
						// Use Hanson Hirschman method for finding B for each segment
                  xf = xx - mfilx[i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1];
                  yf = yy - mfily[i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1];
                  zf = zz - mfilz[i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1];
                  rf = sqrt(xf * xf + yf * yf + zf * zf);
               
                  l = sqrt((xf-xi) * (xf-xi) + (yf-yi) * (yf-yi) + (zf-zi) * (zf-zi));  

                  ex = ((xf-xi)/l);  
                  ey = ((yf-yi)/l);  
                  ez = ((zf-zi)/l);  

						// Prefactor for cross product 
                  coef = cur * l * (ri+rf) / (ri*rf * ((ri+rf) * (ri+rf) - l * l));   

                  bx += coef * (ey * zi - ez * yi);
                  by += coef * (ez * xi - ex * zi);
                  bz += coef * (ex * yi - ey * xi);
              
                  //End of segment k becomes beginning of segment k+1
                  xi = xf;
                  yi = yf;
                  zi = zf;
                  ri = rf;
               }
            }
			}     
       
		   //Rotate back to first field period
         bxx += (bx * rot_cos - by * rot_sin) * symfac;  
         byy += (bx * rot_sin + by * rot_cos); 
         bzz += bz;

			// Reset for next field period calculation
         bx = 0;
         by = 0;
         bz = 0;          

			// Update symmetry factor
         is = is - 1;
      }
   }  

	// Apply scaled permeability factor to get final result    
	*Bx = bxx * factor;
   *By = byy * factor;
   *Bz = bzz * factor;
}

// Driving code for GPU calculations 
// In the future, more kernels will be called after magnetic_field
void CalculateFieldParallelGPU(int start, int end) {
	
   magnetic_field1(mfilx, mfily, mfilz, xsurf, ysurf, zsurf,
                  currents, Ncoil * Nfils, Nfils, Nfp, Nseg+1, size_fp,
                  Bmfilx, Bmfily, Bmfilz, start, end);  
}
