#include "bfield.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <omp.h>

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

double cosnfp(int ip){ 
//----------------------------------------------------------------------------------------------------
// Cosine of vector component at the ip-th field period
//----------------------------------------------------------------------------------------------------
   return cos((ip - 1) * M_PI * 2 / Nfp);
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

double sinnfp(int ip){ 
//----------------------------------------------------------------------------------------------------
// Sine of vector component at the ip-th field period
//----------------------------------------------------------------------------------------------------
   return sin((ip - 1) * M_PI * 2 / Nfp);

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

void CalculateMultiFieldSym(double x, double y, double z, 
                            double* Bx, double* By, double* Bz){
//----------------------------------------------------------------------------------------------------
// Calculates the field due to multi-filament representation at the point (x,y,z)
// Based on the method of Hanson and Hirschman.
// Periodicity is supported, but not stellarator symmetry
//
// The input point is used to find the same point on each period,
// (x,y,z) --> (xxi_1,yy_1,zz_1) ... (xx_nfp-1,yy_nfp-1,zz_nfp-1)
// The field is then calculated at each of the nfp points and rotated back.
// Rotations are handled by the cosnfp,sinnfp functions. 
//----------------------------------------------------------------------------------------------------  
   
   double muc = 1.0e-7; // (mu / 4 * pi)
   int ip, i, j, k, is;
   double factor = muc * 2.0 / Nfils;
   double l, ex, ey, ez; 
   double bx, by, bz, bxx, byy, bzz;
   double xx, yy, zz, ri, xi, yi, zi, rf, xf, yf, zf;
   double cur, eps, eta, coef, symfac, rot_cos, rot_sin;

	bx = 0.0;
	by = 0.0;
	bz = 0.0;

   bxx = 0.0;
   byy = 0.0;
   bzz = 0.0;
   
   zz = z;

   for(ip = 1; ip < Nfp + 1; ip++) { 

      //Find coefficients of transformation matrix
      rot_cos = cosnfp(ip);
      rot_sin = sinnfp(ip);
      is = Ns; 

      while(is>-1) {  

         // For stellarator symmetry case 
         symfac = pow(-1, is);         
      
         // Find symmetric point on other field periods
         xx = (  x * rot_cos + y * rot_sin ); 
         yy = ( -x * rot_sin + y * rot_cos ) * symfac; 
         zz = z * symfac;

			// Loop over all coils
         for(i = 0; i < iCoil; i++) {	
	    		
				// Loop over all filaments in coil
				for(j = 0; j < Nfils; j++) {

               //Store first point of current filament before main loop
               xi = xx - mfilx[i*Nfils*(Nseg+1)+j*(Nseg+1)];
               yi = yy - mfily[i*Nfils*(Nseg+1)+j*(Nseg+1)];
               zi = zz - mfilz[i*Nfils*(Nseg+1)+j*(Nseg+1)];
               ri = sqrt(pow(xi,2) + pow(yi,2) + pow(zi,2));           
            
					// Loop over all segments in a filament
               for(k = 0; k < Nseg; k++) {

                  // Calculate vector from input point to end of segment
                  xf = xx - mfilx[i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1];
                  yf = yy - mfily[i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1];
                  zf = zz - mfilz[i*Nfils*(Nseg+1)+j*(Nseg+1)+k+1];
                  rf = sqrt(pow(xf,2) + pow(yf,2) + pow(zf,2));
               
                  l = sqrt( pow(xf-xi,2) + pow(yf-yi,2) + pow(zf-zi,2) );  

                  ex = (xf - xi) / l;
                  ey = (yf - yi) / l;  
                  ez = (zf - zi) / l;  
       
						// Factor common to dBx, dBy and dBz
						// Note that the factor of two and current are absent here
	               coef = l * (ri+rf) / (ri*rf * (pow(ri + rf,2) - pow(l,2)) );	

						bx += coef * (ey*zi-ez*yi); 
                  by += coef * (ez*xi-ex*zi);
                  bz += coef * (ex*yi-ey*xi);

/* MODIFIED 04/14/21: get rid of unnecessary current multiplication in innermost loop (759) 
 * This was implicit in how coef was defined previously 	
                  bx += coef * (ey*zi-ez*yi);
                  by += coef * (ez*xi-ex*zi);
                  bz += coef * (ex*yi-ey*xi);
*/              
                  //End of segment k becomes beginning of segment k+1
                  xi = xf;
                  yi = yf;
                  zi = zf;
                  ri = rf;
               }
            }
         	bx *= currents[i];
				by *= currents[i];
				bz *= currents[i];
			}     
        
         //Rotate back to first field period
         bxx += (bx * rot_cos - by * rot_sin) * symfac;  
         byy += (bx * rot_sin + by * rot_cos); 
         bzz += bz;

         bx = 0;
         by = 0;
         bz = 0;          
  
         is = is - 1;
      }
   }  
   
   *Bx = bxx * factor;
   *By = byy * factor;
   *Bz = bzz * factor;
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

