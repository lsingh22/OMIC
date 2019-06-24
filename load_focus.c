#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
#include <netcdf.h>
#include <stdlib.h>
#include <math.h>
int main(int argc, char **argv) {

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////  START LOAD_FOCUS.C ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Read the input .nc file//

/*
for (int i=0;i<argc;i++){
   printf("argv[%d]: %s\n", i, argv[i]);
}   
*/

   char* file;
   file = argv[1];

   //Define and allocate pointers// 

   //Decide what to print out
   int print_sfil, print_mfil, print_vect, print_mmfil;
   print_sfil = 1; //print out single filament file if 1
   print_mfil = 0; //print out multifilament file if 1
   print_vect = 0; //print out normal and binormal vectors
   print_mmfil = 0;

   FILE* fp;
   char* line = NULL;
   size_t len = 0;
   int count = 0;
   char* data;
   int* n;
   int* m;
   int i,j,k,N;
   int sum = 0;

   int ncid, varid;
   int* Ncoils;
   int* Nfp;
   int* isSym;
   //TODO: int* NFcoil;
   //int NFcoil[11] = { 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8 };
 
   //int NFcoil[48] = { 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5};
   int NFcoil[48] = { 16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16};
   Ncoils = (int *) malloc(sizeof(int));
//   NFcoil = (int *) malloc(sizeof(int));   
   Nfp = (int *) malloc(sizeof(int));
   isSym = (int *) malloc(sizeof(int));

   //TODO: Notify CZ about NFcoil bug

   //nc_open("example/focus_hsx.m12_07.nc", NC_NOWRITE, &ncid);
   nc_open(file, NC_NOWRITE, &ncid);
   nc_inq_varid(ncid, "Ncoils", &varid);
   nc_get_var_int(ncid, varid, Ncoils);
   nc_inq_varid(ncid, "Nfp", &varid);
   nc_get_var_int(ncid, varid, Nfp);
   nc_inq_varid(ncid, "IsSymmetric", &varid);
   nc_get_var_int(ncid, varid, isSym);

   //TODO: Fix the NFcoil storage after bug is fixed
   //nc_inq_varid(ncid, "NFcoil", &varid);
   //nc_get_var_int(ncid, varid, NFcoil);

   //Determine the total number of FS amplitudes// 
   for (i =0; i < (*Ncoils);i++){
        sum = sum + NFcoil[i] +3;
   }

   //Store currents and FS harmonics for each coil//

   double* coildata;
   double* centroids;
   double* coilamps;
   double* currents;  
   

   coildata = (double *) malloc(59*101*sizeof(double));
   centroids = (double *) malloc((*Ncoils)*3*sizeof(double));
   coilamps = (double *) malloc(sum*6*sizeof(double));
   currents = (double *) malloc((*Ncoils)*sizeof(double));
   nc_inq_varid(ncid, "coilspace", &varid);
   nc_get_var_double(ncid, varid, coildata);
   nc_close(ncid);

   int ind = 0;
   int ind_arr[*Ncoils];

   for(i=0;i<(*Ncoils);i++){
      if (i==0) 
      {
        *(currents + i) = coildata[ind];
      }
      else
      {
        *(currents + i) = coildata[ind+i];
      }
      for(j=0;j<6*(NFcoil[i]+1)-3;j++){
         *(coilamps + ind + j) = coildata[ind + j + i + 1 ];
      }
      ind_arr[i] = ind;
      ind = ind + ((NFcoil[i])*6+3);
   }
   //Store centroids from coilamps array//

   for(i=0;i<(*Ncoils);i++){
        centroids[i*3 + 0] = coilamps[ind_arr[i]];
        centroids[i*3 + 1] = coilamps[ind_arr[i] + 2*(NFcoil[i] + 1)-1];
        centroids[i*3 + 2] = coilamps[ind_arr[i] + 4*(NFcoil[i] + 1)-2];
        //printf("%f  %f  %f\n", centroids[i*3 + 1], centroids[i*3 + 2], centroids[i*3 + 3]);

   }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////  START SINGLE_FIL.C /////////////////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Function that constructs single filament from modes 
//Input should be pointer to vector that has modes for each coil
//Output should be pointer to XYZ points

   
   int Nseg = 128;
   double pi = M_PI;
   double theta; //parameterizes each coil via Fourier Series
   double x,y,z,x0,y0,z0;
   double* sfilx; 
   double* sfily;
   double* sfilz;
  
   //printf("%d\n",(*Ncoils)*Nseg*sizeof(double));
   sfilx = (double *) malloc((*Ncoils)*Nseg*sizeof(double)); 
   sfily = (double *) malloc((*Ncoils)*Nseg*sizeof(double));
   sfilz = (double *) malloc((*Ncoils)*Nseg*sizeof(double));
   
   //Calculate the xyz for the single filament
   if(print_mfil == 1 | print_sfil == 1)
   {
      printf(" periods 1\n begin filament\n mirror NIL\n");
   }
   //printf("%f\n",currents[0]);
   for(i=0;i<(*Ncoils);i++){
      x0=0;y0=0;z0=0;
      for(j=0;j<Nseg;j++){
         theta = ((2*pi)/Nseg)*j;
	 x=0;y=0;z=0;
         for(k=0;k<NFcoil[i]+1;k++){ //add the cosine components
           x = x + coilamps[ ind_arr[i] + k ]*cos(k*theta);
           y = y + coilamps[ ind_arr[i] + k + 2*NFcoil[i] + 1 ]*cos(k*theta);
           z = z + coilamps[ ind_arr[i] + k + 4*NFcoil[i] + 2 ]*cos(k*theta);                
         }
         for(k=1;k<NFcoil[i]+1;k++){ //add the sine components
           x = x + coilamps[ ind_arr[i] +   NFcoil[i] + 0 + k ]*sin(k*theta);
           y = y + coilamps[ ind_arr[i] + 3*NFcoil[i] + 1 + k ]*sin(k*theta);
           z = z + coilamps[ ind_arr[i] + 5*NFcoil[i] + 2 + k ]*sin(k*theta);
         }
         *(sfilx + i*Nseg + j ) = x;
	 *(sfily + i*Nseg + j ) = y;
         *(sfilz + i*Nseg + j ) = z;
         if( j == 0 ) 
         {
	    x0 = x;
	    y0 = y;
            z0 = z;
	 }
         if( print_sfil == 1 )
	 {
            if( j == Nseg - 1 )
            {
               printf(" %.15f  %.15f  %.15f  %.8f \n",x,y,z,coildata[ind_arr[i]+i]);
               printf(" %.15f  %.15f  %.15f  %.8f Mod 1\n",x0,y0,z0,coildata[ind_arr[i]+i]);
            }else{
               printf(" %.15f  %.15f  %.15f  %.8f\n",x,y,z,coildata[ind_arr[i]+i]);
            }
         }  
      } 
   }
      if(print_sfil == 1)
   {
      printf("end\n");
   }

   


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////  START MULTI_FIL.C //////////////////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//This function will take the modes that define the centroid of the coils and will take the alpha vectors and the axis information
//This code will construct the XYZ points for the multifilaments and will return these XYZ points so that they can be plotted
//
//Use centroid and gram schmidt process to define local coordinate frame
//From local coordinate frame construct finite build from specified dimensions

   double hwid = 0.030; //normal half length 0.060
   double hlen = 0.015; //binormal half length 0.030
   double norm;

   double* tx; tx = (double *) malloc((*Ncoils)*Nseg*sizeof(double));
   double* ty; ty = (double *) malloc((*Ncoils)*Nseg*sizeof(double));
   double* tz; tz = (double *) malloc((*Ncoils)*Nseg*sizeof(double));
   double* nx; nx = (double *) malloc((*Ncoils)*Nseg*sizeof(double));
   double* ny; ny = (double *) malloc((*Ncoils)*Nseg*sizeof(double));
   double* nz; nz = (double *) malloc((*Ncoils)*Nseg*sizeof(double));
   double* bx; bx = (double *) malloc((*Ncoils)*Nseg*sizeof(double));
   double* by; by = (double *) malloc((*Ncoils)*Nseg*sizeof(double));
   double* bz; bz = (double *) malloc((*Ncoils)*Nseg*sizeof(double));
   double* nxa; nxa = (double *) malloc((*Ncoils)*Nseg*sizeof(double));
   double* nya; nya = (double *) malloc((*Ncoils)*Nseg*sizeof(double));
   double* nza; nza = (double *) malloc((*Ncoils)*Nseg*sizeof(double));
   double* bxa; bxa = (double *) malloc((*Ncoils)*Nseg*sizeof(double));
   double* bya; bya = (double *) malloc((*Ncoils)*Nseg*sizeof(double));
   double* bza; bza = (double *) malloc((*Ncoils)*Nseg*sizeof(double));

   double* alpha; alpha = (double *) malloc((*Ncoils)*Nseg*sizeof(double));  

  //Assume a rectangular cross section for now
   double* mfilx; mfilx = (double *) malloc((*Ncoils)*Nseg*5*sizeof(double));
   double* mfily; mfily = (double *) malloc((*Ncoils)*Nseg*5*sizeof(double));
   double* mfilz; mfilz = (double *) malloc((*Ncoils)*Nseg*5*sizeof(double));

 
   //Use FS representation to calculate the unit tangent vector at each point on coil 
   //TODO:Hardcode alpha for now; for laziness, set two extra sin/cos weighting parameters
   double a,b,c,d;
   a = 0;
   b = 0;
   c = 1;
   d = 1;
   for(i=0;i<(*Ncoils);i++){
      x0=0;y0=0;z0=0;
      for(j=0;j<Nseg;j++){
         theta = ((2*pi)/Nseg)*j;
         x=0;y=0;z=0;norm=0;
         *(alpha + i*Nseg + j) = a*cos(c*theta) + b*sin(d*theta);
         //*(alpha + i*Nseg + j) = theta;
         //(alpha + i*Nseg + j) = exp(theta);
         for(k=0;k<NFcoil[i]+1;k++){ //add the cosine components
           x = x - k*coilamps[ ind_arr[i] + k ]*sin(k*theta);
           y = y - k*coilamps[ ind_arr[i] + k + 2*NFcoil[i] + 1 ]*sin(k*theta);
           z = z - k*coilamps[ ind_arr[i] + k + 4*NFcoil[i] + 2 ]*sin(k*theta);
         }
         for(k=1;k<NFcoil[i]+1;k++){ //add the sine components
           x = x + k*coilamps[ ind_arr[i] +   NFcoil[i] + 0 + k ]*cos(k*theta);
           y = y + k*coilamps[ ind_arr[i] + 3*NFcoil[i] + 1 + k ]*cos(k*theta);
           z = z + k*coilamps[ ind_arr[i] + 5*NFcoil[i] + 2 + k ]*cos(k*theta);
         }
         norm = sqrt(x*x + y*y + z*z);
         *(tx + i*Nseg + j ) = x/norm;
         *(ty + i*Nseg + j ) = y/norm;
         *(tz + i*Nseg + j ) = z/norm;
         if( j == 0 )
         {
            x0 = x/norm;
            y0 = y/norm;
            z0 = z/norm;
         }
         if( j == Nseg - 1 )
         {
            //printf(" %.15f  %.15f  %.15f \n",x/norm,y/norm,z/norm);
            //printf(" %.15f  %.15f  %.15f Mod 1\n",x0,y0,z0);
         }else{
            //printf(" %.15f  %.15f  %.15f  \n",x/norm,y/norm,z/norm);
         }
      }
   }
   //Carry out the calculation of build directions, rotation of build directions, and finally xyz for the finite build and mfils 
   double dot;
   double sum1 = 0;
   double sum2 = 0;
   double sum3 = 0;
   double* sfilxa;
   double* sfilya;
   double* sfilza;

   //printf("%d\n",(*Ncoils)*Nseg*sizeof(double));
   sfilxa = (double *) malloc((*Ncoils)*Nseg*sizeof(double));
   sfilya = (double *) malloc((*Ncoils)*Nseg*sizeof(double));
   sfilza = (double *) malloc((*Ncoils)*Nseg*sizeof(double));

   //Store the vector pointing from coil centroid to point on coil for each coil
   for(i=0;i<(*Ncoils);i++){
      for(j=0;j<Nseg;j++){
         *(sfilxa + i*Nseg + j ) = *(sfilx + i*Nseg + j ) - *(centroids + i*3 + 0);
         *(sfilya + i*Nseg + j ) = *(sfily + i*Nseg + j ) - *(centroids + i*3 + 1);
         *(sfilza + i*Nseg + j ) = *(sfilz + i*Nseg + j ) - *(centroids + i*3 + 2);
	 //printf("%f %f %f\n", *(sfilxa + i*Nseg + j ), *(sfilya + i*Nseg + j ), *(sfilza + i*Nseg + j ));
      }
   }

   for(i=0;i<(*Ncoils)*Nseg;i++){
      x=0;y=0;z=0;
      dot = *(sfilxa + i) * *(tx + i) + *(sfilya + i) * *(ty + i) + *(sfilza + i) * *(tz + i);
      x = *(sfilxa + i) - dot* *(tx + i);
      y = *(sfilya + i) - dot* *(ty + i);
      z = *(sfilza + i) - dot* *(tz + i);
      norm = sqrt(x*x + y*y + z*z);
      *(nx+i) = x/norm;
      *(ny+i) = y/norm;
      *(nz+i) = z/norm;
      *(bx+i) = *(ty+i) * *(nz+i) - *(tz+i) * *(ny+i);
      *(by+i) = *(tz+i) * *(nx+i) - *(tx+i) * *(nz+i);
      *(bz+i) = *(tx+i) * *(ny+i) - *(ty+i) * *(nx+i);
      //printf("Normal vector: %f %f %f\n",*(nx+ i),*(ny+i),*(nz+i));
      //printf("Binormal vector: %f %f %f\n",*(bx+ i),*(by+i),*(bz+i));
      sum1 = sum1 + *(nx +i )* *(tx + i) + *(ny +i )* *(ty + i) + *(nz +i )* *(tz + i);   
      sum2 = sum2 + *(bx +i )* *(tx + i) + *(by +i )* *(ty + i) + *(bz +i )* *(tz + i);  
      sum3 = sum3 + *(bx +i )* *(nx + i) + *(by +i )* *(ny + i) + *(bz +i )* *(nz + i);
      //printf("Sums are: %f %f %f\n",sum1,sum2,sum3);

   //Rotate the normal and binormal vectors by an angle alpha about the tangent vector
      *(nxa+i) = *(nx+i)*cos(*(alpha+i)) + *(bx+i)*sin(*(alpha+i));
      *(nya+i) = *(ny+i)*cos(*(alpha+i)) + *(by+i)*sin(*(alpha+i));
      *(nza+i) = *(nz+i)*cos(*(alpha+i)) + *(bz+i)*sin(*(alpha+i));

      *(bxa+i) = -*(nx+i)*sin(*(alpha+i)) + *(bx+i)*cos(*(alpha+i));
      *(bya+i) = -*(ny+i)*sin(*(alpha+i)) + *(by+i)*cos(*(alpha+i));
      *(bza+i) = -*(nz+i)*sin(*(alpha+i)) + *(bz+i)*cos(*(alpha+i));
      //sum1 = sum1 + *(nx +i )* *(tx + i) + *(ny +i )* *(ty + i) + *(nz +i )* *(tz + i);
      //sum2 = sum2 + *(bx +i )* *(tx + i) + *(by +i )* *(ty + i) + *(bz +i )* *(tz + i);
      //sum3 = sum3 + *(bx +i )* *(nx + i) + *(by +i )* *(ny + i) + *(bz +i )* *(nz + i);  

     //Using specified half length and width, determine rectangular winding pack xyz points 
      for(j=0;j<5;j++){
         *(mfilx + 5*i) =     *(sfilx+i) + hwid * *(nxa+i) + hlen * *(bxa+i);
         *(mfilx + 5*i + 1) = *(sfilx+i) - hwid * *(nxa+i) + hlen * *(bxa+i);
         *(mfilx + 5*i + 2) = *(sfilx+i) - hwid * *(nxa+i) - hlen * *(bxa+i);
         *(mfilx + 5*i + 3) = *(sfilx+i) + hwid * *(nxa+i) - hlen * *(bxa+i);
         *(mfilx + 5*i + 4) = *(sfilx+i) + hwid * *(nxa+i) + hlen * *(bxa+i);

         *(mfily + 5*i) =     *(sfily+i) + hwid * *(nya+i) + hlen * *(bya+i);
         *(mfily + 5*i + 1) = *(sfily+i) - hwid * *(nya+i) + hlen * *(bya+i);
         *(mfily + 5*i + 2) = *(sfily+i) - hwid * *(nya+i) - hlen * *(bya+i);
         *(mfily + 5*i + 3) = *(sfily+i) + hwid * *(nya+i) - hlen * *(bya+i);
         *(mfily + 5*i + 4) = *(sfily+i) + hwid * *(nya+i) + hlen * *(bya+i);

         *(mfilz + 5*i)     = *(sfilz+i) + hwid * *(nza+i) + hlen * *(bza+i);
         *(mfilz + 5*i + 1) = *(sfilz+i) - hwid * *(nza+i) + hlen * *(bza+i);
         *(mfilz + 5*i + 2) = *(sfilz+i) - hwid * *(nza+i) - hlen * *(bza+i);
         *(mfilz + 5*i + 3) = *(sfilz+i) + hwid * *(nza+i) - hlen * *(bza+i);
         *(mfilz + 5*i + 4) = *(sfilz+i) + hwid * *(nza+i) + hlen * *(bza+i);
	 
	 //printf("%f,%f,%f,%f,%f\n",*(mfilz + 5*i + 1),*(mfilz + 5*i + 1),*(mfilz + 5*i + 1),*(mfilz + 5*i + 1),*(mfilz + 5*i + 1));
      }    
   }
   if(print_mfil == 1){
      for(i=0;i<(*Ncoils);i++){
	 for(j=0;j<Nseg;j++){
            if( j == Nseg - 1) 
	    {
               for(k=0;k<5;k++){
                   printf(" %.15f  %.15f  %.15f  %.8f \n",*(mfilx+i*Nseg*5+j*5+k),*(mfily+i*Nseg*5+j*5+k),*(mfilz+i*Nseg*5+j*5+k),coildata[0]);
               }
               for(k=0;k<4;k++){
                   printf(" %.15f  %.15f  %.15f  %.8f \n",*(mfilx+i*Nseg*5+k),*(mfily+i*Nseg*5+k),*(mfilz+i*Nseg*5+k),coildata[0]);
	       }
               printf(" %.15f  %.15f  %.15f  %.8f Mod 1\n",*(mfilx+i*Nseg*5+4),*(mfily+i*Nseg*5+4),*(mfilz+i*Nseg*5+4),coildata[0]);
	    }else{
                 for(k=0;k<5;k++){
                    printf(" %.15f  %.15f  %.15f  %.8f \n",*(mfilx+i*Nseg*5+j*5+k),*(mfily+i*Nseg*5+j*5+k),*(mfilz+i*Nseg*5+j*5+k),coildata[0]);
		 }
	     }
          }
       }
       printf("end\n");
   }
 
//////////////////////////////////// START MULTIFILAMENT CONST ////////////
   
   int s,t;
   int l=0;
   s = 3;
   t = 2;
   double mtheta;
   double slen = hlen / (s/2);
   double tlen = hwid /(2*t);
   double sloc, tloc;
   double* mmfilx; mmfilx = (double *) malloc( (*Ncoils)*Nseg*s*t*sizeof(double));
   double* mmfily; mmfily = (double *) malloc( (*Ncoils)*Nseg*s*t*sizeof(double));
   double* mmfilz; mmfilz = (double *) malloc( (*Ncoils)*Nseg*s*t*sizeof(double));
   
   /*
   for(i=0;i<(*Ncoils);i++){
      for(j=0;j<t;j++){
         tor_loc = mflen*(-(t-1)/2 + j);
	 for(k=0;k<Nseg;k++){
            mtheta = (2*s*pi / Nseg) * k;
	    rad_loc = mfwid*(s/2 - mtheta/(2*pi));
	    *(mmfilx + i*t*Nseg + j*Nseg + k) = *(sfilx +i*Nseg + k) + rad_loc* *(nxa + i*Nseg + k) + tor_loc* *(bxa + i*Nseg + k);
	    *(mmfily + i*t*Nseg + j*Nseg + k) = *(sfily +i*Nseg + k) + rad_loc* *(nya + i*Nseg + k) + tor_loc* *(bya + i*Nseg + k);
	    *(mmfilz + i*t*Nseg + j*Nseg + k) = *(sfilz +i*Nseg + k) + rad_loc* *(nza + i*Nseg + k) + tor_loc* *(bza + i*Nseg + k);
	 }
      }
   }
   */
   //printf("%.10f %.10f",hwid,tlen);
   for(i=0;i<(*Ncoils);i++){
      for(j=0;j<Nseg;j++){
         for(k=0;k<s;k++){
	    sloc = slen*(-(s-1) + 2*k);
            for(l=0;l<t;l++){
            tloc = tlen*(-(t-1) + 2*l);
            *(mmfilx + i*t*s*Nseg + j*t*s + k*t + l) = *(sfilx +i*Nseg + j) + sloc* *(nxa + i*Nseg + j) + tloc* *(bxa + i*Nseg + j);
            *(mmfily + i*t*s*Nseg + j*t*s + k*t + l) = *(sfily +i*Nseg + j) + sloc* *(nya + i*Nseg + j) + tloc* *(bya + i*Nseg + j);
	    *(mmfilz + i*t*s*Nseg + j*t*s + k*t + l) = *(sfilz +i*Nseg + j) + sloc* *(nza + i*Nseg + j) + tloc* *(bza + i*Nseg + j);
	    //printf("%.15f %.15f %.15f %.8f\n", *(mmfilx+ + i*t*s*Nseg + j*t*s + k*t + l), *(mmfily  + i*t*s*Nseg + j*t*s + k*t + l), *(mmfilz  + i*t*s*Nseg + j*t*s + k*t + l),coildata[0]);
	    }
         }
      }
   }

   if(print_mmfil == 1)
   {
      printf(" periods 1\n begin filament\n mirror NIL\n");

      for(i=0;i<*(Ncoils);i++){ //handle all but the last filament
         for(j=0;j<s*t-1;j++){
	    for(k=0;k<Nseg;k++){
	       printf(" %.15f %.15f %.15f %.8f\n", *(mmfilx+i*s*t*Nseg+k*t*s+j), *(mmfily+i*s*t*Nseg+k*t*s+j), *(mmfilz+i*s*t*Nseg+k*t*s+j),coildata[0]);
	    }
            printf(" %.15f %.15f %.15f %.8f\n", *(mmfilx+i*s*t*Nseg+j), *(mmfily+i*s*t*Nseg+j), *(mmfilz+i*s*t*Nseg+j),coildata[0]);
         }
         
            for(k=0;k<Nseg;k++){ //handle the last filament and then increment coil 
               printf(" %.15f %.15f %.15f %.8f\n", *(mmfilx+i*s*t*Nseg+k*t*s+s*t-1), *(mmfily+i*s*t*Nseg+k*t*s+s*t-1), *(mmfilz+i*s*t*Nseg+k*t*s+s*t-1),coildata[0]);
         }
            printf(" %.15f %.15f %.15f %.8f Mod 1\n", *(mmfilx+i*s*t*Nseg+s*t-1), *(mmfily+i*s*t*Nseg+s*t-1), *(mmfilz+i*s*t*Nseg+s*t-1),coildata[0]);
      }
      printf("end\n");
   }
    




   //For now, this writes two new files stored locally in multi runs directory
   if(print_vect == 1)
   {
      FILE* fn;
      fn = fopen("/home/luquants/multi/runs/norm.txt","w");
      fprintf(fn," periods 1\n begin filament\n mirror NIL\n"); //TODO: For use with T.Kruger's PlotCoils, to be updated
      for(i=0;i<(*Ncoils);i++){
         for(j=0;j<Nseg;j++){
            if( j == Nseg - 1 )
            {
               fprintf(fn," %.15f  %.15f  %.15f \n", *(nx + i*Nseg+j), *(ny + i*Nseg+j), *(nz + i*Nseg+j));
               fprintf(fn," %.15f  %.15f  %.15f Mod 1\n", *(nx + i*Nseg), *(ny + i*Nseg), *(nz + i*Nseg));
            }else{
               fprintf(fn," %.15f  %.15f  %.15f  \n", *(nx + i*Nseg+j), *(ny + i*Nseg+j), *(nz + i*Nseg+j)); 
            }
         }
      }
      fprintf(fn,"end\n"); //TODO: See above TODO
      fclose(fn);
   
  

      FILE* fb;
      fb = fopen("/home/luquants/multi/runs/binorm.txt","w");
      fprintf(fb," periods 1\n begin filament\n mirror NIL\n"); //TODO: For use with T.Kruger's PlotCoils, to be updated
      for(i=0;i<(*Ncoils);i++){
         for(j=0;j<Nseg;j++){
            if( j == Nseg - 1 )
            {
               fprintf(fb," %.15f  %.15f  %.15f \n", *(bx + i*Nseg+j), *(by + i*Nseg+j), *(bz + i*Nseg+j));
               fprintf(fb," %.15f  %.15f  %.15f Mod 1\n", *(bx + i*Nseg), *(by + i*Nseg), *(bz + i*Nseg));
            }else{
               fprintf(fb," %.15f  %.15f  %.15f  \n", *(bx + i*Nseg+j), *(by + i*Nseg+j), *(bz + i*Nseg+j)); 
            }
         }
      }
      fprintf(fb,"end"); //TODO: See above TODO
      fclose(fb);
   }
}
