
#include "multi_fil.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

//GLOBALS SCOPED IN SOURCE FILE

double hwid;
double hlen;
int Nseg;
int Nradfil;
int Ntorfil;

double* tx;
double* ty;
double* tz;
double* nx;
double* ny;
double* nz;
double* bx; 
double* by; 
double* bz; 

double* nxa;
double* nya;
double* nza;
double* bxa;
double* bya;
double* bza;
double* alp;   

double* cx;
double* cy;
double* cz;
double* sfilx;
double* sfily;
double* sfilz;

double* mfilx;
double* mfily;
double* mfilz;

double* finx;
double* finy;
double* finz;


void CalculateBuildDirections(void){

   double norm;
   double theta, x,y,z;
   double pi = M_PI;
   tx = (double *) malloc((Ncoils)*Nseg*sizeof(double));
   ty = (double *) malloc((Ncoils)*Nseg*sizeof(double));
   tz = (double *) malloc((Ncoils)*Nseg*sizeof(double));
   nx = (double *) malloc((Ncoils)*Nseg*sizeof(double));
   ny = (double *) malloc((Ncoils)*Nseg*sizeof(double));
   nz = (double *) malloc((Ncoils)*Nseg*sizeof(double));
   bx = (double *) malloc((Ncoils)*Nseg*sizeof(double));
   by = (double *) malloc((Ncoils)*Nseg*sizeof(double));
   bz = (double *) malloc((Ncoils)*Nseg*sizeof(double));
    
   nxa = (double *) malloc((Ncoils)*Nseg*sizeof(double));
   nya = (double *) malloc((Ncoils)*Nseg*sizeof(double));
   nza = (double *) malloc((Ncoils)*Nseg*sizeof(double));
   bxa = (double *) malloc((Ncoils)*Nseg*sizeof(double));
   bya = (double *) malloc((Ncoils)*Nseg*sizeof(double));
   bza = (double *) malloc((Ncoils)*Nseg*sizeof(double));
   
   int i,j,k; 
   double dot;
   //Calculate unit tangent vector
   for(i=0;i<Ncoils;i++){
      for(j=0;j<Nseg;j++){
         theta = ((2*pi)/Nseg)*j;
         x=0;y=0;z=0;norm=0;
         for(k=0;k<NFcoil+1;k++){ //add the cosine components
           x = x - k*coilamps[ ind_arr[i] + k ]*sin(k*theta);
           y = y - k*coilamps[ ind_arr[i] + k + 2*NFcoil + 1 ]*sin(k*theta);
           z = z - k*coilamps[ ind_arr[i] + k + 4*NFcoil + 2 ]*sin(k*theta);
         }
         for(k=1;k<NFcoil+1;k++){ //add the sine components
           x = x + k*coilamps[ ind_arr[i] +   NFcoil + 0 + k ]*cos(k*theta);
           y = y + k*coilamps[ ind_arr[i] + 3*NFcoil + 1 + k ]*cos(k*theta);
           z = z + k*coilamps[ ind_arr[i] + 5*NFcoil + 2 + k ]*cos(k*theta);
         }
         norm = sqrt(x*x + y*y + z*z);
         *(tx + i*Nseg + j ) = x/norm;
         *(ty + i*Nseg + j ) = y/norm;
         *(tz + i*Nseg + j ) = z/norm;
 
      }
   }

   //printf("%d\n",(*Ncoils)*Nseg*sizeof(double));
   double* sfilxa; sfilxa = (double *) malloc((Ncoils)*Nseg*sizeof(double));
   double* sfilya; sfilya = (double *) malloc((Ncoils)*Nseg*sizeof(double));
   double* sfilza; sfilza = (double *) malloc((Ncoils)*Nseg*sizeof(double));

   //Calculate vector pointing from coil centroid to point on coil for each coil
   for(i=0;i<Ncoils;i++){
      for(j=0;j<Nseg;j++){
         *(sfilxa + i*Nseg + j ) = *(sfilx + i*Nseg + j ) - *(cx + i);
         *(sfilya + i*Nseg + j ) = *(sfily + i*Nseg + j ) - *(cy +i);
         *(sfilza + i*Nseg + j ) = *(sfilz + i*Nseg + j ) - *(cz + i);
	 //printf("%f %f %f\n", *(sfilxa + i*Nseg + j ), *(sfilya + i*Nseg + j ), *(sfilza + i*Nseg + j ));
      }
   }

   alp = (double *) malloc(Ncoils*Nseg*sizeof(double)); 

   for(i=0;i<Ncoils*Nseg;i++){
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
   //Rotate the normal and binormal vectors by an angle alpha about the tangent vector
      *(alp+i)=0.0;     
      *(nxa+i) = *(nx+i)*cos(*(alp+i)) + *(bx+i)*sin(*(alp+i));
      *(nya+i) = *(ny+i)*cos(*(alp+i)) + *(by+i)*sin(*(alp+i));
      *(nza+i) = *(nz+i)*cos(*(alp+i)) + *(bz+i)*sin(*(alp+i));

      *(bxa+i) = -*(nx+i)*sin(*(alp+i)) + *(bx+i)*cos(*(alp+i));
      *(bya+i) = -*(ny+i)*sin(*(alp+i)) + *(by+i)*cos(*(alp+i));
      *(bza+i) = -*(nz+i)*sin(*(alp+i)) + *(bz+i)*cos(*(alp+i));
//printf("%f %f %f %f %f %f\n", *(nx+i),*(ny+i),*(nz+i),*(bx+i),*(by+i),*(alp+i)); 
   }
}


void CalculateMultiFilaments(void){

   int i,j,k,l;
   //Set a length and width scale for placing the filements
   //len and wid are the true length and width of the finite build
   double gridlen = len / (2*Nradfil);
   double gridwid = wid / (2*Ntorfil);
   int Nfils = Nradfil*Ntorfil;

   mfilx = (double*) malloc(Ncoils*Nfils*Nseg*sizeof(double));
   mfily = (double*) malloc(Ncoils*Nfils*Nseg*sizeof(double));
   mfilz = (double*) malloc(Ncoils*Nfils*Nseg*sizeof(double));
 
   for(i=0;i<Ncoils;i++){
      for(j=0;j<Ntorfil;j++){
         for(k=0;k<Nradfil;k++){
            for(l=0;l<Nseg;l++){
            *(mfilx + i*Nfils*Nseg + j*Nseg*Nradfil + k*Nseg + l) = *(sfilx +i*Nseg + l) + (gridlen* (-(Nradfil-1)+2*k))* *(nxa + i*Nseg + l) + (gridwid* (-(Ntorfil-1)+2*j))* *(bxa + i*Nseg + l);
            *(mfily + i*Nfils*Nseg + j*Nseg*Nradfil + k*Nseg + l) = *(sfily +i*Nseg + l) + (gridlen* (-(Nradfil-1)+2*k))* *(nya + i*Nseg + l) + (gridwid* (-(Ntorfil-1)+2*j))* *(bya + i*Nseg + l);
	    *(mfilz + i*Nfils*Nseg + j*Nseg*Nradfil + k*Nseg + l) = *(sfilz +i*Nseg + l) + (gridlen* (-(Nradfil-1)+2*k))* *(nza + i*Nseg + l) + (gridwid* (-(Ntorfil-1)+2*j))* *(bza + i*Nseg + l);
	    }
         }
      }
   }
}

//void CalculateFiniteBuild(void){}

void WriteMultiFilaments(void){

   int i,j,k;
   FILE* fb;
   fb = fopen("./outputfiles/mfil.out","w");
   fprintf(fb, "periods 1\n begin filament\n mirror NIL\n");
   int Nfils = Ntorfil*Nradfil;
   
 /*  for(i=0;i<Ncoils;i++){
      for(j=0;j<Nfils;j++){
         for(k=0;k<Nseg;k++){
         printf("%.15f %.15f %.15f %.8f \n", *(mfilx+i*Nseg*Nfils+j*Nseg+k), *(mfily+i*Nseg*Nfils+j*Nseg+k), *(mfilz+i*Nseg*Nfils+j*Nseg+k), *(currents+i));
         }
      printf("%.15f %.15f %.15f %.8f Mod %d %d\n", *(mfilx+i*Nseg*Nfils+j*Nseg), *(mfily+i*Nseg*Nfils+j*Nseg), *(mfilz+i*Nseg*Nfils+j*Nseg), *(currents+i), i+1,j+
1);
      }
   }  
*/
   for(i=0;i<Ncoils;i++){
      for(j=0;j<Nfils;j++){
         for(k=0;k<Nseg;k++){
         fprintf(fb,"%.15f %.15f %.15f %.8f \n", *(mfilx+i*Nseg*Nfils+j*Nseg+k), *(mfily+i*Nseg*Nfils+j*Nseg+k), *(mfilz+i*Nseg*Nfils+j*Nseg+k), *(currents+i));     
         }
      fprintf(fb,"%.15f %.15f %.15f %.8f Mod %d %d\n", *(mfilx+i*Nseg*Nfils+j*Nseg), *(mfily+i*Nseg*Nfils+j*Nseg), *(mfilz+i*Nseg*Nfils+j*Nseg), *(currents+i), i+1,j+1);         
      }
   }
   fprintf(fb,"end");
}

//void WriteFiniteBuild(void){}


