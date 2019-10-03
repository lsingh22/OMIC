
#include "multi_fil.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "bfield.h"
#include <omp.h>
#include "alpha.h"

int Nthreads;

//GLOBALS SCOPED IN SOURCE FILE

double len_rad;
double len_tor;

int Nseg;
int Nradfil;
int Ntorfil;
int Nfp;

int Nzeta;
int Nteta;

double* tx;
double* ty;
double* tz;
double* nx;
double* ny;
double* nz;
double* bx; 
double* by; 
double* bz; 

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

double* Bmfilx;
double* Bmfily;
double* Bmfilz;
double* Bmfil;
double* Bmfiln;


void CalculateBuildDirections(void){

   double norm;
   double theta, x,y,z;
   double pi = M_PI;
   tx = (double *) malloc((Ncoils)*(Nseg+1)*sizeof(double));
   ty = (double *) malloc((Ncoils)*(Nseg+1)*sizeof(double));
   tz = (double *) malloc((Ncoils)*(Nseg+1)*sizeof(double));
   nx = (double *) malloc((Ncoils)*(Nseg+1)*sizeof(double));
   ny = (double *) malloc((Ncoils)*(Nseg+1)*sizeof(double));
   nz = (double *) malloc((Ncoils)*(Nseg+1)*sizeof(double));
   bx = (double *) malloc((Ncoils)*(Nseg+1)*sizeof(double));
   by = (double *) malloc((Ncoils)*(Nseg+1)*sizeof(double));
   bz = (double *) malloc((Ncoils)*(Nseg+1)*sizeof(double));
      
   int i,j,k; 
   double dot;
   //Calculate unit tangent vector
   for(i=0;i<Ncoils;i++){
      for(j=0;j<Nseg+1;j++){
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
         *(tx + i*(Nseg+1) + j ) = x/norm;
         *(ty + i*(Nseg+1) + j ) = y/norm;
         *(tz + i*(Nseg+1) + j ) = z/norm;
 
      }
   }

   //printf("%d\n",(*Ncoils)*Nseg*sizeof(double));
   double* sfilxa; sfilxa = (double *) malloc((Ncoils)*(Nseg+1)*sizeof(double));
   double* sfilya; sfilya = (double *) malloc((Ncoils)*(Nseg+1)*sizeof(double));
   double* sfilza; sfilza = (double *) malloc((Ncoils)*(Nseg+1)*sizeof(double));

   //Calculate vector pointing from coil centroid to point on coil for each coil
   for(i=0;i<Ncoils;i++){
      for(j=0;j<Nseg+1;j++){
         *(sfilxa + i*(Nseg+1) + j ) = *(sfilx + i*(Nseg+1) + j ) - *(cx + i);
         *(sfilya + i*(Nseg+1) + j ) = *(sfily + i*(Nseg+1) + j ) - *(cy +i);
         *(sfilza + i*(Nseg+1) + j ) = *(sfilz + i*(Nseg+1) + j ) - *(cz + i);
	 //printf("%f %f %f\n", *(sfilxa + i*Nseg + j ), *(sfilya + i*Nseg + j ), *(sfilza + i*Nseg + j ));
      }
   }

   for(i=0;i<Ncoils*(Nseg+1);i++){
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
   }
}


void CalculateMultiFilaments(void){
   
   Unpack_alpha();
   CalculateBuildDirections();
   
   int i,j,k,l;
   //Set a length and width scale for placing the filements
   //len and wid are the true length and width of the finite build
   double gridlen = len_tor / (1*Nradfil);
   double gridwid = len_rad / (4*Ntorfil);
   double hlen_rad = len_rad / 2;
   double hlen_tor = len_tor / 2;


   int Nfils = Nradfil*Ntorfil;
   double* nxa;
   double* nya;
   double* nza;
   double* bxa;
   double* bya;
   double* bza;

   nxa = (double *) malloc((Ncoils)*(Nseg+1)*sizeof(double));
   nya = (double *) malloc((Ncoils)*(Nseg+1)*sizeof(double));
   nza = (double *) malloc((Ncoils)*(Nseg+1)*sizeof(double));
   bxa = (double *) malloc((Ncoils)*(Nseg+1)*sizeof(double));
   bya = (double *) malloc((Ncoils)*(Nseg+1)*sizeof(double));
   bza = (double *) malloc((Ncoils)*(Nseg+1)*sizeof(double));
   
   //alp = (double*) malloc(Ncoils*Nfils*(Nseg+1)*sizeof(double));
 
   mfilx = (double*) malloc(Ncoils*Nfils*(Nseg+1)*sizeof(double));
   mfily = (double*) malloc(Ncoils*Nfils*(Nseg+1)*sizeof(double));
   mfilz = (double*) malloc(Ncoils*Nfils*(Nseg+1)*sizeof(double));
 
   ffilx = (double*) malloc(Ncoils*Nfils*(Nseg+1)*5*sizeof(double));
   ffily = (double*) malloc(Ncoils*Nfils*(Nseg+1)*5*sizeof(double));
   ffilz = (double*) malloc(Ncoils*Nfils*(Nseg+1)*5*sizeof(double));
   
   // Rotate the local basis using alpha parameter

   for(i=0;i<Ncoils*(Nseg+1);i++){
      //*(alp+i) = 0.0;
      *(nxa+i) = *(nx+i)*cos(*(alp+i)) + *(bx+i)*sin(*(alp+i));
      *(nya+i) = *(ny+i)*cos(*(alp+i)) + *(by+i)*sin(*(alp+i));
      *(nza+i) = *(nz+i)*cos(*(alp+i)) + *(bz+i)*sin(*(alp+i));

      *(bxa+i) = -*(nx+i)*sin(*(alp+i)) + *(bx+i)*cos(*(alp+i));
      *(bya+i) = -*(ny+i)*sin(*(alp+i)) + *(by+i)*cos(*(alp+i));
      *(bza+i) = -*(nz+i)*sin(*(alp+i)) + *(bz+i)*cos(*(alp+i));
   
      for(j=0;j<5;j++){
         *(ffilx + 5*i) =     *(sfilx+i) + hlen_rad * *(nxa+i) + hlen_tor * *(bxa+i);
         *(ffilx + 5*i + 1) = *(sfilx+i) - hlen_rad * *(nxa+i) + hlen_tor * *(bxa+i);
         *(ffilx + 5*i + 2) = *(sfilx+i) - hlen_rad * *(nxa+i) - hlen_tor * *(bxa+i);
         *(ffilx + 5*i + 3) = *(sfilx+i) + hlen_rad * *(nxa+i) - hlen_tor * *(bxa+i);
         *(ffilx + 5*i + 4) = *(sfilx+i) + hlen_rad * *(nxa+i) + hlen_tor * *(bxa+i);

         *(ffily + 5*i) =     *(sfily+i) + hlen_rad * *(nya+i) + hlen_tor * *(bya+i);
         *(ffily + 5*i + 1) = *(sfily+i) - hlen_rad * *(nya+i) + hlen_tor * *(bya+i);
         *(ffily + 5*i + 2) = *(sfily+i) - hlen_rad * *(nya+i) - hlen_tor * *(bya+i);
         *(ffily + 5*i + 3) = *(sfily+i) + hlen_rad * *(nya+i) - hlen_tor * *(bya+i);
         *(ffily + 5*i + 4) = *(sfily+i) + hlen_rad * *(nya+i) + hlen_tor * *(bya+i);

         *(ffilz + 5*i)     = *(sfilz+i) + hlen_rad * *(nza+i) + hlen_tor * *(bza+i);
         *(ffilz + 5*i + 1) = *(sfilz+i) - hlen_rad * *(nza+i) + hlen_tor * *(bza+i);
         *(ffilz + 5*i + 2) = *(sfilz+i) - hlen_rad * *(nza+i) - hlen_tor * *(bza+i);
         *(ffilz + 5*i + 3) = *(sfilz+i) + hlen_rad * *(nza+i) - hlen_tor * *(bza+i);
         *(ffilz + 5*i + 4) = *(sfilz+i) + hlen_rad * *(nza+i) + hlen_tor * *(bza+i);	 
      }    
   }
   int iCoils = Ncoils / Nfp;

   for(i=0;i<iCoils;i++){
      for(j=0;j<Ntorfil;j++){
         for(k=0;k<Nradfil;k++){
            for(l=0;l<Nseg+1;l++){
            *(mfilx + i*Nfils*(Nseg+1) + j*(Nseg+1)*Nradfil + k*(Nseg+1) + l) = *(sfilx +i*(Nseg+1) + l) + (gridlen* (-(Nradfil-1)+2*k))* *(nxa + i*(Nseg+1) + l) + (gridwid* (-(Ntorfil-1)+2*j))* *(bxa + i*(Nseg+1) + l);
            *(mfily + i*Nfils*(Nseg+1) + j*(Nseg+1)*Nradfil + k*(Nseg+1) + l) = *(sfily +i*(Nseg+1) + l) + (gridlen* (-(Nradfil-1)+2*k))* *(nya + i*(Nseg+1) + l) + (gridwid* (-(Ntorfil-1)+2*j))* *(bya + i*(Nseg+1) + l);
	    *(mfilz + i*Nfils*(Nseg+1) + j*(Nseg+1)*Nradfil + k*(Nseg+1) + l) = *(sfilz +i*(Nseg+1) + l) + (gridlen* (-(Nradfil-1)+2*k))* *(nza + i*(Nseg+1) + l) + (gridwid* (-(Ntorfil-1)+2*j))* *(bza + i*(Nseg+1) + l);
	    }
         }
      }
   }
   int ip;
   
   for(ip=2;ip<Nfp+1;ip++){
      for(j=0;j<iCoils*Nfils*(Nseg+1);j++){
         *(mfilx + (ip-1)*(iCoils*Nfils*(Nseg+1))+j) = *(mfilx+j)*cosnfp(ip) - *(mfily+j)*sinnfp(ip);
         *(mfily + (ip-1)*(iCoils*Nfils*(Nseg+1))+j) = *(mfilx+j)*sinnfp(ip) + *(mfily+j)*cosnfp(ip);
         *(mfilz + (ip-1)*(iCoils*Nfils*(Nseg+1))+j) = *(mfilz+j);
      }
   }

}


/*
void MultiFilField(void){
  
   Bmfilx = (double*) malloc(Nzeta*Nteta*sizeof(double));
   Bmfily = (double*) malloc(Nzeta*Nteta*sizeof(double));
   Bmfilz = (double*) malloc(Nzeta*Nteta*sizeof(double));
   Bmfiln = (double*) malloc(Nzeta*Nteta*sizeof(double));
    Bmfil = (double*) malloc(Nzeta*Nteta*sizeof(double));

   int i;
   double timefield;
   double startfield, endfield;
   
   omp_set_num_threads(4);
   startfield = omp_get_wtime();
 
   #pragma omp parallel for
   for(i=0;i<Nzeta*Nteta;i++){
      
      CalculateMultiField( *(xsurf+i), *(ysurf+i), *(zsurf+i), \
                            Bmfilx+i, Bmfily+i, Bmfilz+i );    
      *(Bmfiln+i) = *(Bmfilx+i) * *(nsurfx+i) + *(Bmfily+i) * *(nsurfy+i) + \
                    *(Bmfilz+i) * *(nsurfz+i);  
      *(Bmfil+i) = sqrt( pow(*(Bmfilx+i),2) + pow(*(Bmfily+i),2) + pow(*(Bmfilz+i),2) ); 
   }
   endfield = omp_get_wtime();
   //printf("\nTotal time for multi fil field calculation: %f\n\n", endfield-startfield);   

}
*/

void MultiFilFieldSym(void){
   
   Bmfilx = (double*) malloc(Nzeta*Nteta*sizeof(double));
   Bmfily = (double*) malloc(Nzeta*Nteta*sizeof(double));
   Bmfilz = (double*) malloc(Nzeta*Nteta*sizeof(double));
   Bmfiln = (double*) malloc(Nzeta*Nteta*sizeof(double));
    Bmfil = (double*) malloc(Nzeta*Nteta*sizeof(double));

   int i,ip;
   double timefield;
   double startfield, endfield;
   int size_fp = Nzeta*Nteta / Nfp;
   
   omp_set_num_threads(Nthreads);
   startfield = omp_get_wtime();
 
   #pragma omp parallel for
   for(i=0;i<size_fp;i++){
      
      CalculateMultiFieldSym( *(xsurf+i), *(ysurf+i), *(zsurf+i), \
                            Bmfilx+i, Bmfily+i, Bmfilz+i );    
      *(Bmfiln+i) = *(Bmfilx+i) * *(nsurfx+i) + *(Bmfily+i) * *(nsurfy+i) + \
                    *(Bmfilz+i) * *(nsurfz+i);  // Need to add in area element
      *(Bmfil+i) = sqrt( pow(*(Bmfilx+i),2) + pow(*(Bmfily+i),2) + pow(*(Bmfilz+i),2) ); 
   }
   
   for(ip=1;ip<Nfp;ip++){
      for(i=0;i<size_fp;i++){
         *(Bmfil + ip*size_fp+i) = *(Bmfil+i);       
         *(Bmfiln + ip*size_fp+i) = *(Bmfiln+i);       
      }
   }


   endfield = omp_get_wtime();
   
   //Map the field to all periods
   
   printf("\nTotal time for multi fil field calculation: %f\n\n", endfield-startfield);   

}


#define MFILB_FILE_NAME "./outputfiles/mfilB.nc"
   
void WriteMultiB(void){
   //Write to NC
   int ncid, xvarid, yvarid, zvarid, bvarid, xdimid, ydimid;
   int bxvarid, byvarid, bzvarid, bnvarid;
   int dimids[2];
   
   nc_create(MFILB_FILE_NAME, NC_CLOBBER, &ncid); 
   nc_def_dim(ncid, "Nzeta", Nzeta, &xdimid);
   nc_def_dim(ncid, "Nteta", Nteta, &ydimid);
   dimids[0] = xdimid;
   dimids[1] = ydimid;
   
   nc_def_var(ncid, "xsurf", NC_DOUBLE, 2, dimids, &xvarid);
   nc_def_var(ncid, "ysurf", NC_DOUBLE, 2, dimids, &yvarid);
   nc_def_var(ncid, "zsurf", NC_DOUBLE, 2, dimids, &zvarid);  
   nc_def_var(ncid, "B", NC_DOUBLE, 2, dimids, &bvarid);  
   nc_def_var(ncid, "Bx", NC_DOUBLE, 2, dimids, &bxvarid);
   nc_def_var(ncid, "By", NC_DOUBLE, 2, dimids, &byvarid);
   nc_def_var(ncid, "Bz", NC_DOUBLE, 2, dimids, &bzvarid);  
   nc_def_var(ncid, "Bn", NC_DOUBLE, 2, dimids, &bnvarid);

   nc_enddef(ncid);
   nc_put_var_double(ncid, xvarid, &xsurf[0]);
   nc_put_var_double(ncid, yvarid, &ysurf[0]);
   nc_put_var_double(ncid, zvarid, &zsurf[0]);
   nc_put_var_double(ncid, bvarid, &Bmfil[0]);
   nc_put_var_double(ncid, bxvarid, &Bmfilx[0]);
   nc_put_var_double(ncid, byvarid, &Bmfily[0]);
   nc_put_var_double(ncid, bzvarid, &Bmfilz[0]);
   nc_put_var_double(ncid, bnvarid, &Bmfiln[0]);

   nc_close(ncid);
}

 
void WriteMultiFilaments(void){

   int i,j,k;
   FILE* fb;
   //fb = fopen("./outputfiles/mfil.out","w");
   fb = fopen(mfil_output,"w");
   fprintf(fb, "periods 1\n begin filament\n mirror NIL\n");
   int Nfils = Ntorfil*Nradfil;
   
   for(i=0;i<Ncoils;i++){
      for(j=0;j<Nfils;j++){
         for(k=0;k<Nseg;k++){
         fprintf(fb,"%.15f %.15f %.15f %.8f \n", *(mfilx+i*(Nseg+1)*Nfils+j*(Nseg+1)+k), \
                                                 *(mfily+i*(Nseg+1)*Nfils+j*(Nseg+1)+k), \
                                                 *(mfilz+i*(Nseg+1)*Nfils+j*(Nseg+1)+k), *(currents+i) / Nfils);     
         }
      fprintf(fb,"%.15f %.15f %.15f %.8f Mod %d %d\n", *(mfilx+i*(Nseg+1)*Nfils+j*(Nseg+1)), \
                                                       *(mfily+i*(Nseg+1)*Nfils+j*(Nseg+1)), \
                                                       *(mfilz+i*(Nseg+1)*Nfils+j*(Nseg+1)), *(currents+i) / Nfils, 1,1);   
      }
   }
   fprintf(fb,"end");
}


//void WriteFiniteBuild(void){}



