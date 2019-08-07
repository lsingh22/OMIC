
#include "single_fil.h"
#include "bfield.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <unistd.h>

// GLOBALS SCOPED IN SOURCE FILE

int Ncoils;
int Nseg;
int Nfp;
int isSym;
int NFcoil; //TODO

int Nzeta;
int Nteta;

double* coilspace;
double* currents;
double* cx;
double* cy;
double* cz;
double* sfilx;
double* sfily;
double* sfilz;
double* coilamps;
int* ind_arr;

double* Bsfilx;
double* Bsfily;
double* Bsfilz;
double* Bsfiln;
double* Bsfil;
 
double* nsurfx;
double* nsurfy;
double* nsurfz;

double* xsurf;
double* ysurf;
double* zsurf;


void UnpackSingleFilaments(void){

   int ind = 0;
   ind_arr = malloc(Ncoils*sizeof(int));
   int i,j,k;
   double x,y,z;
   double theta; //TODO: insert theta directly to improve speed
   double pi = M_PI;
   currents = (double*) malloc(Ncoils*sizeof(double));
   cx = (double*) malloc(Ncoils*sizeof(double));
   cy = (double*) malloc(Ncoils*sizeof(double));
   cz = (double*) malloc(Ncoils*sizeof(double));

   coilamps = malloc(Ncoils*(NFcoil+3)*6*sizeof(double));
   
   for(i=0;i<Ncoils;i++){
      if (i==0) 
      {
        *(currents + i) = coilspace[0];
      }
      else
      {
        *(currents + i) = coilspace[ind+i];
      }
      for(j=0;j<6*(NFcoil+1)-3;j++){
         *(coilamps + ind + j) = coilspace[ind + j + i + 1 ];
      }
      ind_arr[i] = ind;
      ind = ind + NFcoil*6+3;
   }
   //Store centroids from coilamps array//

   for(i=0;i<Ncoils;i++){
        cx[i] = coilamps[ind_arr[i]];
        cy[i] = coilamps[ind_arr[i] + 2*(NFcoil + 1)-1];
        cz[i] = coilamps[ind_arr[i] + 4*(NFcoil + 1)-2];
   } 

   sfilx = (double*) malloc(Ncoils*(Nseg+1)*sizeof(double));
   sfily = (double*) malloc(Ncoils*(Nseg+1)*sizeof(double));
   sfilz = (double*) malloc(Ncoils*(Nseg+1)*sizeof(double));

   for(i=0;i<Ncoils;i++){   
      for(j=0;j<Nseg+1;j++){
         theta = ((2*pi)/Nseg)*j;
	 x=0;y=0;z=0;
         for(k=0;k<NFcoil+1;k++){ //add the cosine components
           x = x + coilamps[ ind_arr[i] + k ]*cos(k*theta);
           y = y + coilamps[ ind_arr[i] + k + 2*NFcoil + 1 ]*cos(k*theta);
           z = z + coilamps[ ind_arr[i] + k + 4*NFcoil + 2 ]*cos(k*theta);                
         }
         for(k=1;k<NFcoil+1;k++){ //add the sine components
           x = x + coilamps[ ind_arr[i] +   NFcoil + 0 + k ]*sin(k*theta);
           y = y + coilamps[ ind_arr[i] + 3*NFcoil + 1 + k ]*sin(k*theta);
           z = z + coilamps[ ind_arr[i] + 5*NFcoil + 2 + k ]*sin(k*theta);
         }
        
         *(sfilx + i*(Nseg+1) + j ) = x;
	 *(sfily + i*(Nseg+1) + j ) = y;
         *(sfilz + i*(Nseg+1) + j ) = z;
      }
   }
}

void SingleFilField(void){
  
   Bsfilx = (double*) malloc(Nzeta*Nteta*sizeof(double));
   Bsfily = (double*) malloc(Nzeta*Nteta*sizeof(double));
   Bsfilz = (double*) malloc(Nzeta*Nteta*sizeof(double));
   Bsfiln = (double*) malloc(Nzeta*Nteta*sizeof(double));
    Bsfil = (double*) malloc(Nzeta*Nteta*sizeof(double));

   int i;

   double timeBfield;
   double startBfield, endBfield;

   omp_set_num_threads(Nthreads);

   startBfield = omp_get_wtime(); //clock();
   // for(i=0;i<size_surf*size_surf;i++){
   //OpenMP parallelization 
   #pragma omp parallel for
   for(i=0;i<Nzeta*Nteta;i++){

      CalculateSingleField( *(xsurf+i), *(ysurf+i), *(zsurf+i), \
                            Bsfilx+i, Bsfily+i, Bsfilz+i );   

      *(Bsfiln+i) = *(Bsfilx+i) * *(nsurfx+i) + *(Bsfily+i) * *(nsurfy+i) + \
                    *(Bsfilz+i) * *(nsurfz+i);  
      *(Bsfil+i) = sqrt( pow(*(Bsfilx+i),2) + pow(*(Bsfily+i),2) + pow(*(Bsfilz+i),2) ); 
   }

   endBfield = omp_get_wtime(); //clock();
   printf("\nTotal time of single fil field calculation: %f\n\n", endBfield-startBfield);//timeBfield);  
}

#define SFILB_FILE_NAME "./outputfiles/sfilB.nc"
   
void WriteSingleB(void){
   //Write to NC
   int ncid, xvarid, yvarid, zvarid, bvarid, xdimid, ydimid;
   int bxvarid, byvarid, bzvarid;
   int dimids[2];
   
   nc_create(SFILB_FILE_NAME, NC_CLOBBER, &ncid); 
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
 

   nc_enddef(ncid);
   nc_put_var_double(ncid, xvarid, &xsurf[0]);
   nc_put_var_double(ncid, yvarid, &ysurf[0]);
   nc_put_var_double(ncid, zvarid, &zsurf[0]);
   nc_put_var_double(ncid, bvarid, &Bsfil[0]);
   nc_put_var_double(ncid, bxvarid, &Bsfilx[0]);
   nc_put_var_double(ncid, byvarid, &Bsfily[0]);
   nc_put_var_double(ncid, bzvarid, &Bsfilz[0]);
 
   nc_close(ncid);
}


void WriteSingleFilaments(void){
  
   int i,j;
   FILE* fb;
   fb = fopen("./outputfiles/sfil.out","w");
   fprintf(fb, "periods 1\n begin filament\n mirror NIL\n");

   for(i=0;i<Ncoils;i++){
      for(j=0;j<Nseg;j++){
         fprintf(fb,"%.15f %.15f %.15f %.15f \n", *(sfilx+i*(Nseg+1)+j), *(sfily+i*(Nseg+1)+j), *(sfilz+i*(Nseg+1)+j), *(currents+i));         
         }
      fprintf(fb,"%.15f %.15f %.15f %.15f Mod %d\n", *(sfilx+i*(Nseg+1)), *(sfily+i*(Nseg+1)), *(sfilz+i*(Nseg+1)), *(currents+i), i+1);         
   }
   fprintf(fb,"end");

}

void WriteSingleFilamentsNC(void){

//TODO: Write NETCDF out for more general and faster plotting

}
