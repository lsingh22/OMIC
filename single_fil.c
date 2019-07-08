
#include "single_fil.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// GLOBALS SCOPED IN SOURCE FILE

int Ncoils;
int Nseg;
int Nfp;
int isSym;
int NFcoil; //TODO
double* coilspace;
double* currents;
double* cx;
double* cy;
double* cz;
double* sfilx;
double* sfily;
double* sfilz;



void UnpackSingleFilaments(void){

   int ind = 0;
   int ind_arr[Ncoils];
   int i,j,k;
   double x,y,z;
   double theta; //TODO: insert theta directly to improve speed
   double pi = M_PI;
   currents = (double*) malloc(Ncoils*sizeof(double));
   cx = (double*) malloc(Ncoils*sizeof(double));
   cy = (double*) malloc(Ncoils*sizeof(double));
   cz = (double*) malloc(Ncoils*sizeof(double));

   double* coilamps; coilamps = malloc(Ncoils*(NFcoil+3)*6*sizeof(double));
   
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

   sfilx = (double*) malloc(Ncoils*Nseg*sizeof(double));
   sfily = (double*) malloc(Ncoils*Nseg*sizeof(double));
   sfilz = (double*) malloc(Ncoils*Nseg*sizeof(double));

   for(i=0;i<Ncoils;i++){   //TODO: Replace theta for speed
      for(j=0;j<Nseg;j++){
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
        
         *(sfilx + i*Nseg + j ) = x;
	 *(sfily + i*Nseg + j ) = y;
         *(sfilz + i*Nseg + j ) = z;
      }
   }
}

void WriteSingleFilaments(void){
  
   int i,j;
   FILE* fb;
   fb = fopen("./outputfiles/sfil.out","w");
   fprintf(fb, "periods 1\n begin filament\n mirror NIL\n");

   for(i=0;i<Ncoils;i++){
      for(j=0;j<Nseg;j++){
         fprintf(fb,"%.15f %.15f %.15f %.15f \n", *(sfilx+i*Nseg+j), *(sfily+i*Nseg+j), *(sfilz+i*Nseg+j), *(currents+i));         
         }
      fprintf(fb,"%.15f %.15f %.15f %.15f Mod %d\n", *(sfilx+i*Nseg), *(sfily+i*Nseg), *(sfilz+i*Nseg), *(currents+i), i+1);         
   }
   fprintf(fb,"end");

}

void WriteSingleFilamentsNC(void){

//TODO: Write NETCDF out for more general and faster plotting

}
