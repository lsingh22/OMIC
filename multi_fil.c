#include "multi_fil.h"
#include "alpha.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "bfield.h"
#include "bfield_gpu.cuh"
#include <omp.h>
#include <mpi.h>

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

//GLOBALS SCOPED IN SOURCE FILE

int Nthreads; double len_rad; double len_tor; int Nseg; int Nradfil; int Ntorfil; int Nfils; int Nfp;

int Nzeta; int Nteta; int nproc; int pn; int Ncoil; int iCoil; int size_fp; int Ns;

double* tx; double* ty; double* tz; double* nx; double* ny; double* nz; double* bx; double* by; double* bz; 

double* alp; double* cx; double* cy; double* cz; double* sfilx; double* sfily; double* sfilz;

double* mfilx; double* mfily; double* mfilz; double* finx; double* finy; double* finz;

double* Bmfilx; double* Bmfily; double* Bmfilz; double* Bmfil; double* Bmfiln;

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

void CalculateBuildDirections(void){
//----------------------------------------------------------------------------------------------------
// Determine the coil-centroid local coordinate frame from single-filament Fourier representation   
//----------------------------------------------------------------------------------------------------

   double norm; //used for finding unit vectors
   double theta, x,y,z;
   double pi = 3.14159265358979323846;

   tx = (double *) malloc((Ncoil)*(Nseg+1)*sizeof(double));
   ty = (double *) malloc((Ncoil)*(Nseg+1)*sizeof(double));
   tz = (double *) malloc((Ncoil)*(Nseg+1)*sizeof(double));
   nx = (double *) malloc((Ncoil)*(Nseg+1)*sizeof(double));
   ny = (double *) malloc((Ncoil)*(Nseg+1)*sizeof(double));
   nz = (double *) malloc((Ncoil)*(Nseg+1)*sizeof(double));
   bx = (double *) malloc((Ncoil)*(Nseg+1)*sizeof(double));
   by = (double *) malloc((Ncoil)*(Nseg+1)*sizeof(double));
   bz = (double *) malloc((Ncoil)*(Nseg+1)*sizeof(double));
      
   register int i,j;
   int k; 
   double dot; //used to store dot products

   //Calculate unit tangent vector by taking the derivative of the single-filament position vector
   for(i=0;i<Ncoil;i++)
   {
      for(j=0;j<Nseg+1;j++)
      { //consider a coil as Nseg unique points (the Nseg-th point is just the starting point)
         theta = ((2*pi)/Nseg)*j; 
         x=0; y=0; z=0; norm=0;
         for(k=0;k<NFcoil+1;k++)
         { //add the cosine components
           x = x - k*coilamps[ ind_arr[i] + k ]*sin(k*theta);
           y = y - k*coilamps[ ind_arr[i] + k + 2*NFcoil + 1 ]*sin(k*theta);
           z = z - k*coilamps[ ind_arr[i] + k + 4*NFcoil + 2 ]*sin(k*theta);
         }
         for(k=1;k<NFcoil+1;k++)
         { //add the sine components
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

   //Calculate vector pointing from coil centroid to point on coil for each coil
   //In the JPP 2020 publication, this is the delta vector used to determine coil-centroid normal vector
   double* sfilxa; sfilxa = (double *) malloc((Ncoil)*(Nseg+1)*sizeof(double));
   double* sfilya; sfilya = (double *) malloc((Ncoil)*(Nseg+1)*sizeof(double));
   double* sfilza; sfilza = (double *) malloc((Ncoil)*(Nseg+1)*sizeof(double));

   for(i=0;i<Ncoil;i++)
   {
      for(j=0;j<Nseg+1;j++)
      {
         *(sfilxa + i*(Nseg+1) + j ) = *(sfilx + i*(Nseg+1) + j ) - *(cx + i);
         *(sfilya + i*(Nseg+1) + j ) = *(sfily + i*(Nseg+1) + j ) - *(cy + i);
         *(sfilza + i*(Nseg+1) + j ) = *(sfilz + i*(Nseg+1) + j ) - *(cz + i);
	 //printf("%f %f %f\n", *(sfilxa + i*Nseg + j ), *(sfilya + i*Nseg + j ), *(sfilza + i*Nseg + j ));
      }
   }

   for(i=0;i<Ncoil*(Nseg+1);i++)
   {
      x=0;y=0;z=0;
      //Dot product of delta and tangent vector
      dot = *(sfilxa + i) * *(tx + i) + *(sfilya + i) * *(ty + i) + *(sfilza + i) * *(tz + i);
      //Components of normal vector before making it a unit vector
      x = *(sfilxa + i) - dot* *(tx + i);
      y = *(sfilya + i) - dot* *(ty + i);
      z = *(sfilza + i) - dot* *(tz + i);
      norm = sqrt(x*x + y*y + z*z); //normal vector magnitude
      //Unit normal vector components
      *(nx+i) = x/norm;
      *(ny+i) = y/norm;
      *(nz+i) = z/norm;
      //Unit binormal found via croos product b = t X n
      *(bx+i) = *(ty+i) * *(nz+i) - *(tz+i) * *(ny+i);
      *(by+i) = *(tz+i) * *(nx+i) - *(tx+i) * *(nz+i);
      *(bz+i) = *(tx+i) * *(ny+i) - *(ty+i) * *(nx+i);
   }
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

void CalculateMultiFilaments(void){
//----------------------------------------------------------------------------------------------------
// Calculate multi-filament coils using pre-computed single-filament coils, coil-centroid frame, and 
// alpha optimization function for each coil. 
// A rectangular winding pack is assumed
// Filaments are spaced evenly within the winding pack  
//----------------------------------------------------------------------------------------------------

   Unpack_alpha();
   CalculateBuildDirections();
   
   register int i,j,k,l;
   int ip; //used for rotations if periodicity is enforced
   
   //Half the spacing between filaments in radial and toroidal directions
   //TODO: Check this for spacing considerations with different winding packs
   double gridlen = len_tor / (1*Nradfil);
   double gridwid = len_rad / (4*Ntorfil);
   //for locating surface interp points
   double hlen_rad = len_rad / 2;
   double hlen_tor = len_tor / 2;
   
   //The rotated normal and binormal vectors
   double* nxa;
   double* nya;
   double* nza;
   double* bxa;
   double* bya;
   double* bza;

   nxa = (double *) malloc((Ncoil)*(Nseg+1)*sizeof(double));
   nya = (double *) malloc((Ncoil)*(Nseg+1)*sizeof(double));
   nza = (double *) malloc((Ncoil)*(Nseg+1)*sizeof(double));
   bxa = (double *) malloc((Ncoil)*(Nseg+1)*sizeof(double));
   bya = (double *) malloc((Ncoil)*(Nseg+1)*sizeof(double));
   bza = (double *) malloc((Ncoil)*(Nseg+1)*sizeof(double));
    
   mfilx = (double*) malloc(Ncoil*Nfils*(Nseg+1)*sizeof(double));
   mfily = (double*) malloc(Ncoil*Nfils*(Nseg+1)*sizeof(double));
   mfilz = (double*) malloc(Ncoil*Nfils*(Nseg+1)*sizeof(double));
 
   ffilx = (double*) malloc(Ncoil*Nfils*(Nseg+1)*5*sizeof(double));
   ffily = (double*) malloc(Ncoil*Nfils*(Nseg+1)*5*sizeof(double));
   ffilz = (double*) malloc(Ncoil*Nfils*(Nseg+1)*5*sizeof(double));
   
   // Rotate the local basis using alpha parameter
   for(i=0;i<Ncoil*(Nseg+1);i++)
   {
      *(nxa+i) = *(nx+i)*cos(*(alp+i)) + *(bx+i)*sin(*(alp+i));
      *(nya+i) = *(ny+i)*cos(*(alp+i)) + *(by+i)*sin(*(alp+i));
      *(nza+i) = *(nz+i)*cos(*(alp+i)) + *(bz+i)*sin(*(alp+i));

      *(bxa+i) = -*(nx+i)*sin(*(alp+i)) + *(bx+i)*cos(*(alp+i));
      *(bya+i) = -*(ny+i)*sin(*(alp+i)) + *(by+i)*cos(*(alp+i));
      *(bza+i) = -*(nz+i)*sin(*(alp+i)) + *(bz+i)*cos(*(alp+i));
   
      for(j=0;j<5;j++)
      {
         *(ffilx + 5*i) =     *(sfilx+i) + hlen_rad * *(nxa+i) + hlen_tor * *(bxa+i); //top right x-coord
         *(ffilx + 5*i + 1) = *(sfilx+i) - hlen_rad * *(nxa+i) + hlen_tor * *(bxa+i); //top left
         *(ffilx + 5*i + 2) = *(sfilx+i) - hlen_rad * *(nxa+i) - hlen_tor * *(bxa+i); //bottom left
         *(ffilx + 5*i + 3) = *(sfilx+i) + hlen_rad * *(nxa+i) - hlen_tor * *(bxa+i); //bottom right
         *(ffilx + 5*i + 4) = *(sfilx+i) + hlen_rad * *(nxa+i) + hlen_tor * *(bxa+i); //top right again

         *(ffily + 5*i) =     *(sfily+i) + hlen_rad * *(nya+i) + hlen_tor * *(bya+i); //analogous y-coords
         *(ffily + 5*i + 1) = *(sfily+i) - hlen_rad * *(nya+i) + hlen_tor * *(bya+i);
         *(ffily + 5*i + 2) = *(sfily+i) - hlen_rad * *(nya+i) - hlen_tor * *(bya+i);
         *(ffily + 5*i + 3) = *(sfily+i) + hlen_rad * *(nya+i) - hlen_tor * *(bya+i);
         *(ffily + 5*i + 4) = *(sfily+i) + hlen_rad * *(nya+i) + hlen_tor * *(bya+i);

         *(ffilz + 5*i)     = *(sfilz+i) + hlen_rad * *(nza+i) + hlen_tor * *(bza+i); //analogous z-coords
         *(ffilz + 5*i + 1) = *(sfilz+i) - hlen_rad * *(nza+i) + hlen_tor * *(bza+i);
         *(ffilz + 5*i + 2) = *(sfilz+i) - hlen_rad * *(nza+i) - hlen_tor * *(bza+i);
         *(ffilz + 5*i + 3) = *(sfilz+i) + hlen_rad * *(nza+i) - hlen_tor * *(bza+i);
         *(ffilz + 5*i + 4) = *(sfilz+i) + hlen_rad * *(nza+i) + hlen_tor * *(bza+i);	 
      }    
   }

   //Calculates the multi-filament coils
   //Only one field period if enforcing periodicity
   
   for(i=0;i<iCoil;i++)
   {
      for(j=0;j<Ntorfil;j++)
      {
         for(k=0;k<Nradfil;k++)
         {
            for(l=0;l<Nseg;l++)
            {
            *(mfilx + i*Nfils*(Nseg+1) + j*(Nseg+1)*Nradfil + k*(Nseg+1) + l) = *(sfilx +i*(Nseg+1) + l) \
                                       + (gridlen* (-(Nradfil-1)+2*k))* *(nxa + i*(Nseg+1) + l) + (gridwid* (-(Ntorfil-1)+2*j))* *(bxa + i*(Nseg+1) + l);
            *(mfily + i*Nfils*(Nseg+1) + j*(Nseg+1)*Nradfil + k*(Nseg+1) + l) = *(sfily +i*(Nseg+1) + l) \
				       + (gridlen* (-(Nradfil-1)+2*k))* *(nya + i*(Nseg+1) + l) + (gridwid* (-(Ntorfil-1)+2*j))* *(bya + i*(Nseg+1) + l);
	    *(mfilz + i*Nfils*(Nseg+1) + j*(Nseg+1)*Nradfil + k*(Nseg+1) + l) = *(sfilz +i*(Nseg+1) + l) \
				       + (gridlen* (-(Nradfil-1)+2*k))* *(nza + i*(Nseg+1) + l) + (gridwid* (-(Ntorfil-1)+2*j))* *(bza + i*(Nseg+1) + l);
	    		} 
 	    *(mfilx + i*Nfils*(Nseg+1) + j*(Nseg+1)*Nradfil + k*(Nseg+1) + Nseg ) = *(mfilx + i*Nfils*(Nseg+1) + j*(Nseg+1)*Nradfil + k*(Nseg+1) );
            *(mfily + i*Nfils*(Nseg+1) + j*(Nseg+1)*Nradfil + k*(Nseg+1) + Nseg ) = *(mfily + i*Nfils*(Nseg+1) + j*(Nseg+1)*Nradfil + k*(Nseg+1) );
	    *(mfilz + i*Nfils*(Nseg+1) + j*(Nseg+1)*Nradfil + k*(Nseg+1) + Nseg ) = *(mfilz + i*Nfils*(Nseg+1) + j*(Nseg+1)*Nradfil + k*(Nseg+1) );         
         }
      }
   }
   
   //Calculates periodic coils when enforced
   //TODO: this should maybe be flagged or at least print something like "No periodicity enforced -- calculating all coils."
   for(ip=2;ip<Nfp+1;ip++){
      for(j=0;j<iCoil*Nfils*(Nseg+1);j++){
         *(mfilx + (ip-1)*(iCoil*Nfils*(Nseg+1))+j) = *(mfilx+j)*cosnfp(ip) - *(mfily+j)*sinnfp(ip);
         *(mfily + (ip-1)*(iCoil*Nfils*(Nseg+1))+j) = *(mfilx+j)*sinnfp(ip) + *(mfily+j)*cosnfp(ip);
         *(mfilz + (ip-1)*(iCoil*Nfils*(Nseg+1))+j) = *(mfilz+j);
      }
   }
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

void MultifilamentField(void){
//----------------------------------------------------------------------------------------------------
// Calculate the field due to multi-filament coils
// Periodicity is assumed  
//----------------------------------------------------------------------------------------------------
  
   int i, j, nave;
   float start, end, total=0.0; 

   Bmfilx = (double*) malloc(Nzeta * Nteta * sizeof(double));
   Bmfily = (double*) malloc(Nzeta * Nteta * sizeof(double));
   Bmfilz = (double*) malloc(Nzeta * Nteta * sizeof(double));
   Bmfiln = (double*) malloc(Nzeta * Nteta * sizeof(double));
    Bmfil = (double*) malloc(Nzeta * Nteta * sizeof(double));

   // MPI rank chunks
   int first = startind[pn];
   int last  = endind[pn] + 1;

	// Number of times to run integration for averaging (759)
	nave = 5;

	for(j = 0; j < nave; j++) {   
	   start = MPI_Wtime();	
	
   	// Calculate field using serial computation
  		//CalculateFieldSerial(first, last);

		// Calculate field using GPU
		CalculateFieldParallelGPU(first, last);
		
	   end = MPI_Wtime();
		total += (end - start);
	}

	// Print the average time
	if(pn==0){printf("%f\n", total / nave);}
   
   //Calculate |B| and B*n on one field period
   if(nproc == 1) {
      for(i = 0; i < size_fp; i++) {
         Bmfiln[i] = Bmfilx[i] * nsurfx[i] + Bmfily[i] * nsurfy[i] + Bmfilz[i] * nsurfz[i];  
         Bmfil[i] =  sqrt( pow(Bmfilx[i], 2) + pow(Bmfily[i], 2) + pow(Bmfilz[i], 2) ); 
      }
      
		// Reflect magnetics to symmetric field periods
		ReflectFieldPeriod();
   }
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

void ReflectFieldPeriod(void) {
//----------------------------------------------------------------------------------------------------
// Write magnetics to symmetric field periods
// Periodicity assumed, but not stellarator symmetry
//----------------------------------------------------------------------------------------------------

	register int i, ip;
   double rot_cos, rot_sin;

   // Reflect magnetics to all field periods
   for(ip = 2; ip < Nfp + 1; ip++) {
      rot_cos = cosnfp(ip);
      rot_sin = sinnfp(ip);
       
      for(i = 0; i < size_fp; i++) {
         Bmfil [(ip-1) * size_fp +i] = Bmfil[i];       
         Bmfiln[(ip-1) * size_fp +i] = Bmfiln[i];
         Bmfilx[(ip-1) * size_fp +i] = Bmfilx[i] * rot_cos - Bmfily[i] * rot_sin;
         Bmfily[(ip-1) * size_fp +i] = Bmfilx[i] * rot_sin + Bmfily[i] * rot_cos;
         Bmfilz[(ip-1) * size_fp +i] = Bmfilz[i];
      }
   }
}

void GatherFieldData(void){
//----------------------------------------------------------------------------------------------------
// Store MPI data on the field due to multi-filament coils
// Gathers data from each node and stores it in the global Bx,By,Bz,Bn,B
//----------------------------------------------------------------------------------------------------
   
   int pri, sizepn, first; 
   register int ip, i, j, k;
   double* temp_Bmfilx;
   double* temp_Bmfily;
   double* temp_Bmfilz;
   double rot_cos, rot_sin;
   MPI_Status status; 

   MPI_Barrier(MPI_COMM_WORLD);   

   if(pn==0) //Parent process
   {
      for(int i=2;i<nproc+1;i++)
      {
         pri = i - 1;
         sizepn = 1 + *(endind+pri) - *(startind+pri);

         temp_Bmfilx = (double*) malloc(sizepn*sizeof(double));
         temp_Bmfily = (double*) malloc(sizepn*sizeof(double));
         temp_Bmfilz = (double*) malloc(sizepn*sizeof(double));
 
         first = *(startind+pri);

         MPI_Recv(temp_Bmfilx, sizepn, MPI_DOUBLE, pri, 10+100*pri, MPI_COMM_WORLD, &status);
         for(j=0;j<sizepn; *(Bmfilx+first+j) = *(temp_Bmfilx+j),j++);

         MPI_Recv(temp_Bmfily, sizepn, MPI_DOUBLE, pri, 11+100*pri, MPI_COMM_WORLD, &status);
         for(j=0;j<sizepn; *(Bmfily+first+j) = *(temp_Bmfily+j),j++);

         MPI_Recv(temp_Bmfilz, sizepn, MPI_DOUBLE, pri, 12+100*pri, MPI_COMM_WORLD, &status);
         for(j=0;j<sizepn; *(Bmfilz+first+j) = *(temp_Bmfilz+j),j++);

         free(temp_Bmfilx);
         free(temp_Bmfily);
         free(temp_Bmfilz);
      }
      
      //Reflect to each period of the plasma boundary
      for(ip=2;ip<Nfp+1;ip++)
      {
         rot_cos = cosnfp(ip);
         rot_sin = sinnfp(ip);
 
         for(k=0;k<size_fp;k++)
         {
            *(Bmfilx + (ip-1)*size_fp+k) = *(Bmfilx+k) * rot_cos - * (Bmfily+k) * rot_sin;
            *(Bmfily + (ip-1)*size_fp+k) = *(Bmfilx+k) * rot_sin + * (Bmfily+k) * rot_cos;
            *(Bmfilz + (ip-1)*size_fp+k) = *(Bmfilz+k);
         }
      }
      
      //Finally, calculate the field and normal field strength.
      for(i=0;i<size_fp;i++)
      {
         *(Bmfiln+i) = *(Bmfilx+i) * *(nsurfx+i) + *(Bmfily+i) * *(nsurfy+i) + \
                    *(Bmfilz+i) * *(nsurfz+i); 
         *(Bmfil+i) = sqrt( pow(*(Bmfilx+i),2) + pow(*(Bmfily+i),2) + pow(*(Bmfilz+i),2) ); 
      }

   }else //Not parent process
   {
      sizepn = 1 + *(endind+pn) - *(startind+pn);  
         
      temp_Bmfilx = (double*) malloc(sizepn*sizeof(double));
      temp_Bmfily = (double*) malloc(sizepn*sizeof(double));
      temp_Bmfilz = (double*) malloc(sizepn*sizeof(double));

      first = *(startind+pn);

      for(j=0;j<sizepn; *(temp_Bmfilx+j) = *(Bmfilx + first + j), j++);
      MPI_Send(temp_Bmfilx, sizepn, MPI_DOUBLE, 0, 10+100*pn, MPI_COMM_WORLD);
      
      for(j=0;j<sizepn; *(temp_Bmfily+j) = *(Bmfily + first + j), j++);
      MPI_Send(temp_Bmfily, sizepn, MPI_DOUBLE, 0, 11+100*pn, MPI_COMM_WORLD);
 
      for(j=0;j<sizepn; *(temp_Bmfilz+j) = *(Bmfilz + first + j), j++);
      MPI_Send(temp_Bmfilz, sizepn, MPI_DOUBLE, 0, 12+100*pn, MPI_COMM_WORLD);

      free(temp_Bmfilx);
      free(temp_Bmfily);
      free(temp_Bmfilz);
   }  

   MPI_Barrier(MPI_COMM_WORLD);

}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

void WriteMultiFilaments(void){
//----------------------------------------------------------------------------------------------------
// Writes a coils file containing multi-filament coordinates to a txt file
//----------------------------------------------------------------------------------------------------

   register int i,j,k;
   FILE* fb;
   //fb = fopen("./outputfiles/mfil.out","w");
   fb = fopen(mfil_output,"w");
   fprintf(fb, "periods %d \nbegin filament\nmirror NIL\n", Nfp);
   int Nfils = Ntorfil*Nradfil;
   
   for(i=0;i<Ncoil;i++)
   {
      for(j=0;j<Nfils;j++)
      {
         for(k=0;k<Nseg;k++)
         {
            fprintf(fb,"%.15f %.15f %.15f %.15f \n", *(mfilx+i*(Nseg+1)*Nfils+j*(Nseg+1)+k), \
                                                     *(mfily+i*(Nseg+1)*Nfils+j*(Nseg+1)+k), \
                                                     *(mfilz+i*(Nseg+1)*Nfils+j*(Nseg+1)+k), *(currents+i) / Nfils);     
         }
         fprintf(fb,"%.15f %.15f %.15f %.15f 1 Mod\n", *(mfilx+i*(Nseg+1)*Nfils+j*(Nseg+1)), \
                                                       *(mfily+i*(Nseg+1)*Nfils+j*(Nseg+1)), \
                                                       *(mfilz+i*(Nseg+1)*Nfils+j*(Nseg+1)), 0.0 );   
      }
   }  
   fprintf(fb,"end");
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

