#include "solvers.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "multi_fil.h"
#include "single_fil.h"
#include "alpha.h"
#include "single_fil.h"
#include <mpi.h>

#define MY_PI 3.14159265358979323846 
//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
 
// GLOBALS SCOPED IN SOURCE FILE

double* alpampsinit; double* derivs; double* alpamps; double* derivs; double* descent_dir; double* Bmfiln; double* nsurfn; double* fbn;
double weight_comp; double nvals_scaling; double comp_penalty_init;

int nproc; int pn;
int case_alpha; int Ncoil; int Nfp; int NFalpha; int size_alpamp; double alp_const; int size_fp; int iCoil;

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
 
double SingleFieldError(void){
//----------------------------------------------------------------------------------------------------
// Returns the single-filament normal error objective function value
//----------------------------------------------------------------------------------------------------

   register int i;
   double dsfactor = 4*pow(MY_PI,2) / (Nteta*Nzeta);
   double feval = 0.0;

   for(i=0;i<size_fp;i++)
   {
      feval += (0.5)*pow(*(Bsfiln+i),2) * (1/(pow(*(Bsfilx+i),2)+pow(*(Bsfily+i),2)+pow(*(Bsfilz+i),2))) * (*(nsurfn+i))*dsfactor;
   }

   feval = feval*Nfp; //Multiply by Nfp since only calculating for one fp
   return feval;

}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
 
double MultiFieldError(void){ 
//----------------------------------------------------------------------------------------------------
// Returns the multi-filament normal error objective function value
//----------------------------------------------------------------------------------------------------

   register int i;
   double dsfactor = 4*pow(MY_PI,2) / (Nteta*Nzeta);
   double feval = 0.0;
 
   for(i=0;i<size_fp;i++)
   {
      feval += (0.5)*pow(*(Bmfiln+i),2) * (1/(pow(*(Bmfilx+i),2)+pow(*(Bmfily+i),2)+pow(*(Bmfilz+i),2))) * (*(nsurfn+i))*dsfactor;
   }

   feval = feval*Nfp;
   return feval;
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
 
double ComplexityPenalty(void){
//----------------------------------------------------------------------------------------------------
// Returns the spectral weighting objective function value
// TODO: make the name SpectralWeight per discussion on 06/23//20  
//----------------------------------------------------------------------------------------------------

   register int i,j;
   double feval = 0.0;
   double* nvals = (double*) malloc( size_alpamp*sizeof(double) );
   int alp_per_coil = 2 * NFalpha + 1;
   
   *(nvals+0) = 0.0;
 
   for(i=1;i<NFalpha+1;i++)
   {
      *(nvals+i) = (double) i;
      *(nvals+i+NFalpha) = (double) i;
   }
      
   for(i=0;i<alp_per_coil;i++)
   {
      for(j=1;j<iCoil;j++)
      {
         *(nvals+j*alp_per_coil+i) = *(nvals+i);
      }
   }
  
   for(i=0;i<size_alpamp;i++)
   {
      feval += pow(*(alpamps+i),2) * pow( *(nvals+i), nvals_scaling);
   }   
   //printf("The complexity function is %.9f\n", feval);
   return feval;
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
 
double CostFunction(int case_objfun, double fb_init) {
//----------------------------------------------------------------------------------------------------
// Returns the value of the total objective function
// case_objfun: 0-fb only  1-fb and fc //TODO: toggle by value of weights to simplify
//----------------------------------------------------------------------------------------------------
   double fb, fc;
   double feval = 0.0;

   //Calculate the objective functions
   fb = MultiFieldError();
   fc = ComplexityPenalty();  

   if(case_objfun==0)
   {   
      feval = fb;
   }
   else if(case_objfun==1)
   {  
      feval = fb / fb_init + weight_comp * fc;
   }
   
   return feval;

}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
 
double SurfaceArea(void){
//----------------------------------------------------------------------------------------------------
// Returns the surface area of the magnetic boundary
// TODO: is this correct? might be off by a factor of nfp 
//----------------------------------------------------------------------------------------------------

   register int i;
   double dsfactor = 4 * pow(MY_PI,2) / (Nteta*Nzeta);
   double area = 0.0;

   for(i=0;i<size_fp;i++)
   {
      area += *(nsurfn+i);   
   }
   area = area*dsfactor; 
 
   return area;
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
 
void Central_diff( double *dof, double fb_init ){
//----------------------------------------------------------------------------------------------------
// Calculates the numerical derivative of the objective function with respect to alpha harmonics
//----------------------------------------------------------------------------------------------------

   int i,j;
   double h = 0.000001;
   
   double* nvals = (double*) malloc( size_alpamp*sizeof(double) );
   derivs = (double*) malloc( size_alpamp*sizeof(double) );
   
   double minus_bn;
   double plus_bn;
   double init_bn = 0.0;
   int alp_per_coil = 2*NFalpha+1;

   if(init_bn > 100000){return;} 

   for(i=0;i<size_alpamp;i++)
   {
      *(alpamps+i) = dof[i];
   }

   if(case_objfun==1)
   {
      *(nvals+0) = 0.0;
      for(i=1;i<NFalpha;i++)
      {
         *(nvals+i) = (double) i;
         *(nvals+i+NFalpha) = (double) i;
      }
      
      for(i=0;i<alp_per_coil;i++)
      {
         for(j=1;j<iCoil;j++)
         {
            *(nvals+j*alp_per_coil+i) = *(nvals+i);
         }
      }
   }
   
   if(pn==0){init_bn = MultiFieldError();}
   
   for(i=0;i<size_alpamp;i++)
   {
      minus_bn = 0.0;
      plus_bn = 0.0;
      // Move alpha positive and redo x,y,z coil calc
      
      *(alpamps+i) += h;
      CalculateMultiFilaments();
      MultifilamentField();
      if(nproc > 1){GatherFieldData();}

      if(pn==0){plus_bn = MultiFieldError();}

      *(alpamps+i) -= 2*h;
      CalculateMultiFilaments();
      MultifilamentField();
      if(nproc > 1){GatherFieldData();}

      if(pn==0){minus_bn = MultiFieldError();} 

      *(alpamps+i) += h;
      
      //if(pn==0)
      //{
      //   printf("The value of multi_errror_init is %.9f and weight_comp is %.9f.\n", fb_init, weight_comp);
      //   printf("The value of plus_bn is %.9f and minus_bn is %.9f.\n", plus_bn, minus_bn);
      //   printf("The value of alp_amps+i is %.9f and nvals+i is %.9f.\n", *(alpamps+i), *(nvals+i));
      //}
      
      if(pn==0)
      {
         if(case_objfun==1)
         {
            *(derivs+i) = (1/fb_init)* ((plus_bn-minus_bn)/(2*h)) + 2 * weight_comp * *(alpamps+i) * pow(*(nvals+i),nvals_scaling);
         }
         else if(case_objfun==0)
         {
            *(derivs+i) = (plus_bn-minus_bn)/(2*h);
         }
      }   
   }
}
//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
 
void Steepest_descent( void ){
//----------------------------------------------------------------------------------------------------
// Calculates the negative gradient of the objective function
// TODO: this should just be flagged in central_diff, no need for another function for adding a - sign
//----------------------------------------------------------------------------------------------------
   
   register int i;
   
   descent_dir = (double*) malloc( size_alpamp*sizeof(double) );
  
   if(pn==0)
   {     
      for(i=0;i<size_alpamp;i++)
      {
            *(descent_dir+i) = -1.0 * (*(derivs+i));
      //printf("The descent of amp %d is %.12f \n",i,*(descent_dir+i));
      }    
   }
}
//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
 
void Forward_track(double fb_init ){
//----------------------------------------------------------------------------------------------------
// Optimizes alpha harmonics using a forward-tracking line search
// Maybe just look at k+step*(1/2) where k+1 is the exit iteration
//----------------------------------------------------------------------------------------------------

   int j, k = 0;

   double step = .000000001; 
   double init_bn = 0;
   double fb_now=0.0, fc_now=0.0;
   double search_bn;
   double hold_bn; 
 
   CalculateMultiFilaments();
   MultifilamentField();
   if(nproc > 1){GatherFieldData();}
   
   if(pn==0)
   {
      init_bn = CostFunction(case_objfun,fb_init);
      hold_bn = init_bn;
      search_bn = 0.0;
   }

   MPI_Bcast(&hold_bn, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
   MPI_Bcast(&search_bn, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
 
   while( hold_bn - search_bn > 0.0 )
   {
      if(pn==0)
      {
         if(k==0)
         {
            search_bn = init_bn;
         }

         hold_bn = search_bn;
         search_bn = 0.0;
         printf("Step size is %f\n", 1000.0*step);
      
         for(j=0;j<size_alpamp;j++)
         {   
         //printf("NF:   %dAlphas   %f\n",j,*(alpamps+j));
            *(alpamps+j) += step*(*(descent_dir+j));      
         }
      }
      
      MPI_Bcast(alpamps,size_alpamp, MPI_DOUBLE, 0, MPI_COMM_WORLD);    
      CalculateMultiFilaments();
      MultifilamentField();
      if(nproc > 1){GatherFieldData();}

      if(pn==0)
      {
         search_bn = CostFunction(case_objfun,fb_init);
         fb_now = MultiFieldError();
         fc_now = ComplexityPenalty();
 
         printf("\nTotal cost function value, tracking iter: %.9f   %d\n", search_bn, k);
         printf("The fB value is: %.9f   \n", fb_now);
         printf("The fC value is: %.9f   \n", fc_now);
      } 
      
      step = step * 2.0;
      k++;
      
      // These two lines allow each processor to exit the loop
      MPI_Bcast(&hold_bn, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
      MPI_Bcast(&search_bn, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
   }
   
   if(pn==0)
   {
      for(j=0;j<size_alpamp;j++)
      {
         *(alpamps+j) -= step*(*(descent_dir+j))/2.0;
      }        
   }

   MPI_Bcast(alpamps,size_alpamp, MPI_DOUBLE, 0, MPI_COMM_WORLD);  
   
   CalculateMultiFilaments();
   MultifilamentField();
   if(nproc > 1){GatherFieldData();}
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
