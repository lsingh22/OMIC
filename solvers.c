//This is for initializing user inputs

#include "solvers.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "multi_fil.h"
#include "single_fil.h"
#include "alpha.h"
// GLOBALS SCOPED IN SOURCE FILE

double* alpampsinit;
double* derivs;
double* alpamps;
double* derivs;
double* descent_dir;
double* Bmfiln;
double* nsurfn;
double* fbn;

int case_alpha;
int Ncoils;
int Nfp;
int NFalpha;
int size_alpamp;
double alp_const;

//TODO: In the future, will need to change indexing if want NFalpha to differ for each coil

double CostFunction(int case_opt, double* dof){

int iCoils = Ncoils / Nfp;
int size_fp = Nteta*Nzeta / Nfp;
int size_alpamp = iCoils*(2*NFalpha+1);
int i;
double dsfactor = 4*pow(M_PI,2) / (Nteta*Nzeta);
double feval = 0.0;

   if(case_opt==0)
   {
      for(i=0;i<size_alpamp;i++)
      {
         *(alpamps+i) = dof[i];
      }
      CalculateMultiFilaments();
      MultiFilFieldSym();
      
      for(i=0;i<size_fp;i++){
         feval += (0.5) * (*(Bmfiln+i) * *(Bmfiln+i))*(*(nsurfn+i))*dsfactor;
      } 
   }
 
   return feval;
}


double SingleFieldError(void){

int iCoils = Ncoils / Nfp;
int size_fp = Nteta*Nzeta / Nfp;
int i;
double dsfactor = 4*pow(M_PI,2) / (Nteta*Nzeta);
double feval = 0.0;

UnpackSingleFilaments();
SingleFilField();
    
   for(i=0;i<size_fp;i++){
      feval += (0.5)*( *(fbn+i)* *(fbn+i) )* (*(nsurfn+i))*dsfactor;   
   } 
return feval;
}
 

void Central_diff( double *dof ){
   
   int iCoils = Ncoils / Nfp;
   int size_alpamp = iCoils*(2*NFalpha+1);   
   //int iCoils = Ncoils;
   int size_fp = Nteta*Nzeta / Nfp;
   int i,j;

   derivs = (double*) malloc( size_alpamp*sizeof(double) );
   double h = 1e-13;
   
   double minus_bn;
   double plus_bn;
   double init_bn = 0.0;

   for(i=0;i<size_alpamp;i++)
   {
      *(alpamps+i) = dof[i];
   }

   
   MultiFilFieldSym(); // Calculates B dot n for the multifilaments might not be needed here 
   init_bn = CostFunction(0,alpamps);

   for(i=0;i<size_alpamp;i++){

      minus_bn = 0.0;
      plus_bn = 0.0;
      // Move alpha positive and redo x,y,z coil calc
      
      *(alpamps+i) += h;
      plus_bn = CostFunction(0,alpamps);

      *(alpamps+i) -= 2*h;
      minus_bn = CostFunction(0,alpamps);  

      *(alpamps+i) += h;
      *(derivs+i) = (plus_bn-minus_bn)/(2*h);
   }
}


void Steepest_descent( void ){
   
   int iCoils = Ncoils / Nfp;
   //int iCoils = Ncoils;
   int size_alpamp = iCoils*(2*NFalpha+1);  
   descent_dir = (double*) malloc( size_alpamp*sizeof(double) );
   int i,j;

   for(i=0;i<size_alpamp;i++){
      *(descent_dir+i) = -1.0 * (*(derivs+i));
      printf("The descent of amp %d is %.12f \n",i,*(descent_dir+i));
   }   

}

void Forward_track( void ){

   //int iCoils = Ncoils;
   int iCoils = Ncoils / Nfp;
   int size_alpamp = iCoils*(2*NFalpha+1);   
   int i,j;
   
   double step = .00000001; // There is small error, I fix later
   double init_bn = 0.0;
   double search_bn;
   double hold_bn;

   init_bn = CostFunction(0,alpamps);
   hold_bn = init_bn;
   search_bn = 0.0;
   
   int k=0;
   while( hold_bn - search_bn > 0.0 ){
      if(k==0){
         search_bn = init_bn;
      }
      hold_bn = search_bn;
      search_bn = 0.0;
      printf("Step size is %f\n", 1000.0*step);
      
      for(j=0;j<size_alpamp;j++){
	 printf("NF:   %dAlphas   %.12f\n",j,*(alpamps+j));
         *(alpamps+j) += step*(*(descent_dir+j));
      
      }
      search_bn = CostFunction(0,alpamps);
      printf("Total field error, tracking iter: %.9f   %d\n",search_bn,k);
      step = step*2.0;
      k++;
   }
   
   for(j=0;j<size_alpamp;j++){
      *(alpamps+j) -= step*(*(descent_dir+j))/2.0;
   }
   CalculateMultiFilaments();
   MultiFilFieldSym();

}


