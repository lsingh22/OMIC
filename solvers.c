//This is for initializing user inputs

#include "solvers.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "multi_fil.h"
#include "alpha.h"
// GLOBALS SCOPED IN SOURCE FILE

double* alpampsinit;
double* derivs;
double* alpamps;
double* derivs;
double* descent_dir;

int case_alpha;
int Ncoils;
int Nfp;
int NFalpha;
double alp_const;

//TODO: In the future, will need to change indexing if want NFalpha to differ for each coil

void Central_diff( void ){

 //  int iCoils = Ncoils / Nfp;
   int iCoils = Ncoils;
   int size_alpamp = iCoils*(2*NFalpha+1); // Maybe make this a global  
   int i,j;

   derivs = (double*) malloc( size_alpamp*sizeof(double) );
   double h = .00000001;
   
   double minus_bn;
   double plus_bn;
   double init_bn = 0.0;

   MultiFilField(); // Calculates B dot n for the multifilaments might not be needed here 
   for(i=0;i<Nzeta*Nteta;i++){
      init_bn += sqrt( (*(Bmfiln+i)) * (*(Bmfiln+i)) ); // Maybe make this global
   }
   for(i=0;i<size_alpamp;i++){

      minus_bn = 0.0;
      plus_bn = 0.0;
      // Move alpha positive and redo x,y,z coil calc
      *(alpamps+i) += h;
      Unpack_alpha( 0 );
      CalculateBuildDirections();
      CalculateMultiFilaments();
      MultiFilField();
      for(j=0;j<Nzeta*Nteta;j++){
         plus_bn += sqrt( (*(Bmfiln+j)) * (*(Bmfiln+j)) );
      }

      *(alpamps+i) -= 2*h;
      Unpack_alpha( 0 );
      CalculateBuildDirections();
      CalculateMultiFilaments();
      MultiFilField();
      for(j=0;j<Nzeta*Nteta;j++){
         minus_bn += sqrt( (*(Bmfiln+j)) * (*(Bmfiln+j)) );
      }
      
      *(alpamps+i) += h;

      *(derivs+i) = (plus_bn-minus_bn)/(2*h);
   }
}


void Steepest_descent( void ){
   
   //int iCoils = Ncoils / Nfp;
   int iCoils = Ncoils;
   int size_alpamp = iCoils*(2*NFalpha+1);  
   descent_dir = (double*) malloc( size_alpamp*sizeof(double) );
   int i,j;

   for(i=0;i<size_alpamp;i++){
      *(descent_dir+i) = -1.0*(*(derivs+i));
   }
   
}

void Forward_track( void ){

   int iCoils = Ncoils;
   int size_alpamp = iCoils*(2*NFalpha+1);   
   int i,j;

   double step = .0000001; // There is small error, I fix later
   double init_bn = 0.0;
   double search_bn;
   double hold_bn;

   for(i=0;i<Nzeta*Nteta;i++){
      init_bn += sqrt( (*(Bmfiln+i)) * (*(Bmfiln+i)) );
   }
   hold_bn = init_bn;
   search_bn = 0.0;
   int k = 0;
   while( hold_bn - search_bn > 0.0 ){
      if(k==0){
         search_bn = init_bn;
      }
      hold_bn = search_bn;
      search_bn = 0.0;
      printf("Step size is %f\n", 1000.0*step);
      for(j=0;j<size_alpamp;j++){
	 printf("NF:   %dAlphas   %f\n",j,*(alpamps+j));
         *(alpamps+j) += step*(*(descent_dir+j));
      }

      Unpack_alpha( 0 );
      CalculateBuildDirections();
      CalculateMultiFilaments();
      MultiFilField();
      for(j=0;j<Nzeta*Nteta;j++){
         search_bn += sqrt( (*(Bmfiln+j)) * (*(Bmfiln+j)) );
      }
      printf("Total field error, tracking iter: %f   %d\n",search_bn,k);
      step = step*2.0;
      k++;
   }
   
   for(j=0;j<size_alpamp;j++){
      *(alpamps+j) -= step*(*(descent_dir+j))/2.0;
   }
   Unpack_alpha( 0 );  //Dont think this is neccessary
   CalculateBuildDirections();
   CalculateMultiFilaments();
   MultiFilField();

}


