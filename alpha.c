//This is for initializing user inputs

#include "alpha.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
// GLOBALS SCOPED IN SOURCE FILE

char* alpha_input;
double* alpampsinit;
double* alpamps;
int case_alpha;
int Ncoils;
int Nfp;
int NFalpha;
double alp_const;

//TODO: In the future, will need to change indexing if want NFalpha to differ for each coil

void Init_alpha( int option){

 //  int iCoils = Ncoils / Nfp;
   int iCoils = Ncoils;
   int size_alpamp = iCoils*(2*NFalpha+1);  
   int i;

   alpampsinit = (double*) malloc( size_alpamp*sizeof(double) );
   alpamps = (double*) malloc( size_alpamp*sizeof(double) );
   
   if(option == 0)
   {
      for(i=0;i<size_alpamp;i++)
      {   
         *(alpampsinit+i) = 0.0;
      }
   }
   else if(option == 1)
   {
      for(i=0;i<size_alpamp;i++)
      {   
        *(alpampsinit+i) = alp_const;
        //printf("%.8f\n", *(alpampsinit+i));
	*(alpamps+i) = 0.0; //alp_const;
      }
   }
   else if(option == 2)
   {
      //TODO: Read namelist or separate file
   }
   else
   {
      printf("Alphas not stored... \n Exiting!!");
      //TODO: Throw error 
   }
}


void Unpack_alpha( int isInit){
   
   //int iCoils = Ncoils / Nfp;
   int iCoils = Ncoils;
   int size_alpamp = iCoils*(2*NFalpha+1);  
   int size_alp = iCoils*(Nseg+1);
   int i,j,k;
   double theta, a;  
   double pi = M_PI;
    
   alp = (double*) malloc(iCoils*(Nseg+1)*sizeof(double));
   
   if(isInit == 0) // Optimization branch
   {
      for(i=0;i<iCoils;i++){
         for(j=0;j<Nseg+1;j++){
            theta = ((2*pi)/Nseg)*j;
            a = 0;
            for(k=0;k<NFalpha+1;k++){
               a = a + alpamps[ (2*NFalpha+1)*i + k ]*cos(k*theta);
            }
            for(k=1;k<NFalpha+1;k++){
               a = a + alpamps[ (2*NFalpha+1)*i + NFalpha + k ]*sin(k*theta);
            }

            *(alp +i*(Nseg+1) + j ) = a;
         }
      }     
   }
   else if(isInit == 1) // Initialization branch
   {
 
      for(i=0;i<iCoils;i++){
         for(j=0;j<Nseg+1;j++){
            theta = ((2*pi)/Nseg)*j;
            a = 0;
            for(k=0;k<NFalpha+1;k++){
               a = a + alpampsinit[ (2*NFalpha+1)*i + k ]*cos(k*theta);
            }
            for(k=1;k<NFalpha+1;k++){
               a = a + alpampsinit[ (2*NFalpha+1)*i + NFalpha + k ]*sin(k*theta);
            }

            *(alp +i*(Nseg+1) + j ) = a;
         }
      } 
   }
   else
   {
      //Throw error
   }

}


