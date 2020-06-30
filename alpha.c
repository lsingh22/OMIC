//This is for initializing user inputs

#include "alpha.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
// GLOBALS SCOPED IN SOURCE FILE

#define ERRCODE 2  
#define ERR(e) {printf("Error: %s\n",nc_strerror(e)); exit(ERRCODE);}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

char* multi_output;
double* alpampsinit;
double* alpamps;
int case_alpha;
int Ncoils;
int Nfp;
int NFalpha;
double alp_const;
double* alp;

//TODO: In the future, will need to change indexing if want NFalpha to differ for each coil

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

void Init_alpha( int option ){
//----------------------------------------------------------------------------------------------------
// Initializes the alpha optimization function Fourier amplitudes 
// Depending on the value of 'option', the following are possible:
//  0 - initialize all amplitudes to zero 
//  1 - set the amplitudes to a constant value 
//  2 - set the amplitudes to those given in a previous .nc output file
// Be careful with option 2, as the output is overridden (best to copy output to another file name)
//----------------------------------------------------------------------------------------------------

   int iCoils = Ncoils/ Nfp; //Nfp here accounts for periodicity
   int size_alpamp = iCoils*(2*NFalpha+1); //each coil has N_alpha+1 and N_alpha non-zero cos, sin terms in its representation  
   int i;

   alpamps = (double*) malloc( size_alpamp*sizeof(double) ); //stores all Fourier components of alpha angle for all coils
   
      if(option == 0)
   {
      for(i=0;i<size_alpamp;i++)
      {   
         *(alpamps+i) = 0.0;
      }
   }
   else if(option == 1)
   {
      for(i=0;i<size_alpamp;i++)
      {   
	*(alpamps+i) = alp_const; //alp_const is specified in the input by user
        printf("%.8f\n", *(alpamps+i));

      }
   }
   else if(option == 2)
   {
      int size_amp = iCoils*(2*NFalpha+1)/Nfp;
      int ncid, varid, dimid,retval;
      nc_open(multi_output, NC_NOWRITE, &ncid); //multi_output is the path to the .nc output file
      nc_inq_varid(ncid, "alpha", &varid);
      nc_get_var_double(ncid, varid, alpamps);
  
      for(i=0;i<size_amp;i++)
      {   
        printf("%.8f\n", *(alpamps+i) );
      } 

      if(retval=nc_inq_varid(ncid,"alpha",&varid))
      ERR(retval);
      
   }
   else
   {
      printf("Alphas not stored... \n Exiting!!");
      //TODO: Throw error 
   }
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

void Unpack_alpha( void ){
//----------------------------------------------------------------------------------------------------
// Calculate the coil rotation angle alpha based on its Fourier series
//----------------------------------------------------------------------------------------------------
  
   int iCoils = Ncoils / Nfp;
   int size_alpamp = iCoils*(2*NFalpha+1); //total number of Fourier amplitudes for all coils  
   int size_alp = iCoils*(Nseg+1); //total number of points for all coils
   int i,j,k;
   double theta, a;  
   double pi = M_PI;
    
   alp = (double*) malloc(iCoils*(Nseg+1)*sizeof(double)); //all alpha angles at each position on each coil

   for(i=0;i<iCoils;i++){
      for(j=0;j<(Nseg+1);j++){
         theta = ((2*pi)/Nseg)*j; //parameterizing variable
         a = 0; //The sum starts with a=0...
         for(k=0;k<NFalpha+1;k++){
            a = a + alpamps[ (2*NFalpha+1)*i + k ]*cos(k*theta); //then adds the cosine terms...
         }
         for(k=1;k<NFalpha+1;k++){
            a = a + alpamps[ (2*NFalpha+1)*i + NFalpha + k ]*sin(k*theta); //and then the sine terms...
         }
            *(alp +i*(Nseg+1) + j ) = a; //before storing the angle in 'alp' array.
      }
   }

}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
