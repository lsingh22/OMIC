#include "alpha.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define ERRCODE 2  
#define ERR(e) {printf("Error: %s\n",nc_strerror(e)); exit(ERRCODE);}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

// GLOBALS SCOPED IN SOURCE FILE
char* multi_output; double* alpampsinit; double* alpamps; int case_alpha; int iCoil;
int Nfp; int NFalpha; double alp_const; double* alp; int size_alpamp;

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

void Init_alpha( int option ){
//----------------------------------------------------------------------------------------------------
// Initializes the alpha optimization function Fourier amplitudes 
// Depending on the value of 'option', the following are possible:
//  0 - initialize all amplitudes to zero 
//  1 - set the amplitudes to a constant value 
//  2 - set the amplitudes to those given in a previous .nc output file
// Be careful with option 2, as the output is overridden (best to copy output to another file first)
//----------------------------------------------------------------------------------------------------

   register int i;

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
      }
   }
   else if(option == 2)
   {
      int size_amp = iCoil * (2*NFalpha+1); //TODO: should probably override if NFalpha doesn't match previous optimization!
      int ncid, varid, dimid,retval;
      nc_open(multi_output, NC_NOWRITE, &ncid); //multi_output is the path to the .nc output file
      nc_inq_varid(ncid, "alpha", &varid);
      nc_get_var_double(ncid, varid, alpamps);
  
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
  
   int size_alp = iCoil*(Nseg+1); //total number of points for all coils
   register int i,j;
   int k;
   double theta, a;  
   double pi = M_PI;
    
   alp = (double*) malloc(iCoil*(Nseg+1)*sizeof(double)); //all alpha angles at each position on each coil

   for(i=0;i<iCoil;i++)
   {
      for(j=0;j<(Nseg+1);j++)
      {
         theta = ((2*pi)/Nseg)*j; //parameterizing variable
         a = 0; 
         for(k=0;k<NFalpha+1;k++)
         {
            a = a + alpamps[(2*NFalpha+1)*i+k] * cos(k*theta); 
         }
         for(k=1;k<NFalpha+1;k++)
         {
            a = a + alpamps[(2*NFalpha+1)*i+NFalpha+k] * sin(k*theta);
         }
         *(alp + i*(Nseg+1) + j ) = a; 
      }
   }
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
