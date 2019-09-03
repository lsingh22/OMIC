//This is for initializing user inputs

#include "alpha.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
// GLOBALS SCOPED IN SOURCE FILE

#define ERRCODE 2  
#define ERR(e) {printf("Error: %s\n",nc_strerror(e)); exit(ERRCODE);}


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

void Init_alpha( int option ){

   int iCoils = Ncoils/ Nfp;
   int size_alpamp = iCoils*(2*NFalpha+1);  
   int i;

   alpamps = (double*) malloc( size_alpamp*sizeof(double) );
   
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
	*(alpamps+i) = alp_const;
        printf("%.8f\n", *(alpamps+i));

      }
   }
   else if(option == 2)
   {
      int ncid, varid, dimid,retval;
      nc_open(multi_output, NC_NOWRITE, &ncid);
      nc_inq_varid(ncid, "alpha", &varid);
      nc_get_var_double(ncid, varid, alpamps);
  
      for(i=0;i<size_alpamp;i++)
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

void Unpack_alpha( void ){
   
   int iCoils = Ncoils / Nfp;
   //int iCoils = Ncoils;
   int size_alpamp = iCoils*(2*NFalpha+1);  
   int size_alp = iCoils*(Nseg+1);
   int i,j,k;
   double theta, a;  
   double pi = M_PI;
    
   alp = (double*) malloc(iCoils*(Nseg+1)*sizeof(double));

   for(i=0;i<iCoils;i++){
      for(j=0;j<(Nseg+1);j++){
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



