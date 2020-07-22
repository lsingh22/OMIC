
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include "read_namelist.h"
#include "alpha.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "libfyaml.h"

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
 
// GLOBALS SCOPED IN SOURCE FILE
double len_rad;
double len_tor;
double Nturns;

char* focus_output;
char* multi_output;
char* sfil_output;
char* mfil_output;
char* omic_input;

double* alp;
int Nradfil;
int Ntorfil;
int Nseg;
int isVaryRotation;
int case_alpha;
int case_objfun;
int NFalpha;
double alp_const;
int niter;
double surface_area;
int nproc;
int* startind;
int* endind;
int isStellSym;
int case_optimize;
int iCoil;
int size_alpamp;
int size_fp;
int Nfils; 
int Ns;

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
 
int CFileExists(char *filename){
   /* try to open and read a file */
   FILE *file;
   if(file = fopen(filename, "r"))
   {
      fclose(file);
      return 1;
   }
   return 0;
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

void Startup(char *ext){
/* Handle errors on startup */
   char *ext_omic_in = malloc(strlen(ext) + 1);
   strcpy(ext_omic_in, ext);
   char *ext_omic_out = malloc(strlen(ext) + 1);
   strcpy(ext_omic_out, ext);
   char *ext_mfil_out = malloc(strlen(ext) + 1);
   strcpy(ext_mfil_out, ext);
   char *ext_sfil_out= malloc(strlen(ext) + 1);
   strcpy(ext_sfil_out, ext);

//   char *prefix = "omic_";
//   char *suffix = ".nc";
   
//   char *omic_temp_out = malloc(strlen(prefix) + strlen(suffix) + strlen(ext_omic_out) + 1);
//   omic_temp_out = strcat(strcat(prefix,ext_omic_out),suffix);// = malloc(strlen(prefix) + strlen(suffix) + strlen(ext_omic_out) + 1);

//   strcpy(omic_temp_out,prefix);
//   strcpy(omic_temp_out,ext_omic_out);
//   strcpy(omic_temp_out,suffix);
//   printf("The string I got is %s \n", omic_temp_out); 
   char *focus_temp    = "focus.output";
   char *omic_temp_in  = strcat(ext_omic_in,  ".input" );
   char *omic_temp_out = strcat(ext_omic_out, ".nc"    );
   char *mfil_temp_out = strcat(ext_mfil_out, ".mcoils");
   char *sfil_temp_out = strcat(ext_sfil_out, ".scoils");

   int exist_focus, exist_omic_in, exist_omic_out; //TODO: exist_boundary 
   exist_focus    = CFileExists(focus_temp);
   exist_omic_in  = CFileExists(omic_temp_in);

   //TODO: in case where we read in previous input   exist_omic_out = CFileExists(omic_temp_in);
   //Would require fy_document scan prior to reading most inputs, which is probably fine

   if(!exist_focus)
   {
      printf("\n1 No file '%s' found. Exiting...\n\n ",focus_temp); 
      exit(0);
   }else
   {
      printf("\nFOCUS output file '%s' loaded.",focus_temp);
      if(!exist_omic_in)
      {
         printf("\n2 No file '%s' found. Exiting...\n\n ",omic_temp_in); 
         exit(0);
      }else
      {
         printf("\nOMIC input file '%s' found.\n",omic_temp_in);    
      }     
   }
  
   focus_output = focus_temp; 
   omic_input   = omic_temp_in; //TODO: not currently supported
   multi_output = omic_temp_out;
   mfil_output  = mfil_temp_out;
   sfil_output  = sfil_temp_out;
   
//   free(ext_omic_in);
//   free(ext_omic_out);
/*   free(ext_mfil_out);
   free(ext_sfil_out); */
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

void ReadInputFile(void){
   /* Read YAML input file */
   int count;
   struct fy_document *fyd;
//   fyd = fy_document_build_from_file(NULL, multi_output);

//   if(!fyd)
//   {
//      printf("ERROR: Failed to build input file.\n");
//      exit(0);
//   }  
     
//   Try this before writing everything incorrectly...
     
//   count = fy_document_scanf(fyd,"/case_optimize %d ", &case_optimize); 
//   if(count!=1){printf("ERROR: Please specify case_optimize input parameter.\n");


//   Rest of the parameters to assign to globals
/* 
   niter = fy_document_scanf(fyd, "/niter %d ");
   case_alpha = 0;
   NFalpha = 1;
   alp_const = 0.000;
 
   isStellSym = 0;
   Nseg = 128;

   weight_comp = 0.01; //complexity weighting
   case_objfun = 1; //0 for fbn , 1 for both fbn and fc
   nvals_scaling = 2; // the beta in the complexity formulation   

   Nradfil = 2;
   Ntorfil = 2;

   len_rad = 0.120; // 1/2 hsx 
   len_tor = 0.060;
*/ 

//   count = EXIT_SUCCESS;

}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

/*
void Initialize(void){

   if(isStellSym == 1)
   {
      iCoil = Ncoil / ( 2 * Nfp );
      size_fp = Nteta * Nzeta / ( 2 * Nfp ); 
      Ns = 1; 
   }
   else if(isStellSym == 0)
   {
      iCoil = Ncoil / Nfp;
      size_fp = Nteta * Nzeta / Nfp; 
      Ns = 0;
   }

   size_alpamp = iCoil * ( 2 * NFalpha + 1 ); 
   Nfils = Nradfil * Ntorfil; 
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

void MPInit(void){
//----------------------------------------------------------------------------------------------------
// Allocate and write some useful arrays for mpi implementation
//----------------------------------------------------------------------------------------------------
 
   startind = (int*) malloc(nproc*sizeof(int));
     endind = (int*) malloc(nproc*sizeof(int)); 

   int size_fp = Nteta*Nzeta / Nfp;
   int floor_points_per_proc = (size_fp - (size_fp % 2)) / nproc; //returns even average of points per processor
   int endpn = nproc - 1; //the index value of last processor

   for(int i=0;i<nproc;i++)
   {
      *(startind+i) = 0 + i*floor_points_per_proc;
      *(endind+i)   = (i+1)*floor_points_per_proc - 1; 
   }
   
   //increment the last index if there are an odd number of points
   if(size_fp % 2 == 1)
   {
      *(endind+endpn) = *(endind+endpn) + 1;
   }
}
*/
//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----


