#include "read_namelist.h"
#include "alpha.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
 
// GLOBALS SCOPED IN SOURCE FILE
double len_rad; double len_tor; double Nturns;

char* focus_output; char* multi_output; char* sfil_output; char* mfil_output; char* omic_input;
double* alp; int Nradfil; int Ntorfil; int Nseg; int isVaryRotation; int case_alpha; int case_objfun;

int NFalpha; double alp_const; int niter; double surface_area; int nproc; int* startind; int* endind;
int isStellSym; int iCoil; int size_alpamp; int size_fp; int Nfils; int Ns;

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
 
int CFileExists(char *filename){
   // try to open and read a file
   FILE *file;
   if((file = fopen(filename, "r")))
   {
      fclose(file);
      return 1;
   }
   return 0;
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

void OMICStartup(char *ext){
/* Handle errors on startup */
   char *ext_omic_in = malloc(strlen(ext) + 1);
   strcpy(ext_omic_in, ext);
   char *ext_omic_out = malloc(strlen(ext) + 1);
   strcpy(ext_omic_out, ext);
   char *ext_mfil_out = malloc(strlen(ext) + 1);
   strcpy(ext_mfil_out, ext);
   char *ext_sfil_out= malloc(strlen(ext) + 1);
   strcpy(ext_sfil_out, ext);

   char *focus_temp    = "focus.output";
   char *omic_temp_in  = strcat(ext_omic_in,  ".input" );
   char *omic_temp_out = strcat(ext_omic_out, ".nc"    );
   char *mfil_temp_out = strcat(ext_mfil_out, ".mcoils");
   char *sfil_temp_out = strcat(ext_sfil_out, ".scoils");

   int exist_focus, exist_omic_in;  
   exist_focus    = CFileExists(focus_temp);
   exist_omic_in  = CFileExists(omic_temp_in);

   if(!exist_focus)
   {
      printf("\n1 No file '%s' found. Exiting...\n\n ",focus_temp); 
      exit(0);
   }else
   {
      //printf("\nFOCUS output file '%s' loaded.",focus_temp);
      if(!exist_omic_in)
      {
         printf("\n2 No file '%s' found. Exiting...\n\n ",omic_temp_in); 
         exit(0);
      }else
      {
      //   printf("\nOMIC input file '%s' found.\n",omic_temp_in);    
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

