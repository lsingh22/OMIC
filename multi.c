//This will serve as the main method for the multifilament code
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
   if( argc != 1 ){
      printf("Wrong number of inputs!\n");
      return 1;
   }
   if( argc != 2){
      printf("Luquant can debug in c now!\n");
      return 2;
   }
   //Allocate pointers and call function to load wout file, load_wout.c
   //
   //Allocate pointers and call function to load .input file, load_input.c
   //
   //Allocate pointers and call function to load .focus file, load_focus.c
   //
   //Call function to produce XYZ points for single filament, single_fil.c
   //
   //Write python code to plot filaments, single or multi, plotFil.py
   //
   //Hardcode alphas to 0 and call function to produce XYZ points for multifilament, multi_fil.c
   //
   //write python code to plot finite build coils, plotBuild.py 
   
   return 0;
}
