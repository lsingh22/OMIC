#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
#include <netcdf.h>
#include <stdlib.h>

int main(int argc, char **argv) {

//	char* file;
//	file = "/home/luquants/multi/dev/test.nc"; 
/*	
	if (argc == 1){
		file = "/home/luquants/multi/dev/test.nc";	
		}
	else{
		printf("Wrong number of inputs! \n");		
		}
*/
   FILE* fp;
   char* line = NULL;
   size_t len = 0;
   int count = 0;
   char* data;
   int* n;
   int* m;
   int i,j,k,N;
   
   int ncid, varid;
   double* bsupuu;
   double* bsupvv;
   int* Ncoils;
   int* NFcoil;
   //bsupuu = (double *) malloc(59*101*sizeof(double));
   //bsupvv = (double *) malloc(59*101*sizeof(double));
   Ncoils = (int *) malloc(sizeof(int));
   NFcoil = (int *) malloc(sizeof(int));   

   nc_open("dev/test.nc", NC_NOWRITE, &ncid);
   //nc_inq_varid(ncid, "bsupumnc", &varid);
   //nc_get_var_double(ncid, varid, bsupuu);
   //nc_inq_varid(ncid, "bsupvmnc", &varid);
   //nc_get_var_double(ncid, varid, bsupvv);
   nc_inq_varid(ncid, "Ncoils", &varid);
   nc_get_var_int(ncid, varid, Ncoils);
   nc_inq_varid(ncid, "NFcoil", &varid);
   nc_get_var_int(ncid, varid, NFcoil);
   nc_close(ncid);
   
   printf("%d\n", *NFcoil);
   printf("%d\n", *Ncoils);

   for( i = 0; i < 59; i++ ){
      //printf("%f   %d\n",bsupuu[i],i);
   }





   
}
















