#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
#include <netcdf.h>
#include <stdlib.h>

int main(int argc, char **argv) {

//NOTE: Due to bug in FOCUS, this script is not completely general yet -- NFcoil needs to be fixed.

//Using /home/luquants/focusruns/boxport/hsx.bp_00/focus_hsx.bp_00.h5 as a test case.


// Read the input .nc file//

//      char* file;
//      file = "/home/luquants/multi/dev/test.nc"; 
/*      
        if (argc == 1){
                file = "/home/luquants/multi/dev/test.nc";      
                }
        else{
                printf("Wrong number of inputs! \n");           
                }
*/

//Define and allocate pointers// 

   FILE* fp;
   char* line = NULL;
   size_t len = 0;
   int count = 0;
   char* data;
   int* n;
   int* m;
   int i,j,N;
   int sum = 0;

   int ncid, varid;
   int* Ncoils;
   int* Nfp;
   int* isSym;
   //TODO: int* NFcoil;
   int NFcoil[11] = { 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8};

   Ncoils = (int *) malloc(sizeof(int));
//   NFcoil = (int *) malloc(sizeof(int));   
   Nfp = (int *) malloc(sizeof(int));
   isSym = (int *) malloc(sizeof(int));

   //TODO: Notify CZ about NFcoil bug

   nc_open("dev/test.nc", NC_NOWRITE, &ncid);
   nc_inq_varid(ncid, "Ncoils", &varid);
   nc_get_var_int(ncid, varid, Ncoils);
   nc_inq_varid(ncid, "Nfp", &varid);
   nc_get_var_int(ncid, varid, Nfp);
   nc_inq_varid(ncid, "IsSymmetric", &varid);
   nc_get_var_int(ncid, varid, isSym);

   //printf("%d\n", *Nfp);
   //printf("%d\n", *isSym);


   //TODO: Fix the NFcoil storage after bug is fixed
   //nc_inq_varid(ncid, "NFcoil", &varid);
   //nc_get_var_int(ncid, varid, NFcoil);

   //Determine the total number of FS amplitudes// 
   for (i =0; i < (*Ncoils);i++){
        sum = sum + NFcoil[i];
}

   //Store currents and FS harmonics for each coil//

   double* coildata;
   double* centroids;
   double* coilamps;
   double* currents;

   coildata = (double *) malloc(59*101*sizeof(double));
   centroids = (double *) malloc((*Ncoils)*3*sizeof(double));
   coilamps = (double *) malloc(sum*6*sizeof(double));
   currents = (double *) malloc((*Ncoils)*sizeof(double));
   nc_inq_varid(ncid, "coilspace", &varid);
   nc_get_var_double(ncid, varid, coildata);
   nc_close(ncid);

   int ind = 0;
   int ind_arr[*Ncoils];

   for(i=0;i<(*Ncoils);i++){
        if (i==0){
                currents[i] = coildata[ind];}
        else{
                currents[i] = coildata[ind-1];
        }

        for(j=0;j<6*(NFcoil[i]+1)-3;j++){
                        coilamps[ind + j] = coildata[ind + j + i + 1 ];
                        printf("%f   %d\n",coilamps[ind+j],ind+j);
                }
                ind_arr[i] = ind;
                //printf("%d\n",ind_arr[i]);
                ind = ind + ((NFcoil[i])*6+3);
        }

   //Store centroids from coilamps array//


   for(i=0;i<(*Ncoils);i++){
        centroids[i*3 + 1] = coilamps[ind_arr[i]];
        centroids[i*3 + 2] = coilamps[ind_arr[i] + 2*(NFcoil[i] + 1)-1];
        centroids[i*3 + 3] = coilamps[ind_arr[i] + 4*(NFcoil[i] + 1)-2];
        printf("%f  %f  %f\n", centroids[i*3 + 1], centroids[i*3 + 2], centroids[i*3 + 3]);

   }

