//This file reads the FOCUS output 

#include "read_focus.h"
#include <stdlib.h>
#include <stdio.h>
// GLOBALS SCOPED IN SOURCE FILE

int Nseg;
int Ncoils;
int Nfp;
int isSym;
int NFcoil; //TODO: Update when bug is fixed 
int Nteta;
int Nzeta;
size_t size_coilspace;
size_t size_surf;
double* coilspace;
double* xsurf;
double* ysurf;
double* zsurf;
double* nsurfx;
double* nsurfy;
double* nsurfz;
double* fbn;
double* fbx;
double* fby;
double* fbz;

// READS AND STORES VARIOUS FOCUS COILSET PARAMETERS
// DETERMINES SIZE TO ALLOCATE OUTPUT FOCUS DATA
void ReadFocusInts(char* output_file){
   
   int ncid, varid, dimid;  
   nc_open(output_file, NC_NOWRITE, &ncid);
   nc_inq_varid(ncid, "Ncoils", &varid);
   nc_get_var_int(ncid, varid, &Ncoils);
   nc_inq_varid(ncid, "Nfp", &varid);
   nc_get_var_int(ncid, varid, &Nfp);
   nc_inq_varid(ncid, "IsSymmetric", &varid);
   nc_get_var_int(ncid, varid, &isSym);
   nc_inq_varid(ncid, "Nseg", &varid);
   nc_get_var_int(ncid, varid, &Nseg);
   nc_inq_varid(ncid, "NFcoil", &varid);
   nc_get_var_int(ncid, varid, &NFcoil);
   nc_inq_varid(ncid, "Nteta", &varid);
   nc_get_var_int(ncid, varid, &Nteta);
   nc_inq_varid(ncid, "Nzeta", &varid);
   nc_get_var_int(ncid, varid, &Nzeta);
   nc_inq_varid(ncid,"coilspace",&varid);
   nc_inq_vardimid(ncid,varid,&dimid);
   nc_inq_dimlen(ncid,dimid,&size_coilspace);  
   nc_inq_varid(ncid,"xsurf",&varid);
   nc_inq_vardimid(ncid,varid,&dimid);
   nc_inq_dimlen(ncid,dimid,&size_surf);
   nc_close(ncid);
}

// ALLOCATES AND STORES FOCUS DATA
void ReadFocusArrays(char* output_file){
   
   int ncid, varid, dimid;
   coilspace = (double*) malloc(size_coilspace*sizeof(double));   
   xsurf = (double*) malloc(size_surf*size_surf*sizeof(double)); 
   ysurf = (double*) malloc(size_surf*size_surf*sizeof(double)); 
   zsurf = (double*) malloc(size_surf*size_surf*sizeof(double)); 
   nsurfx = (double*) malloc(size_surf*size_surf*sizeof(double)); 
   nsurfy = (double*) malloc(size_surf*size_surf*sizeof(double)); 
   nsurfz = (double*) malloc(size_surf*size_surf*sizeof(double)); 
   fbx = (double*) malloc(size_surf*size_surf*sizeof(double)); 
   fby = (double*) malloc(size_surf*size_surf*sizeof(double)); 
   fbz = (double*) malloc(size_surf*size_surf*sizeof(double)); 
   fbn = (double*) malloc(size_surf*size_surf*sizeof(double)); 
     
   nc_open(output_file, NC_NOWRITE, &ncid);
   nc_inq_varid(ncid, "coilspace", &varid);
   nc_get_var_double(ncid, varid, coilspace);
   
   nc_inq_varid(ncid, "xsurf", &varid);
   nc_get_var_double(ncid, varid, xsurf);
   nc_inq_varid(ncid, "ysurf", &varid);
   nc_get_var_double(ncid, varid, ysurf);
   nc_inq_varid(ncid, "zsurf", &varid);
   nc_get_var_double(ncid, varid, zsurf);

   nc_inq_varid(ncid, "nx", &varid);
   nc_get_var_double(ncid, varid, nsurfx);
   nc_inq_varid(ncid, "ny", &varid);
   nc_get_var_double(ncid, varid, nsurfy);
   nc_inq_varid(ncid, "nz", &varid);
   nc_get_var_double(ncid, varid, nsurfz);

   nc_inq_varid(ncid, "Bn", &varid);
   nc_get_var_double(ncid, varid, fbn);
   nc_inq_varid(ncid, "Bx", &varid);
   nc_get_var_double(ncid, varid, fbx);
   nc_inq_varid(ncid, "By", &varid);
   nc_get_var_double(ncid, varid, fby);
   nc_inq_varid(ncid, "Bz", &varid);
   nc_get_var_double(ncid, varid, fbz);
   nc_close(ncid);
}

void WriteBoundary(void){
  
   int i,j;
   FILE* fb;
   fb = fopen("./outputfiles/boundary.out","w");

   for(i=0;i<size_surf*size_surf;i++){
      fprintf(fb,"%.15f %.15f %.15f \n", *(xsurf+i),*(ysurf+i) ,*(zsurf+i));         
   }
}




