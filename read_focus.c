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
   //nc_inq_vardimid(ncid,varid,&dimid);
   //nc_inq_dimlen(ncid,dimid,&size_surf);
   nc_close(ncid);
}

// ALLOCATES AND STORES FOCUS DATA
void ReadFocusArrays(char* output_file){
   
   int ncid, varid, dimid;
   coilspace = (double*) malloc(size_coilspace*sizeof(double));   
   xsurf = (double*) malloc(Nteta*Nzeta*sizeof(double)); 
   ysurf = (double*) malloc(Nteta*Nzeta*sizeof(double)); 
   zsurf = (double*) malloc(Nteta*Nzeta*sizeof(double)); 
   nsurfx = (double*) malloc(Nteta*Nzeta*sizeof(double)); 
   nsurfy = (double*) malloc(Nteta*Nzeta*sizeof(double)); 
   nsurfz = (double*) malloc(Nteta*Nzeta*sizeof(double)); 
   fbx = (double*) malloc(Nteta*Nzeta*sizeof(double)); 
   fby = (double*) malloc(Nteta*Nzeta*sizeof(double)); 
   fbz = (double*) malloc(Nteta*Nzeta*sizeof(double)); 
   fbn = (double*) malloc(Nteta*Nzeta*sizeof(double)); 
     
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

   for(i=0;i<Nzeta*Nteta;i++){
      fprintf(fb,"%.15f %.15f %.15f \n", *(xsurf+i),*(ysurf+i) ,*(zsurf+i));         
   }
}

#define FILE_NAME "./outputfiles/boundary.nc"
   
void WriteBoundaryNC(void){

   int ncid, xvarid, yvarid, zvarid, xdimid, ydimid;
   int dimids[2];
   
   nc_create(FILE_NAME, NC_CLOBBER, &ncid); 
   nc_def_dim(ncid, "size_surf1", Nzeta, &xdimid);
   nc_def_dim(ncid, "size_surf2", Nteta, &ydimid);
   dimids[0] = xdimid;
   dimids[1] = ydimid;
   nc_def_var(ncid, "xsurf", NC_DOUBLE, 2, dimids, &xvarid);
   nc_def_var(ncid, "ysurf", NC_DOUBLE, 2, dimids, &yvarid);
   nc_def_var(ncid, "zsurf", NC_DOUBLE, 2, dimids, &zvarid);  

   nc_enddef(ncid);
   nc_put_var_double(ncid, xvarid, &xsurf[0]);
   nc_put_var_double(ncid, yvarid, &ysurf[0]);
   nc_put_var_double(ncid, zvarid, &zsurf[0]);
   nc_close(ncid);
}


