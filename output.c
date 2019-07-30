/*

#include <stdio.h>
#include "globals.h"


#define SFILB_FILE_NAME "./outputfiles/sfilB.nc"
   
void WriteSingleB(void){
   //Write to NC
   int ncid, xvarid, yvarid, zvarid, bvarid, xdimid, ydimid;
   int bxvarid, byvarid, bzvarid;
   int dimids[2];
   
   nc_create(SFILB_FILE_NAME, NC_CLOBBER, &ncid); 
   nc_def_dim(ncid, "Nzeta", Nzeta, &xdimid);
   nc_def_dim(ncid, "Nteta", Nteta, &ydimid);
   dimids[0] = xdimid;
   dimids[1] = ydimid;
   
   nc_def_var(ncid, "xsurf", NC_DOUBLE, 2, dimids, &xvarid);
   nc_def_var(ncid, "ysurf", NC_DOUBLE, 2, dimids, &yvarid);
   nc_def_var(ncid, "zsurf", NC_DOUBLE, 2, dimids, &zvarid);  
   nc_def_var(ncid, "B", NC_DOUBLE, 2, dimids, &bvarid);  
   nc_def_var(ncid, "Bx", NC_DOUBLE, 2, dimids, &bxvarid);
   nc_def_var(ncid, "By", NC_DOUBLE, 2, dimids, &byvarid);
   nc_def_var(ncid, "Bz", NC_DOUBLE, 2, dimids, &bzvarid);  
 

   nc_enddef(ncid);
   nc_put_var_double(ncid, xvarid, &xsurf[0]);
   nc_put_var_double(ncid, yvarid, &ysurf[0]);
   nc_put_var_double(ncid, zvarid, &zsurf[0]);
   nc_put_var_double(ncid, bvarid, &Bsfil[0]);
   nc_put_var_double(ncid, bxvarid, &Bsfilx[0]);
   nc_put_var_double(ncid, byvarid, &Bsfily[0]);
   nc_put_var_double(ncid, bzvarid, &Bsfilz[0]);
 
   nc_close(ncid);
}


void WriteSingleFilaments(void){
  
   int i,j;
   FILE* fb;
   fb = fopen("./outputfiles/sfil.out","w");
   fprintf(fb, "periods 1\n begin filament\n mirror NIL\n");

   for(i=0;i<Ncoils;i++){
      for(j=0;j<Nseg;j++){
         fprintf(fb,"%.15f %.15f %.15f %.15f \n", *(sfilx+i*(Nseg+1)+j), *(sfily+i*(Nseg+1)+j), *(sfilz+i*(Nseg+1)+j), *(currents+i));         
         }
      fprintf(fb,"%.15f %.15f %.15f %.15f Mod %d\n", *(sfilx+i*(Nseg+1)), *(sfily+i*(Nseg+1)), *(sfilz+i*(Nseg+1)), *(currents+i), i+1);         
   }
   fprintf(fb,"end");













#define MFILB_FILE_NAME "./outputfiles/mfilB.nc"
   
void WriteMultiB(void){
   //Write to NC
   int ncid, xvarid, yvarid, zvarid, bvarid, xdimid, ydimid;
   int bxvarid, byvarid, bzvarid;
   int dimids[2];
   
   nc_create(MFILB_FILE_NAME, NC_CLOBBER, &ncid); 
   nc_def_dim(ncid, "Nzeta", Nzeta, &xdimid);
   nc_def_dim(ncid, "Nteta", Nteta, &ydimid);
   dimids[0] = xdimid;
   dimids[1] = ydimid;
   
   nc_def_var(ncid, "xsurf", NC_DOUBLE, 2, dimids, &xvarid);
   nc_def_var(ncid, "ysurf", NC_DOUBLE, 2, dimids, &yvarid);
   nc_def_var(ncid, "zsurf", NC_DOUBLE, 2, dimids, &zvarid);  
   nc_def_var(ncid, "B", NC_DOUBLE, 2, dimids, &bvarid);  
   nc_def_var(ncid, "Bx", NC_DOUBLE, 2, dimids, &bxvarid);
   nc_def_var(ncid, "By", NC_DOUBLE, 2, dimids, &byvarid);
   nc_def_var(ncid, "Bz", NC_DOUBLE, 2, dimids, &bzvarid);  
 

   nc_enddef(ncid);
   nc_put_var_double(ncid, xvarid, &xsurf[0]);
   nc_put_var_double(ncid, yvarid, &ysurf[0]);
   nc_put_var_double(ncid, zvarid, &zsurf[0]);
   nc_put_var_double(ncid, bvarid, &Bmfil[0]);
   nc_put_var_double(ncid, bxvarid, &Bmfilx[0]);
   nc_put_var_double(ncid, byvarid, &Bmfily[0]);
   nc_put_var_double(ncid, bzvarid, &Bmfilz[0]);

   nc_close(ncid);
}

 
void WriteMultiFilaments(void){

   int i,j,k;
   FILE* fb;
   fb = fopen("./outputfiles/mfil.out","w");
   fprintf(fb, "periods 1\n begin filament\n mirror NIL\n");
   int Nfils = Ntorfil*Nradfil;
   
   for(i=0;i<Ncoils;i++){
      for(j=0;j<Nfils;j++){
         for(k=0;k<Nseg;k++){
         fprintf(fb,"%.15f %.15f %.15f %.8f \n", *(mfilx+i*(Nseg+1)*Nfils+j*(Nseg+1)+k), \
                                                 *(mfily+i*(Nseg+1)*Nfils+j*(Nseg+1)+k), \

*/
