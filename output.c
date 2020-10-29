#include <stdio.h>
#include "globals.h"

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
 
double* xsurf; double* ysurf; double* zsurf; double* sfilx; double* sfily; double* sfilz;
 
double* mfilx; double* mfily; double* mfilz; double* ffilx; double* ffily; double* ffilz;
  
double* alpamps; /* double* Bsfilx; double* Bsfily; double* Bsfilz; double* Bsfil; double* Bsfiln;*/

double* Bmfilx; double* Bmfily; double* Bmfilz; double* Bmfil; double* Bmfiln;

double* cx; double* cy; double* cz; double* bx; double* by; double* bz;

double* tx; double* ty; double* tz; double* nx; double* ny; double* nz;

int Nseg, Ntorfil, Nradfil, Nfils, Nzeta, Nteta, Ncoil, iCoil, NFalpha, Nfp;

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

#define FILE_NAME multi_output

void WriteOutputNC(void){
//----------------------------------------------------------------------------------------------------
// Writes some relevant outputs to a netcdf file format
//----------------------------------------------------------------------------------------------------

   int ncid;
   int xvarid, yvarid, zvarid;
   int xdimid, ydimid, coildimid, symcoildimid, nsegdimid, nsegfildimid, p1dimid, ffildimid;
   int sxvarid, syvarid, szvarid;
   int cxvarid, cyvarid, czvarid;
   int txvarid, tyvarid, tzvarid;
   int nxvarid, nyvarid, nzvarid;
   int bxvarid, byvarid, bzvarid;
   int mxvarid, myvarid, mzvarid;
   int fxvarid, fyvarid, fzvarid;
   int alpvarid, alpampdimid;
//   int sbxvarid, sbyvarid, sbzvarid, sbvarid, sbnvarid;
   int mbxvarid, mbyvarid, mbzvarid, mbvarid, mbnvarid;
   int NFalphavarid, len_radvarid, len_torvarid, nitervarid, Ntorfilvarid, Nradfilvarid, Nfpvarid;
   int dimids[2], sfildims[2], mfildims[2], ffildims[2], alpdims[2];
   int singledims[2], ncoildims[2];
 
   Nseg = Nseg+1;
   
   nc_create(FILE_NAME, NC_CLOBBER, &ncid); 
   nc_def_dim(ncid, "phony_dim_1", 1, &p1dimid); //Just 1
   nc_def_dim(ncid, "phony_dim_2", Nfils*Nseg, &nsegfildimid);
   nc_def_dim(ncid, "phony_dim_3", 2*NFalpha+1, &alpampdimid);
   nc_def_dim(ncid, "phony_dim_4", Nseg*5, &ffildimid);

   nc_def_dim(ncid, "Nzeta", Nzeta, &xdimid);
   nc_def_dim(ncid, "Nteta", Nteta, &ydimid);
   nc_def_dim(ncid, "Nseg", Nseg, &nsegdimid);
   nc_def_dim(ncid, "Ncoil", Ncoil, &coildimid);
   nc_def_dim(ncid, "iCoil", iCoil, &symcoildimid);

   dimids[0] = xdimid;
   dimids[1] = ydimid;
   
   sfildims[0] = coildimid;
   sfildims[1] = nsegdimid;

   mfildims[0] = coildimid;
   mfildims[1] = nsegfildimid;

   alpdims[0] = symcoildimid;
   alpdims[1] = alpampdimid;

   ffildims[0] = coildimid;
   ffildims[1] = ffildimid;

   singledims[0] = p1dimid;
   singledims[1] = p1dimid;

   ncoildims[0] = coildimid;
   ncoildims[1] = p1dimid;

   //Output the surface
   nc_def_var(ncid, "xsurf", NC_DOUBLE, 2, dimids, &xvarid);
   nc_def_var(ncid, "ysurf", NC_DOUBLE, 2, dimids, &yvarid);
   nc_def_var(ncid, "zsurf", NC_DOUBLE, 2, dimids, &zvarid);  

   //Output the single filaments xyz
   nc_def_var(ncid, "sx", NC_DOUBLE, 2, sfildims, &sxvarid); 
   nc_def_var(ncid, "sy", NC_DOUBLE, 2, sfildims, &syvarid); 
   nc_def_var(ncid, "sz", NC_DOUBLE, 2, sfildims, &szvarid); 

   nc_def_var(ncid, "cx", NC_DOUBLE, 2, ncoildims, &cxvarid); 
   nc_def_var(ncid, "cy", NC_DOUBLE, 2, ncoildims, &cyvarid); 
   nc_def_var(ncid, "cz", NC_DOUBLE, 2, ncoildims, &czvarid); 

   nc_def_var(ncid, "tx", NC_DOUBLE, 2, sfildims, &txvarid); 
   nc_def_var(ncid, "ty", NC_DOUBLE, 2, sfildims, &tyvarid); 
   nc_def_var(ncid, "tz", NC_DOUBLE, 2, sfildims, &tzvarid); 

   nc_def_var(ncid, "nx", NC_DOUBLE, 2, sfildims, &nxvarid); 
   nc_def_var(ncid, "ny", NC_DOUBLE, 2, sfildims, &nyvarid); 
   nc_def_var(ncid, "nz", NC_DOUBLE, 2, sfildims, &nzvarid); 

   nc_def_var(ncid, "bx", NC_DOUBLE, 2, sfildims, &bxvarid); 
   nc_def_var(ncid, "by", NC_DOUBLE, 2, sfildims, &byvarid); 
   nc_def_var(ncid, "bz", NC_DOUBLE, 2, sfildims, &bzvarid); 

   //Output the multi filaments xyz
   nc_def_var(ncid, "mx", NC_DOUBLE, 2, mfildims, &mxvarid); 
   nc_def_var(ncid, "my", NC_DOUBLE, 2, mfildims, &myvarid); 
   nc_def_var(ncid, "mz", NC_DOUBLE, 2, mfildims, &mzvarid); 

   //Output finite build coil surf xyz
   nc_def_var(ncid, "fx", NC_DOUBLE, 2, ffildims, &fxvarid); 
   nc_def_var(ncid, "fy", NC_DOUBLE, 2, ffildims, &fyvarid); 
   nc_def_var(ncid, "fz", NC_DOUBLE, 2, ffildims, &fzvarid); 


   //Output the alphss
   nc_def_var(ncid, "alpha", NC_DOUBLE, 2, alpdims, &alpvarid);

   //Output the centroid single filament field
/*   nc_def_var(ncid, "sB", NC_DOUBLE, 2, dimids, &sbvarid);  
   nc_def_var(ncid, "sBx", NC_DOUBLE, 2, dimids, &sbxvarid);
   nc_def_var(ncid, "sBy", NC_DOUBLE, 2, dimids, &sbyvarid);
   nc_def_var(ncid, "sBz", NC_DOUBLE, 2, dimids, &sbzvarid);  
   nc_def_var(ncid, "sBn", NC_DOUBLE, 2, dimids, &sbnvarid);  
*/

   //Output the final multi filament field
   nc_def_var(ncid, "mB", NC_DOUBLE, 2, dimids, &mbvarid);  
   nc_def_var(ncid, "mBx", NC_DOUBLE, 2, dimids, &mbxvarid);
   nc_def_var(ncid, "mBy", NC_DOUBLE, 2, dimids, &mbyvarid);
   nc_def_var(ncid, "mBz", NC_DOUBLE, 2, dimids, &mbzvarid);  
   nc_def_var(ncid, "mBn", NC_DOUBLE, 2, dimids, &mbnvarid);  

   //Output some other useful stuff
   nc_def_var(ncid, "NFalpha", NC_INT, 2, singledims, &NFalphavarid); 
   nc_def_var(ncid, "niter", NC_INT, 2, singledims, &nitervarid);
   nc_def_var(ncid, "len_rad", NC_DOUBLE, 2, singledims, &len_radvarid);
   nc_def_var(ncid, "len_tor", NC_DOUBLE, 2, singledims, &len_torvarid);
   nc_def_var(ncid, "Nradfil", NC_INT, 2, singledims, &Nradfilvarid); 
   nc_def_var(ncid, "Ntorfil", NC_INT, 2,singledims, &Ntorfilvarid);
   nc_def_var(ncid, "Nfp", NC_INT, 2,singledims, &Nfpvarid);
 
 
   
   nc_enddef(ncid);
 
   nc_put_var_double(ncid, xvarid, &xsurf[0]);
   nc_put_var_double(ncid, yvarid, &ysurf[0]);
   nc_put_var_double(ncid, zvarid, &zsurf[0]);

   nc_put_var_double(ncid, sxvarid, &sfilx[0]);
   nc_put_var_double(ncid, syvarid, &sfily[0]);
   nc_put_var_double(ncid, szvarid, &sfilz[0]);

   nc_put_var_double(ncid, cxvarid, &cx[0]);
   nc_put_var_double(ncid, cyvarid, &cy[0]);
   nc_put_var_double(ncid, czvarid, &cz[0]);

   nc_put_var_double(ncid, txvarid, &tx[0]);
   nc_put_var_double(ncid, tyvarid, &ty[0]);
   nc_put_var_double(ncid, tzvarid, &tz[0]);

   nc_put_var_double(ncid, nxvarid, &nx[0]);
   nc_put_var_double(ncid, nyvarid, &ny[0]);
   nc_put_var_double(ncid, nzvarid, &nz[0]);

   nc_put_var_double(ncid, bxvarid, &bx[0]);
   nc_put_var_double(ncid, byvarid, &by[0]);
   nc_put_var_double(ncid, bzvarid, &bz[0]);

   nc_put_var_double(ncid, mxvarid, &mfilx[0]);
   nc_put_var_double(ncid, myvarid, &mfily[0]);
   nc_put_var_double(ncid, mzvarid, &mfilz[0]);

   nc_put_var_double(ncid, fxvarid, &ffilx[0]);
   nc_put_var_double(ncid, fyvarid, &ffily[0]);
   nc_put_var_double(ncid, fzvarid, &ffilz[0]);

   nc_put_var_double(ncid, alpvarid, &alpamps[0]);
/*
   nc_put_var_double(ncid, sbvarid,  &Bsfil[0]);
   nc_put_var_double(ncid, sbxvarid, &Bsfilx[0]);
   nc_put_var_double(ncid, sbyvarid, &Bsfily[0]);
   nc_put_var_double(ncid, sbzvarid, &Bsfilz[0]);
   nc_put_var_double(ncid, sbnvarid, &Bsfiln[0]);
*/ 
   nc_put_var_double(ncid, mbvarid,  &Bmfil[0]);
   nc_put_var_double(ncid, mbxvarid, &Bmfilx[0]);
   nc_put_var_double(ncid, mbyvarid, &Bmfily[0]);
   nc_put_var_double(ncid, mbzvarid, &Bmfilz[0]);
   nc_put_var_double(ncid, mbnvarid, &Bmfiln[0]);
   
   nc_put_var_int(ncid, NFalphavarid, &NFalpha); 
   nc_put_var_double(ncid, len_radvarid, &len_rad); 
   nc_put_var_double(ncid, len_torvarid, &len_tor); 
   nc_put_var_int(ncid, Nradfilvarid, &Nradfil); 
   nc_put_var_int(ncid, Ntorfilvarid, &Ntorfil); 
   nc_put_var_int(ncid, nitervarid, &niter); 
   nc_put_var_int(ncid, Nfpvarid, &Nfp); 

   nc_close(ncid);
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
