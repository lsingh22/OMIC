#ifndef _BFIELD_GPU_CUH
#define _BFIELD_GPU_CUH

__global__ void field_kernel(double* mx, double* my, double* mz, double* coilcurrents, double* dbx, double* dby, double* dbz);

__global__ void hillis_steele(double* dBx, double* dBy, double* dBz, double* bx, double* by, double* bz);

__global__ void magnetic_field(double* bx, double* by, double* bz, double* currents, const double* mfx, const double* mfy, const double* mfz, const int ncoil, const int nseg);

#endif
