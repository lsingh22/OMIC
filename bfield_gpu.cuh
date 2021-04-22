#ifndef _BFIELD_GPU_CUH
#define _BFIELD_GPU_CUH

__global__ void field_kernel(void);

__global__ void magnetic_field(double* bx, double* by, double* bz, double* currents, const double* mfx, const double* mfy, const double* mfz, const int ncoil, const int nseg);

#endif
