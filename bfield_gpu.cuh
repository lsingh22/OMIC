#ifndef BFIELD_GPU_CUH
#define BFIELD_GPU_CUH

#include <cuda.h>
#include <cuda_runtime.h>
#include "globals.h"

__global__ void field_kernel0(const double* mx, const double* my, const double* mz,
                             const double* coilcurrents, double x, double y, double z, 
                             double* dbx, double* dby, double* dbz, 
                             const double* rot,  const int nfp, const int maxsegs);

__global__ void field_kernel1(const double* mx, const double* my, const double* mz,
                             const double* coilcurrents, double x, double y, double z, 
                             double* dbx, double* dby, double* dbz, 
                             const double* rot,  const int nfp, const int maxsegs);

__global__ void parallel_reduction(double* dBx, double* dBy, double* dBz, double* bx, double* by, double* bz);

#ifdef __cplusplus
extern "C" {
#endif

void magnetic_field(const double* mfx, const double* mfy, const double* mfz,
                    double* xs, double* ys, double* zs, const double* currents, 
                    const int ncoil, const int nfil, const int nfp, 
                    const int nseg, const int size_fp);

//void CalculateFieldParallelGPU(void); 

#ifdef __cplusplus
}
#endif

//void CalculateFieldParallelGPU(void);

#endif
