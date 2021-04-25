#ifndef BFIELD_GPU_CUH
#define BFIELD_GPU_CUH

#include <cuda.h>
#include <cuda_runtime.h>
#include "globals.h"

__global__ void field_kernel(const double* mx, const double* my, const double* mz,
                             const double* coilcurrents, double* dbx, double* dby, double* dbz);

__global__ void hillis_steele(double* dBx, double* dBy, double* dBz, double* bx, double* by, double* bz);

__host__ void magnetic_field(const double* mfx, const double* mfy, const double* mfz,
									  const double* currents, const int ncoil, const int nseg, const int size_fp);

__host__ void CalculateFieldParallelGPU(void);

#endif
