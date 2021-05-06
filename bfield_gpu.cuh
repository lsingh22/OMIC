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

#ifdef __cplusplus
extern "C" {
#endif

void magnetic_field0(const double* mfx, const double* mfy, const double* mfz,
                    double* xs, double* ys, double* zs, const double* currents, 
                    const int ncoil, const int nfil, const int nfp, 
                    const int nseg, const int size_fp, 
                    double* Bx, double* By, double* Bz, const int start, const int end);

void magnetic_field1(const double* mfx, const double* mfy, const double* mfz,
                    double* xs, double* ys, double* zs, const double* currents, 
                    const int ncoil, const int nfil, const int nfp, 
                    const int nseg, const int size_fp, 
                    double* Bx, double* By, double* Bz, const int start, const int end);

#ifdef __cplusplus
}
#endif

#endif
