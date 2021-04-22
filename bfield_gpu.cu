
#include <cuda.h>
#include "bfield_gpu.cuh"

__global__ void field_kernel(void) {

	// Loads from constant memory

	// Each thread loads in a filament position to shMem

	// __syncthreads();

	// Calculate start and end vectors and write to shMem

	// __syncthreads();

	// Calculate segment length, unit vector, and cross product prefactor

	// Calculate magnetic field contribution and store in shMem

}

__host__ void magnetic_field(double* bx, double* by, double* bz, const double* currents, const double* mfx, const double* mfy, const double* mfz, const int ncoil, const int nseg) {

	// TODO: for purposes of parallel implementation, ncoil should be Ncoil * Nfils

	// Define execution configuration for field_kernel call
	// A 2D block configuration enables easy tracking of coils, segments
	// A 1D grid configuration is determined by ty, the number of coils / block
	// Field kernel computes B-field due to tx-th segment of ty-th coil

	unsigned int ty_dim = (1024 + nseg - 1) / nseg; // number of coils per block
	unsigned int tx_dim = ty_dim * nseg; // number of segments per coil
	unsigned int DimGrid = (ncoil + ty_dim - 1) / ty_dim;

	dim3 DimBlock(tx_dim, ty_dim);

	// Create helper currents array that stores each filament's current
	double* Icoil;
	
	// Fill with currents of each filament

	cudaMallocManaged((void**)&Icoil, ncoil * sizeof(double));

	// For each point on magnetic surface,

		// Call field kernel to get magnetic field
		// field_kernel<<< >>>(mfx[i], mfy[i], mfz[i], Icoil);
		// cudaDeviceSynchronize(); // synchronize to prepare for reduction

	// Call a parallel reduction to get B-total at each point

	// Call rotation kernel to reflect results to all field periods

	// Cleanup
	cudaFree(Icoil);
}


