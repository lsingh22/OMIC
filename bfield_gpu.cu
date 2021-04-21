
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

__global__ void magnetic_field(double* bx, double* by, double* bz, const double* currents, const double* mfx, const double* mfy, const int nseg) {

	// Dynamically allocate shared memory

	// Set up execution configuration

	// For each point on magnetic surface,

		// Call field kernel to get magnetic field 

	// Call rotation kernel to reflect results to all field periods

	// Cleanup
}

