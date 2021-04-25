#include <cuda.h>
#include <cuda_runtime.h>
#include "bfield_gpu.cuh"

__global__ void field_kernel(const double* mx, const double* my, const double* mz, 
									  const double* curr, double* dbx, double* dby, double* dbz) {

	// Dynamically allocate shared memory

	// Loads from constant memory

	// Each thread loads in a filament position to shMem

	// __syncthreads();

	// Calculate start and end vectors and write to shMem

	// __syncthreads();

	// Calculate segment length, unit vector, and cross product prefactor

	// Calculate magnetic field contribution and store in shMem

}

__global__ void hillis_steele(double* dBx, double* dBy, double* dBz, 
								  double* bx,  double* by,  double* bz);

__host__ void magnetic_field(const double* mfx, const double* mfy, const double* mfz, 
                             const double* currents, const int ncoil, const int nseg, const int size_fp) {

	// TODO: for purposes of parallel implementation, ncoil should be Ncoil * Nfils

	// Define execution configuration for field_kernel call
	// A 2D block configuration enables easy tracking of coils, segments
	// A 1D grid configuration is determined by ty, the number of coils / block
	// Field kernel computes B-field due to tx-th segment of ty-th coil
	// Shared memory holds coil positions, currents, and two reduction buffers

	unsigned int coils_per_block = (1024 + nseg - 1) / nseg;
	unsigned int segs_per_block  = coils_per_block * nseg;   

	dim3 DimBlock(coils_per_block, segs_per_block);
	unsigned int DimGrid = (ncoil + coils_per_block - 1) / coils_per_block;
	unsigned int field_shMem_size = (segs_per_block * 3 + 2 * (segs_per_block - 1) * 3 \
										      + coils_per_block) * sizeof(double);
	unsigned int scan_shMem_size = 2 * DimGrid * 3 * sizeof(double); 

	// Create helper currents and reduction buffer unified arrays
	double* Icoil;
	double* dBx;
	double* dBy;
	double* dBz;

	// TODO: fill Icoil array with values from currents
	int flags = 0;
	 
	cudaMallocManaged((void**)&Icoil, ncoil * sizeof(double), flags);
	cudaMallocManaged((void**)&dBx, DimGrid * sizeof(double), flags);
	cudaMallocManaged((void**)&dBy, DimGrid * sizeof(double), flags);
	cudaMallocManaged((void**)&dBz, DimGrid * sizeof(double), flags);

	// For each point on magnetic surface,
	for(int i = 0; i < size_fp; i++) {
		
		// Call field kernel to get magnetic field
		field_kernel<<<DimGrid, DimBlock, field_shMem_size>>>(mfx, mfy, mfz, Icoil, dBx, dBy, dBz);
		cudaDeviceSynchronize(); // synchronize to prepare for reduction
		
		//Call Hillis_Steele kernel to scan the partial sums and get bx, by, bz
		//hillis_steele<<<1, DimGrid, scan_shMem_size>>>(dBx, dBy, dBz, bx+i, by+i, bz+i);
	}

	// TODO: Call rotation kernel to reflect results to all field periods
	
	// Cleanup
	cudaFree(Icoil);
	cudaFree(dBx);
	cudaFree(dBy);
	cudaFree(dBz);
}

__host__ void MagneticFieldGPU(void) {

	// Allocate unified memory arrays for coil segs/currents and magnetic surface
	cudaMallocManaged((void**)&mfilx, Ncoil * Nfils * (Nseg+1) * sizeof(double));
 	cudaMallocManaged((void**)&mfily, Ncoil * Nfils * (Nseg+1) * sizeof(double));
	cudaMallocManaged((void**)&mfilz, Ncoil * Nfils * (Nseg+1) * sizeof(double));
	cudaMallocManaged((void**)&currents, Ncoil * sizeof(double));

	cudaMallocManaged((void**)&xsurf, size_fp * sizeof(double));
	cudaMallocManaged((void**)&ysurf, size_fp * sizeof(double));
	cudaMallocManaged((void**)&zsurf, size_fp * sizeof(double));

	// Set up timing events
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	float ms;

	// Time call to magnetic field function using CUDA events
	cudaEventRecord(start, 0);
	magnetic_field(mfilx, mfily, mfilz, currents,  
						Ncoil * Nfils, Nseg+1, size_fp);  
	cudaEventRecord(stop, 0);

	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&ms, start, stop);

	// Cleanup
	cudaEventDestroy(start);
	cudaEventDestroy(stop);		
}
