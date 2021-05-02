#include <cuda.h>
#include <cuda_runtime.h>
#include "bfield_gpu.cuh"

__global__ void field_kernel(const double* mx, const double* my, const double* mz, 
									  const double* curr, const double x, const double y, 
                             const double z, double* dbx, double* dby, double* dbz) {

	// Dynamically allocate shared memory for inputs, outputs
   // Inputs are coil positions, currents, and rotation matrix
   // Outputs are field contributions
   extern volatile __shared__ double posx[];
   extern volatile __shared__ double posy[];
   extern volatile __shared__ double posz[];
   extern volatile __shared__ double icoil[];
   extern volatile __shared__ double rot_matrix[];
   extern volatile __shared__ double shdbx[];
   extern volatile __shared__ double shdby[];
   extern volatile __shared__ double shdbz[];

   // Define indexing variables 
   int bid = blockIdx.x;
   int seg = threadIdx.x; // for keeping track of segments
   int cid = threadIdx.y; // for indexing current 
   int idx = bid * seg * cid + seg * cid + seg; // position index

   // Global load of start of segment
   double xi = mx[cid * seg + seg]; 
   double yi = my[cid * seg + seg];
   double zi = mz[cid * seg + seg];

   // Write starting segments to shared memory 
   posx[cid * seg + seg] = xi;
   posy[cid * seg + seg] = yi;
   posz[cid * seg + seg] = zi;

   // Share segment starts to hide latency in loading neighbor segment
	__syncthreads();

   // Shared load of end of segment
   double xf = posx[cid * seg + seg + 1]; 
   double yf = posy[cid * seg + seg + 1];
   double zf = posz[cid * seg + seg + 1]; 
	
   // Using start and end segments, and kernel arguments, compute cartesian dB
   // Parameters are defined as in Hanson & Hirschman (TODO:ref)

   // for loop over nfp
   // for(int ip = 1; i < nfp; i++) {


      // update local bx,by,bz
   // }
      
   // Load segment dB to shared memory to perform block reduction
   // Swapped double buffer method is performed for each cartesian component
   //shdbx[cid * seg + seg] = bx;
   //shdby[cid * seg + seg] = by;
   //shdbz[cid * seg + seg] = bz;

	// __syncthreads();

   // Perform Hillis-Steele reduction and store in global db (x,y,z)   

}

__global__ void hillis_steele(double* dBx, double* dBy, double* dBz, 
								  double* bx,  double* by,  double* bz);

extern "C" void magnetic_field(const double* mfx, const double* mfy, const double* mfz, 
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
// unsigned int scan_shMem_size = 2 * DimGrid * 3 * sizeof(double); 

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
		field_kernel<<<DimGrid, DimBlock, field_shMem_size>>>(mfilx, mfily, mfilz, Icoil, 
                                                            xsurf[i], ysurf[i], zsurf[i],
                                                            dBx, dBy, dBz);
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
