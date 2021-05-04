#include <cstdio>
#include <cuda.h>
#include <cuda_runtime.h>
#include "bfield_gpu.cuh"

#define MY_PI 3.14159265358979323846

__global__ void field_kernel(const double* mx, const double* my, const double* mz, 
									  const double* curr, const double x, const double y, 
                             const double z, double* dBx, double* dBy, double* dBz,
                             const double* Rot, const int nfp, const int MAXSEG) {

	// Dynamically allocate shared memory, positions and rotaiton matrix
   extern volatile __shared__ double posx[];
   extern volatile __shared__ double posy[];
   extern volatile __shared__ double posz[];
   extern volatile __shared__ double Rcs[];

   // Define indexing variables 
   int tx = threadIdx.x;
   int tid    = threadIdx.y * blockDim.x + threadIdx.x;
   int posIdx = blockIdx.x * blockDim.y * blockDim.x + tid;

   // Doubles needed for integration
   double x1, y1, z1, x2, y2, z2;
   double xx, yy, zz, ri, rf, l, ex, ey, ez;
   double bx, by, bz, bxx, byy, bzz, coef;

   // Store current of the relevant coil (filament)
   double cur = curr[threadIdx.y];

   // If segment is valid, calculate contribution to B field
   if(posIdx < MAXSEG) {

      // DEBUG:
      //printf("bid: %d ty: %d tx: %d posIdx: %d bid*ty*tx+tx\n", blockIdx.x, threadIdx.y, threadIdx.x, posIdx);

      // Write rotation matrix to shared memory
      if(tx < nfp) {
         Rcs[tx] = Rot[tx];
         Rcs[nfp + tx] = Rot[tx + nfp];          
      }
      
      // Global load of start of segment
      x1 = mx[posIdx]; 
      y1 = my[posIdx];
      z1 = mz[posIdx];

      // Write starting segments to shared memory 
      posx[tid] = x1;
      posy[tid] = y1;
      posz[tid] = z1;

      // Share segment starts to hide latency in loading neighbor segment
      // Also share the rotation matrix
      __syncthreads();

      // Prevent segment index from reaching out of bounds at end of coil
      if(tx < blockDim.x - 1) {
   
         // Shared load of end of segment
         x2 = posx[tid + 1]; 
         y2 = posy[tid + 1];
         z2 = posz[tid + 1]; 
         
         // Accumulators over symmetric field periods
         bxx = 0.; byy = 0.; bzz = 0.;

         // Integration loop: one thread sums contributions from all points
         // For further development, this could be split up over nfp threads
         for(int ip = 0; ip < nfp; ip++) {

            // Locate stellarator symmetric point on magnetic boundary
            xx =  x * Rcs[ip] + y * Rcs[nfp + ip];
            yy = -x * Rcs[nfp + ip] + y * Rcs[ip];
            zz = z; // not needed here, but needed for future speedup

            // Use H & H method to find field contribution due to segment
            x1 = xx - x1;
            y1 = yy - y1;
            z1 = zz - z1;

            x2 = xx - x2;
            y2 = yy - y2;
            z2 = zz - z2;

            ri = sqrt(x1 * x1 + y1 * y1 + z1 * z1);  
            rf = sqrt(x2 * x2 + y2 * y2 + z2 * z2);

            l = sqrt((x2-x1) * (x2-x1) + (y2-y1) * (y2-y1) + (z2-z1) * (z2-z1));  

            ex = (x2-x1)/l;  
            ey = (y2-y1)/l;  
            ez = (z2-z1)/l;  

            // Prefactor for cross product 
            coef = cur * l * (ri+rf) / (ri*rf * ((ri+rf) * (ri+rf) - l * l));   

            bx = coef * (ey * z1 - ez * y1);
            by = coef * (ez * x1 - ex * z1);
            bz = coef * (ex * y1 - ey * x1);
        
            bxx += bx * Rcs[ip] - by * Rcs[nfp + ip];
            byy += bx * Rcs[nfp + ip] + by * Rcs[ip];
            bzz += bz;    
         }

         // Write segment dB to output   
         dBx[posIdx] = bxx;
         dBy[posIdx] = byy;
         dBz[posIdx] = bzz;
      }
   }
}

__global__ void parallel_reduction(double* dBx, double* dBy, double* dBz, 
								           double* bx,  double* by,  double* bz);

double cosnfp(int ip) { 

   // Cosine of vector component at the ip-th field period
   return cos((ip - 1) * MY_PI * 2 / Nfp);
}

double sinnfp(int ip) { 
   
   // Sine of vector component at the ip-th field period
   return sin((ip - 1) * MY_PI * 2 / Nfp);
}

extern "C" void magnetic_field(const double* mfx, const double* mfy, const double* mfz, 
                               const double* current, const int ncoil, const int nfil,
                               const int nfp, const int nseg, const int size_fp) {

	// Define execution configuration for field_kernel call
	// A 2D block configuration enables easy tracking of coils, segments
	// A 1D grid configuration is determined by ty, the number of coils / block
	// Field kernel computes B-field due to tx-th segment of ty-th coil
	// Shared memory holds coil positions, currents, and two reduction buffers
  
   // Change this for lower CUDA CC hardware
   unsigned int MAX_BLOCK_THREADS = 1024;

	unsigned int coils_per_block = MAX_BLOCK_THREADS / nseg;
	unsigned int segs_per_block  = coils_per_block * nseg;   
	dim3 DimBlock(coils_per_block, nseg);
	unsigned int DimGrid = (ncoil + coils_per_block - 1) / coils_per_block;
	unsigned int field_shmem_size = (segs_per_block * 3 + 2 * nfp) * sizeof(double);
   unsigned int ndb = ncoil * (nseg - 1);

	// Create helper currents and reduction buffer unified arrays
	double* Ic = (double*) malloc(ncoil * sizeof(double));
   double* R = (double*) malloc(nfp * 2 * sizeof(double));
	double* dBx = (double*) malloc(ndb * sizeof(double));
	double* dBy = (double*) malloc(ndb * sizeof(double));
	double* dBz = (double*) malloc(ndb * sizeof(double));

   // Define maximum number of segments for kernel memory safety
   unsigned int max = nseg * ncoil;

	// Fill Ic helper array with values from currents
   for(int i = 0; i < ncoil / nfil; i++) {
      for(int j = 0; j < nfil; j++) {
         Ic[i * nfil + j] = current[i];
      }
   }

   // Fill shared rotation matrix
   for(int i = 0; i < nfp; i++) {
      R[i] = cosnfp(i + 1); 
      R[nfp + i] = sin(i + 1);
   }
      
   //printf("coils_per_block and nseg are: %d and %d\n", coils_per_block, nseg);
   // Allocate managed memory for coil segs/currents and magnetic surface
   cudaMallocManaged((void**)&mfilx, ncoil * nseg * sizeof(double));
 	cudaMallocManaged((void**)&mfily, ncoil * nseg * sizeof(double));
	cudaMallocManaged((void**)&mfilz, ncoil * nseg * sizeof(double));
	cudaMallocManaged((void**)&xsurf, size_fp * sizeof(double));
	cudaMallocManaged((void**)&ysurf, size_fp * sizeof(double));
	cudaMallocManaged((void**)&zsurf, size_fp * sizeof(double));

   // Allocate managed memory for the helper arrays
	cudaMallocManaged((void**)&Ic,  ncoil   * sizeof(double));
   cudaMallocManaged((void**)&R,   nfp * 2 * sizeof(double));
	cudaMallocManaged((void**)&dBx, ndb * sizeof(double));
	cudaMallocManaged((void**)&dBy, ndb * sizeof(double));
	cudaMallocManaged((void**)&dBz, ndb * sizeof(double));

   field_kernel<<<DimGrid, dim3(129,7), field_shmem_size>>>(mfilx, mfily, mfilz, Ic, 
                                                     xsurf[0], ysurf[0], zsurf[0],
                                                     dBx, dBy, dBz, R, nfp, max);

   // Exit if there is a problem with kernel execution	
   cudaError_t err = cudaGetLastError();
   
   if(err != cudaSuccess) {
      printf("CUDA Error: %s\n", cudaGetErrorString(err));
      exit(-1);
   }

   cudaDeviceSynchronize();

   // Calculate segment dB contributions (field_kernel)
   // Perform parallel reduction to superpose dB (parallel_reduction)
	for(int i = 0; i < size_fp; i++) {
		// Call field kernel to get magnetic field
/*		field_kernel<<<DimGrid, dim3(6,129), field_shmem_size>>>(mfilx, mfily, mfilz, Ic, 
                                                            xsurf[i], ysurf[i], zsurf[i],
                                                            dBx, dBy, dBz, R, nfp, max);

      cudaDeviceSynchronize(); // synchronize to prepare for reduction
*/		//parallel_scan<<<1, DimGrid, segs_per_block>>>(dBx, dBy, dBz, bx+i, by+i, bz+i);
	}

	// Call rotation kernel to reflect results to all field periods
	
	// Cleanup
	cudaFree(Ic);
	cudaFree(dBx);
	cudaFree(dBy);
	cudaFree(dBz);
}

// Driving code for GPU calculation 
extern "C" void CalculateFieldParallelGPU(void) {

	// Set up timing events
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	float ms;

	// Time call to magnetic field function using CUDA events
	cudaEventRecord(start, 0);
	magnetic_field(mfilx, mfily, mfilz, currents, Ncoil * Nfils, Nfils, Nfp, Nseg+1, size_fp);  
   cudaEventRecord(stop, 0);

	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&ms, start, stop);

	// Cleanup
	cudaEventDestroy(start);
	cudaEventDestroy(stop);		
}
