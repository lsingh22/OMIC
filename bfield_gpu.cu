#include <cstdio>
#include <cuda.h>
#include <cuda_runtime.h>
#include "bfield_gpu.cuh"
#include "globals.h"
#include <stdio.h>

#define MY_PI 3.14159265358979323846

__global__ void field_kernel0(const double* mx, const double* my, const double* mz, 
									   const double* curr, const double x, const double y, 
                              const double z, double* dBx, double* dBy, double* dBz,
                              const double* Rot, const int nfp, const int MAXSEG) {

	// Dynamically allocate shared memory, positions and rotaiton matrix
   extern __shared__ double posx[];
   extern __shared__ double posy[];
   extern __shared__ double posz[];
   extern __shared__ double Rcs[];

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
      
         // DEBUG:
/*         if(tx==0 && threadIdx.y==0){
            for(int i = 0; i < blockDim.x * blockDim.y; i++) {
               printf("mx: %f my: %f mz: %f\n", mx[i], my[i], mz[i]);
               printf("bx: %f by: %f bz: %f\n", dBx[i], dBy[i], dBz[i]);
            }
         }      
*/    }
   }
}

// their size in the looped kernel0 implementation.
__global__ void field_kernel1(const double* mx, const double* my, const double* mz, 
									   const double* curr, double x, double y, double z, 
                              double* dBx, double* dBy, double* dBz,
                              const double* Rot, const int nfp, const int MAXSEG) {

	// Dynamically allocate shared memory, positions and rotaiton matrix
   extern __shared__ double posx[];
   extern __shared__ double posy[];
   extern __shared__ double posz[];
   extern __shared__ double Rcs[];

   // Define indexing variables 
   int tx = threadIdx.x;
   int tz = threadIdx.z;
   int tid    = threadIdx.y * blockDim.x + threadIdx.x;
   int posIdx = blockIdx.x * blockDim.y * blockDim.x + tid;

   // Doubles needed for integration
   double x1, y1, z1, x2, y2, z2;
   double ri, rf, l, ex, ey, ez;
   double bx, by, bz, coef;

   // Store current of the relevant coil (filament)
   double cur = curr[threadIdx.y];

   // If segment is valid, calculate contribution to B field
   if(posIdx < MAXSEG) {

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
         
         // Integration loop: one thread sums contributions from all points

         // Locate stellarator symmetric point on magnetic boundary
         double tempx = x;
         x =  x * Rcs[tz] + y * Rcs[nfp + tz];
         y = -tempx * Rcs[nfp + tz] + y * Rcs[tz];
         z = z; // not needed here, but needed for future speedup

         // Use H & H method to find field contribution due to segment
         x1 = x - x1;
         y1 = y - y1;
         z1 = z - z1;

         x2 = x - x2;
         y2 = y - y2;
         z2 = z - z2;

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

         dBx[nfp * posIdx + tz] = bx * Rcs[tz] - by * Rcs[nfp + tz];
         dBy[nfp * posIdx + tz] = bx * Rcs[nfp + tz] + by * Rcs[tz];
         dBz[nfp * posIdx + tz] = bz;
      }
   }
}

__global__ void reduce_kernel(double* g_idata, double* g_odata, unsigned int n) {

   // Dynamically allocate shared memory 
   extern __shared__ double sdata[];
   
   // Useful indexing variables adapted from slides   
   unsigned int tid = threadIdx.x;
   unsigned int i   = blockIdx.x * (blockDim.x * 2) + threadIdx.x;
   
   // First add during load 
   // Check bounds before writing sum to shared memory 

   sdata[tid] = 0;
   if(i < n) {
      sdata[tid] = (i + blockDim.x < n) ? g_idata[i] + g_idata[i + blockDim.x] : g_idata[i];
   }  
   
   // Give all threads access to the shared memory  
   __syncthreads();
   
   // For the input length of the array, determine sum for given thread 
   // Reduce in the shared memory 
   if (i < n) {
      for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
         if (tid < s) { 
            sdata[tid] += sdata[tid + s];
         }      
         __syncthreads();
      }  
      
      // Write block result to global memory 
      if (tid == 0) {
         g_odata[blockIdx.x] = sdata[0];
      }  
   }  
}  

double cosnfp(int ip, int nfp) { 

   // Cosine of vector component at the ip-th field period
   return cos((ip - 1) * MY_PI * 2 / nfp);
}

double sinnfp(int ip, int nfp) { 
   
   // Sine of vector component at the ip-th field period
   return sin((ip - 1) * MY_PI * 2 / nfp);
}

// Helper function for parallel reduction
void swap(double** a, double** b) {

   double* temp = *a;
   *a = *b;
   *b = temp;
}

extern "C" void magnetic_field(const double* mfx, const double* mfy, const double* mfz, 
                               double* xsurf, double* ysurf, double* zsurf, 
                               const double* current, const int ncoil, const int nfil, 
                               const int nfp, const int nseg, const int size_fp) {

	// Define execution configuration for field_kernel call
	// A 2D block configuration enables easy tracking of coils, segments
	// A 1D grid configuration is determined by ty, the number of coils / block
	// Field kernel computes B-field due to tx-th segment of ty-th coil
	// Shared memory holds coil positions, currents, and two reduction buffers
  
   // Change this for lower CUDA CC hardware
   unsigned int MAX_BLOCK_THREADS = 1024;

	unsigned int coils_per_block = MAX_BLOCK_THREADS / nseg - 1; // TODO: for some reason, this -1 is needed right now
	unsigned int segs_per_block  = coils_per_block * nseg;   
	dim3 DimBlock(nseg, coils_per_block);
	unsigned int DimGrid = (ncoil + coils_per_block - 1) / coils_per_block;
	unsigned int field_shmem_size = (segs_per_block * 3 + 2 * nfp) * sizeof(double);
   unsigned int ndb = ncoil * nseg - 1; // TODO: there will be ncoil empty entries with this definition

   // Define maximum number of segments for kernel memory safety
   unsigned int max = nseg * ncoil;

	// Create helper currents and reduction buffer unified arrays
	double* dBx = (double*) malloc(ndb * sizeof(double));
	double* dBy = (double*) malloc(ndb * sizeof(double));
	double* dBz = (double*) malloc(ndb * sizeof(double));

   double *mx, *my, *mz, *xs, *ys, *zs, *Ic, *R;
 
   // Allocate managed memory for coil segs/currents and magnetic surface
   cudaMallocManaged((void**)&mx, ncoil * nseg * sizeof(double));
 	cudaMallocManaged((void**)&my, ncoil * nseg * sizeof(double));
	cudaMallocManaged((void**)&mz, ncoil * nseg * sizeof(double));
	cudaMallocManaged((void**)&xs, size_fp * sizeof(double));
	cudaMallocManaged((void**)&ys, size_fp * sizeof(double));
	cudaMallocManaged((void**)&zs, size_fp * sizeof(double));

	cudaMallocManaged((void**)&Ic,  ncoil   * sizeof(double));
   cudaMallocManaged((void**)&R,   nfp * 2 * sizeof(double));
	cudaMallocManaged((void**)&dBx, ndb * sizeof(double));
	cudaMallocManaged((void**)&dBy, ndb * sizeof(double));
	cudaMallocManaged((void**)&dBz, ndb * sizeof(double));
   
   // Copy CPU globals to managed memory arrays
   memcpy(mx, mfx, ncoil * nseg * sizeof(double));
   memcpy(my, mfy, ncoil * nseg * sizeof(double));
   memcpy(mz, mfz, ncoil * nseg * sizeof(double));

   memcpy(xs, xsurf, size_fp * sizeof(double));
   memcpy(ys, ysurf, size_fp * sizeof(double));
   memcpy(zs, zsurf, size_fp * sizeof(double));

   // Fill Ic helper array with values from currents
   // Holds the multifilament current scaled by mu_0/pi
   for(int i = 0; i < ncoil / nfil; i++) {
      for(int j = 0; j < nfil; j++) {
         Ic[i * nfil + j] = 1.0e-7 * current[i] / nfil;
      }
   }

   // Fill shared rotation matrix
   for(int i = 0; i < nfp; i++) {
      R[i] = cosnfp(i + 1, nfp); 
      R[nfp + i] = sinnfp(i + 1, nfp);
   }
   
   // Helper arrays and variables for parallel reduction
   unsigned int N;
   unsigned int threads_per_block;
   unsigned int blocks_per_grid;
   double* outputx;
   double* outputy;
   double* outputz;
   cudaMallocManaged((void**)&outputx, ncoil * nseg * sizeof(double));
   cudaMallocManaged((void**)&outputy, ncoil * nseg * sizeof(double));
   cudaMallocManaged((void**)&outputz, ncoil * nseg * sizeof(double));

   // Calculate segment dB contributions (field_kernel),
   // Perform parallel reduction to superpose dB (parallel_reduction)
   // TODO: openmp may help here, too
   for(int i = 0; i < size_fp; i++) {
      N = ncoil * nseg - 1;
      threads_per_block = MAX_BLOCK_THREADS;
      blocks_per_grid = (N + 2 * threads_per_block - 1) / (2 * threads_per_block);

		// Call field kernel to get magnetic field
      field_kernel0<<<DimGrid, DimBlock, field_shmem_size>>>(mx, my, mz, Ic, 
                                                             xs[i], ys[i], zs[i],
                                                             dBx, dBy, dBz, R, nfp, max);
      cudaDeviceSynchronize();
      while (N > 1) {
         reduce_kernel<<<blocks_per_grid, threads_per_block, threads_per_block * sizeof(double)>>>(dBx, outputx, N);
         reduce_kernel<<<blocks_per_grid, threads_per_block, threads_per_block * sizeof(double)>>>(dBy, outputy, N);
         reduce_kernel<<<blocks_per_grid, threads_per_block, threads_per_block * sizeof(double)>>>(dBz, outputz, N);

         N = blocks_per_grid;
         blocks_per_grid = (N + threads_per_block * 2 - 1) / (2 * threads_per_block);
         
         // Swap input and output
         if(N != 1) {
            swap(&dBx, &outputx);
            swap(&dBy, &outputy);
            swap(&dBz, &outputz);
         }
      }
      cudaDeviceSynchronize();
   }
	
	// Cleanup
	cudaFree(Ic);
	cudaFree(dBx);
	cudaFree(dBy);
	cudaFree(dBz);
}
