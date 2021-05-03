#include <cstdio>
#include <cuda.h>
#include <cuda_runtime.h>
#include "bfield_gpu.cuh"

#define MY_PI 3.14159265358979323846

__global__ void field_kernel(const double* mx, const double* my, const double* mz, 
									  const double* curr, const double x, const double y, 
                             const double z, double* dbx, double* dby, double* dbz,
                             const double* Rot, const int nfp, const int MAXSEG) {

	// Dynamically allocate shared memory for inputs, outputs
   // Inputs are coil positions, currents, and rotation matrix
   // Outputs are field contributions
   extern volatile __shared__ double posx[];
   extern volatile __shared__ double posy[];
   extern volatile __shared__ double posz[];
   extern volatile __shared__ double Rcs[];
   extern volatile __shared__ double shdbx[];
   extern volatile __shared__ double shdby[];
   extern volatile __shared__ double shdbz[];

   printf("EXECUTED\n");

   // Define indexing variables 
   int bid = blockIdx.x;
   int seg = threadIdx.x; // for keeping track of segments
   int cid = threadIdx.y; // for indexing current 
   int idx = bid * seg * cid + seg * cid + seg; // position index

   // Doubles needed for integration
   double x1, y1, z1, x2, y2, z2;
   double xx, yy, zz, ri, xi, yi, zi, rf, xf, yf, zf;
   double l, ex, ey, ez, bx, by, bz, bxx, byy, bzz;
   double coef;

   // Store current of the relevant coil (filament)
   double cur = curr[cid];

   // TODO: probably include in args the max size for if-else safety net here
   // if(idx < SEG_MAX) {

   // TODO: pass rotation matrix as global, then put in shared memory
   // Number of elements is 2 * Nfp, Nfp are cosines, Nfp are sines
   // if(idx < Nfp) {shared memory write from RMatrix to Rc, Rs} 

   // Global load of start of segment
   x1 = mx[cid * seg + seg]; 
   y1 = my[cid * seg + seg];
   z1 = mz[cid * seg + seg];

   // Write starting segments to shared memory 
   posx[cid * seg + seg] = x1;
   posy[cid * seg + seg] = y1;
   posz[cid * seg + seg] = z1;

   // Share segment starts to hide latency in loading neighbor segment
	__syncthreads();

   // Shared load of end of segment
   x2 = posx[cid * seg + seg + 1]; 
   y2 = posy[cid * seg + seg + 1];
   z2 = posz[cid * seg + seg + 1]; 
	
   // Using start and end segments, and kernel arguments, compute cartesian dB
   // Parameters are defined as in Hanson & Hirschman (TODO:ref)

   // Accumulators over symmetric field periods
   bxx = 0.; byy = 0.; bzz = 0.;

   for(int ip = 0; ip < nfp; ip++) {

      // Locate stellarator symmetric point on magnetic boundary
      xx =  x * Rcs[ip] + y * Rcs[nfp + ip];
      yy = -x * Rcs[nfp + ip] + y * Rcs[ip];
      zz = z; // not needed here, but needed for future speedup

      // Use H & H method to find field contribution due to segment
      xi = xx - x1;
      yi = yy - y1;
      zi = zz - z1;

      xf = xx - x2;
      yf = yy - y2;
      zf = zz - z2;

      ri = sqrt(xi * xi + yi * yi + zi * zi);  
      rf = sqrt(xf * xf + yf * yf + zf * zf);

      l = sqrt((xf-xi) * (xf-xi) + (yf-yi) * (yf-yi) + (zf-zi) * (zf-zi));  

      ex = ((xf-xi)/l);  
      ey = ((yf-yi)/l);  
      ez = ((zf-zi)/l);  

      // Prefactor for cross product 
      coef = cur * l * (ri+rf) / (ri*rf * ((ri+rf) * (ri+rf) - l * l));   

      bx = coef * (ey * zi - ez * yi);
      by = coef * (ez * xi - ex * zi);
      bz = coef * (ex * yi - ey * xi);
  
      bxx += bx * Rcs[ip] - by * Rcs[nfp + ip];
      byy += bx * Rcs[nfp + ip] + by * Rcs[ip];
      bzz += bz;    
   }

   // Write segment contribution to shared memory to prepare for reduction   
   // Swapped double buffer method is performed for each cartesian component
   shdbx[cid * seg + seg] = bxx;
   shdby[cid * seg + seg] = byy;
   shdbz[cid * seg + seg] = bzz;

	__syncthreads();

   // Perform Hillis-Steele reduction and store in global db (x,y,z)   
}

__global__ void hillis_steele(double* dBx, double* dBy, double* dBz, 
								  double* bx,  double* by,  double* bz);

double cosnfp(int ip) { 
//----------------------------------------------------------------------------------------------------
// Cosine of vector component at the ip-th field period
//----------------------------------------------------------------------------------------------------
   return cos((ip - 1) * MY_PI * 2 / Nfp);
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

double sinnfp(int ip) { 
//----------------------------------------------------------------------------------------------------
// Sine of vector component at the ip-th field period
//----------------------------------------------------------------------------------------------------
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
	unsigned int shmem_size = (segs_per_block * 3 + 2 * (segs_per_block - 1) * 3 \
										      + 2 * nfp) * sizeof(double);
   //unsigned int shmem_size = 2 * DimGrid * 3 * sizeof(double); 

   //printf("First point of multifilaments is (%f,%f,%f)\n", mfx[0],mfy[0],mfz[0]);
   //printf("nseg: %d, size_fp: %d\n", nseg, size_fp);
   //printf("First value of current: %f Number of coils: %d Number of fils %d\n", current[0], ncoil, nfil);

	// Create helper currents and reduction buffer unified arrays
	double* Ic = (double*) malloc(ncoil * sizeof(double));
   double* R = (double*) malloc(nfp * 2 * sizeof(double));
	double* dBx;
	double* dBy;
	double* dBz;

	// Fill Ic helper array with values from currents
   for(int i = 0; i < ncoil / nfil; i++) {
      for(int j = 0; j < nfil; j++) {
         Ic[i * nfil + j] = current[i];
      }
   }

   // Fill rotation matrix
   for(int i = 0; i < nfp; i++) {
      R[i] = cosnfp(i + 1); 
      R[nfp + i] = sin(i + 1);
   }
      
   // Define maximum number of segments for kernel memory safety
   int max = nseg * ncoil;

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
	cudaMallocManaged((void**)&dBx, DimGrid * sizeof(double));
	cudaMallocManaged((void**)&dBy, DimGrid * sizeof(double));
	cudaMallocManaged((void**)&dBz, DimGrid * sizeof(double));

   // TODO:DEBUG
   field_kernel<<<DimGrid, dim3(6,129), 10000>>>(mfilx, mfily, mfilz, Ic, 
                                                     xsurf[0], ysurf[0], zsurf[0],
                                                     dBx, dBy, dBz, R, nfp, max);

   // Exit if there is a problem with kernel execution	
   cudaError_t err = cudaGetLastError();
   
   if(err != cudaSuccess) {
      printf("CUDA Error: %s\n", cudaGetErrorString(err));
      exit(-1);
   }

   cudaDeviceSynchronize();

	for(int i = 0; i < size_fp; i++) {
		// Call field kernel to get magnetic field
/*		field_kernel<<<DimGrid, DimBlock, field_shMem_size>>>(mfilx, mfily, mfilz, Ic, 
                                                            xsurf[i], ysurf[i], zsurf[i],
                                                            dBx, dBy, dBz, R, nfp, max);


      cudaDeviceSynchronize(); // synchronize to prepare for reduction
*/		
		//Call Hillis_Steele kernel to scan the partial sums and get bx, by, bz
		//hillis_steele<<<1, DimGrid, scan_shMem_size>>>(dBx, dBy, dBz, bx+i, by+i, bz+i);
	}

	// Call rotation kernel to reflect results to all field periods
	
	// Cleanup
	cudaFree(Ic);
	cudaFree(dBx);
	cudaFree(dBy);
	cudaFree(dBz);
}
