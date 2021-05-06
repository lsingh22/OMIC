# at PPPL, these modules are needed
# module load gcc szip hdf hdf5-serial curl netcdf-c

PROGRAM = OMIC

FILES.c = startup.c read_namelist.c omic.c read_focus.c single_fil.c multi_fil.c bfield.c alpha.c output.c solvers.c 

FILES.h = startup.h read_namelist.h read_focus.h single_fil.h multi_fil.h bfield.h globals.h alpha.h output.h solvers.h

FILES.o = ${FILES.c:.c=.o} bfield_gpu.o

CC      = mpicc

#CC		  = tau_cc

NVCC	  = nvcc

SFLAGS  = -std=c11

GFLAGS  = -g

OFLAGS  = -O3

WFLAG1  = -Wall

WFLAG2  = -Wextra

WFLAG3  = -Werror

WFLAG4  = -Wstrict-prototypes

WFLAG5  = -Wmissing-prototypes

WFLAGS  = ${WFLAG1} ${WFLAG2} #${WFLAG3} ${WFLAG4} ${WFLAG5}

UFLAGS  = # Set on command line only

#For use on PPPL cluster
#NETCDF_HOME = ${NETCDF_C_HOME}
#NETCDF = -I ${NETCDF_HOME}/include -L ${NETCDF_HOME}/lib -lnetcdf

NETCDF = -I /usr/local/netcdf/x86_64/4.8.0/include

NETCDF_LIB = -L /usr/local/netcdf/x86_64/4.8.0/lib64 -lnetcdf

#NETCDF = /usr/local/netcdf/x86_64/4.8.0 -L/lib -lnetcdf

#DEBUG	  = -g

#CUDADEBUG = -g -lineinfo

CFLAGS  = ${DEBUG} ${SFLAGS} ${OFLAGS} ${WFLAGS} ${UFLAGS} ${NETCDF} -fopenmp 

CUDAFLAGS = -lcudadevrt -lcudart

LDFLAGS = ${CUDAFLAGS} ${NETCDF_LIB} -lgomp

LDLIBS  =

#GENCODE_SM70 := -gencode arch=compute_70,code=sm_70 # CC 7 on Travis

#GENCODE_SM80 := -gencode arch=compute_80,code=sm_80 # CC 8 on Euler

all:    ${PROGRAM} 

${PROGRAM}: ${FILES.o}

	${CC} -o $@ ${LDFLAGS} $^ ${NETCDF} -lm 

read_namelist.o: read_namelist.c globals.h read_namelist.h
	${CC} ${CFLAGS} -c $< -o $@
omic.o: omic.c read_namelist.h 
	${CC} ${CFLAGS} -c $< -o $@
read_focus.o: read_focus.c read_focus.h globals.h
	${CC} ${CFLAGS} -c $< -o $@
single_fil.o: single_fil.c single_fil.h globals.h
	${CC} ${CFLAGS} -c $< -o $@
bfield_gpu.o: bfield_gpu.cu bfield_gpu.cuh globals.h 
	${NVCC} ${CUDADEBUG} -Xcompiler -O3 -Xcompiler -Wall -Xptxas -O3 -c $< -o $@
bfield.o: bfield.c bfield.h bfield_gpu.cuh globals.h
	${CC} ${CFLAGS} -c $< -o $@ 
multi_fil.o: multi_fil.c multi_fil.h bfield_gpu.cuh globals.h
	${CC} ${CFLAGS} -c $< -o $@
alpha.o: alpha.c alpha.h globals.h
	${CC} ${CFLAGS} -c $< -o $@
solvers.o: solvers.c solvers.h globals.h
	${CC} ${CFLAGS} -c $< -o $@
startup.o: startup.c startup.h globals.h
	${CC} ${CFLAGS} -c $< -o $@

DEBRIS = a.out core *~ *.dSYM

RM_FR  = rm -fr

clean:
	${RM_FR} ${FILES.o} ${PROGRAM} ${DEBRIS}
