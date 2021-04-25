# at PPPL, these modules are needed
# module load gcc szip hdf hdf5-serial curl netcdf-c

PROGRAM = OMIC

FILES.c = startup.c read_namelist.c omic.c read_focus.c single_fil.c multi_fil.c bfield.c alpha.c output.c solvers.c 

FILES.h = startup.h read_namelist.h read_focus.h single_fil.h multi_fil.h bfield.h globals.h alpha.h output.h solvers.h

FILES.o = ${FILES.c:.c=.o} bfield_gpu.o

CC      = mpicc

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

NETCDF_HOME = ${NETCDF_C_HOME}

NETCDF = -I ${NETCDF_HOME}/include -L ${NETCDF_HOME}/lib -lnetcdf

CFLAGS  = ${SFLAGS} ${GFLAGS} ${OFLAGS} ${WFLAGS} ${UFLAGS} -fopenmp 

CUDAFLAGS = -lcudart

LDFLAGS = 

LDLIBS  =

all:    ${PROGRAM}

${PROGRAM}: ${FILES.o}

	${CC} -o $@ ${CFLAGS} ${LDFLAGS} ${CUDAFLAGS} ${LDLIBS} $^ ${NETCDF} -lm

read_namelist.o: read_namelist.c globals.h read_namelist.h
	${CC} ${NETCDF} -c $< -o $@
omic.o: omic.c read_namelist.h 
	${CC}  -c $< -o $@
read_focus.o: read_focus.c read_focus.h globals.h
	${CC} ${NETCDF} -c $< -o $@
single_fil.o: single_fil.c single_fil.h globals.h
	${CC} -fopenmp ${NETCDF} -c $< -o $@
multi_fil.o: multi_fil.c multi_fil.h globals.h
	${CC} -fopenmp ${NETCDF} -c $< -o $@
bfield_gpu.o: bfield_gpu.cu bfield_gpu.h
	${NVCC} -lcudart ${NETCDF} -Xcompiler -O3 -Xcompiler -Wall -Xptxas -O3 -c $< -o $@
bfield.o: bfield.c bfield.h  bfield_gpu.h globals.h
	${NVCC} -lcudart ${NETCDF} -Xcompiler -O3 -Xcompiler -Wall -Xptxas -O3 -c $< -o $@
alpha.o: alpha.c alpha.h globals.h
	${CC} -c $< -o $@
#output.o: output.c output.h globals.h
#	${CC} -c $< -o $@
solvers.o: solvers.c solvers.h globals.h
	${CC} -c $< -o $@
startup.o: startup.c startup.h globals.h
	${CC} -c $< -o $@

DEBRIS = a.out core *~ *.dSYM

RM_FR  = rm -fr

clean:
	${RM_FR} ${FILES.o} ${PROGRAM} ${DEBRIS}
