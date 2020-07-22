# at PPPL, these modules are needed
# module load gcc szip hdf hdf5-serial curl netcdf-c

PROGRAM = multi

FILES.c = startup.c read_namelist.c multi.c read_focus.c single_fil.c multi_fil.c bfield.c alpha.c output.c solvers.c

FILES.h = libfyaml.h config.h startup.h read_namelist.h read_focus.h single_fil.h multi_fil.h bfield.h globals.h alpha.h output.h solvers.h

FILES.o = ${FILES.c:.c=.o}

CC      = mpicc

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

CFLAGS  = ${SFLAGS} ${GFLAGS} ${OFLAGS} ${WFLAGS} ${UFLAGS} -fopenmp #'pkg-config --cflag libfyaml'

LDFLAGS = #'pkg-config --libs libfyaml'

LDLIBS  =

all:    ${PROGRAM}

${PROGRAM}: ${FILES.o}

	${CC} -o $@ ${CFLAGS} ${LDFLAGS} ${LDLIBS} $^ ${NETCDF} -lm

read_namelist.o: read_namelist.c globals.h read_namelist.h
	${CC} ${NETCDF} -c $< -o $@
multi.o: multi.c read_namelist.h 
	${CC}  -c $< -o $@
read_focus.o: read_focus.c read_focus.h globals.h
	${CC} ${NETCDF} -c $< -o $@
single_fil.o: single_fil.c single_fil.h globals.h
	${CC} -fopenmp ${NETCDF} -c $< -o $@
multi_fil.o: multi_fil.c multi_fil.h globals.h
	${CC} -fopenmp ${NETCDF} -c $< -o $@
bfield.o: bfield.c bfield.h globals.h
	${CC} ${NETCDF} -c $< -o $@
alpha.o: alpha.c alpha.h globals.h
	${CC} -c $< -o $@
#output.o: output.c output.h globals.h
#	${CC} -c $< -o $@
solvers.o: solvers.c solvers.h globals.h
	${CC} -c $< -o $@
startup.o: startup.c startup.h globals.h config.h
	${CC} -c $< -o $@

DEBRIS = a.out core *~ *.dSYM

RM_FR  = rm -fr



clean:
	${RM_FR} ${FILES.o} ${PROGRAM} ${DEBRIS}
