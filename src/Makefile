
include ../make.inc

LINKS = \
-I$(OBJ_DIRECTORY) \
-I$(LIB_DIRECTORY)

EXECUTABLE = ./build_lib.exe

fishpack_precision.o: \
fishpack_precision.f90
	$(FC) $(FFLAGS) -c fishpack_precision.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

fftpack_routines.o: \
fishpack_precision.o \
fftpack_routines.f90
	$(FC) $(FFLAGS) -c fftpack_routines.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_FishpackWorkspace.o: \
fishpack_precision.o \
fftpack_routines.o \
type_FishpackWorkspace.f90
	$(FC) $(FFLAGS) -c type_FishpackWorkspace.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_ComfAux.o: \
fishpack_precision.o \
type_ComfAux.f90
	$(FC) $(FFLAGS) -c type_ComfAux.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_BlktriAux.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
type_ComfAux.o \
type_BlktriAux.f90
	$(FC) $(FFLAGS) -c type_BlktriAux.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

real_block_tridiagonal_linear_systems_solver.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
type_BlktriAux.o \
real_block_tridiagonal_linear_systems_solver.f90
	$(FC) $(FFLAGS) -c real_block_tridiagonal_linear_systems_solver.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

complex_block_tridiagonal_linear_systems_solver.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
type_ComfAux.o \
complex_block_tridiagonal_linear_systems_solver.f90
	$(FC) $(FFLAGS) -c complex_block_tridiagonal_linear_systems_solver.f90 $(LINKS) \
	-o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

cmgnbn.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
cmgnbn.f90
	$(FC) $(FFLAGS) -c cmgnbn.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_GenbunAux.o: \
fishpack_precision.o \
type_GenbunAux.f90
	$(FC) $(FFLAGS) -c type_GenbunAux.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_FFTpack.o: \
fishpack_precision.o \
type_FFTpack.f90
	$(FC) $(FFLAGS) -c type_FFTpack.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

genbun.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
type_GenbunAux.o \
genbun.f90
	$(FC) $(FFLAGS) -c genbun.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

poistg.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
type_GenbunAux.o \
poistg.f90
	$(FC) $(FFLAGS) -c poistg.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

pois3d.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
type_FFTpack.o \
pois3d.f90
	$(FC) $(FFLAGS) -c pois3d.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

hw3crt.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
pois3d.o \
hw3crt.f90
	$(FC) $(FFLAGS) -c hw3crt.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

staggered_helmholtz_solvers.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
genbun.o \
poistg.o \
type_BlktriAux.o \
staggered_helmholtz_solvers.f90
	$(FC) $(FFLAGS) -c staggered_helmholtz_solvers.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)
	 
staggered_cartesian_solver.o: \
staggered_helmholtz_solvers.o \
staggered_cartesian_solver.f90
	$(FC) $(FFLAGS) -c staggered_cartesian_solver.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

staggered_axisymmetric_spherical_solver.o: \
staggered_helmholtz_solvers.o \
staggered_axisymmetric_spherical_solver.f90
	$(FC) $(FFLAGS) -c staggered_axisymmetric_spherical_solver.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

staggered_cylindrical_solver.o: \
staggered_helmholtz_solvers.o \
staggered_cylindrical_solver.f90
	$(FC) $(FFLAGS) -c staggered_cylindrical_solver.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

staggered_polar_solver.o: \
staggered_helmholtz_solvers.o \
staggered_polar_solver.f90
	$(FC) $(FFLAGS) -c staggered_polar_solver.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

staggered_spherical_solver.o: \
staggered_helmholtz_solvers.o \
staggered_spherical_solver.f90
	$(FC) $(FFLAGS) -c staggered_spherical_solver.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

centered_helmholtz_solvers.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
genbun.o \
poistg.o \
type_BlktriAux.o \
centered_helmholtz_solvers.f90
	$(FC) $(FFLAGS) -c centered_helmholtz_solvers.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

centered_cartesian_solver.o: \
centered_helmholtz_solvers.o \
centered_cartesian_solver.f90
	$(FC) $(FFLAGS) -c centered_cartesian_solver.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

centered_axisymmetric_spherical_solver.o: \
centered_helmholtz_solvers.o \
centered_axisymmetric_spherical_solver.f90
	$(FC) $(FFLAGS) -c centered_axisymmetric_spherical_solver.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

centered_cylindrical_solver.o: \
centered_helmholtz_solvers.o \
centered_cylindrical_solver.f90
	$(FC) $(FFLAGS) -c centered_cylindrical_solver.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

centered_polar_solver.o: \
centered_helmholtz_solvers.o \
centered_polar_solver.f90
	$(FC) $(FFLAGS) -c centered_polar_solver.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

centered_spherical_solver.o: \
centered_helmholtz_solvers.o \
centered_spherical_solver.f90
	$(FC) $(FFLAGS) -c centered_spherical_solver.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_SepAux.o: \
fishpack_precision.o \
type_SepAux.f90
	$(FC) $(FFLAGS) -c type_SepAux.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

sepeli.o: \
fishpack_precision.o \
type_FishpackWorkspace.o\
type_BlktriAux.o \
type_SepAux.o \
sepeli.f90
	$(FC) $(FFLAGS) -c sepeli.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

sepx4.o: \
fishpack_precision.o \
type_FishpackWorkspace.o\
genbun.o \
type_SepAux.o \
sepx4.f90
	$(FC) $(FFLAGS) -c sepx4.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_RectangularDomain.o: \
fishpack_precision.o \
type_RectangularDomain.f90
	$(FC) $(FFLAGS) -c type_RectangularDomain.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_FishpackGrid.o: \
fishpack_precision.o \
type_FishpackGrid.f90
	$(FC) $(FFLAGS) -c type_FishpackGrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_Grid.o: \
fishpack_precision.o \
type_RectangularDomain.o \
type_Grid.f90
	$(FC) $(FFLAGS) -c type_Grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_CenteredGrid.o: \
fishpack_precision.o \
type_Grid.o \
type_CenteredGrid.f90
	$(FC) $(FFLAGS) -c type_CenteredGrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_StaggeredGrid.o: \
fishpack_precision.o \
type_Grid.o \
type_StaggeredGrid.f90
	$(FC) $(FFLAGS) -c type_StaggeredGrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_MacGrid.o:\
fishpack_precision.o \
type_Grid.o \
type_MacGrid.f90
	$(FC) $(FFLAGS) -c type_MacGrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_HelmholtzData.o: \
fishpack_precision.o \
type_RectangularDomain.o \
type_Grid.o\
type_HelmholtzData.f90
	$(FC) $(FFLAGS) -c type_HelmholtzData.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_HelmholtzSolver.o: \
fishpack_precision.o \
staggered_helmholtz_solvers.o \
centered_helmholtz_solvers.o \
type_HelmholtzData.o \
type_Grid.o\
type_StaggeredGrid.o\
type_CenteredGrid.o\
type_HelmholtzSolver.f90
	$(FC) $(FFLAGS) -c type_HelmholtzSolver.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_PoissonSolver.o: \
fishpack_precision.o \
type_HelmholtzSolver.o \
type_HelmholtzSolver.f90
	$(FC) $(FFLAGS) -c type_PoissonSolver.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_TridiagonalData.o: \
fishpack_precision.o \
type_Grid.o\
type_TridiagonalData.f90
	$(FC) $(FFLAGS) -c type_TridiagonalData.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_TridiagonalSolver.o: \
fishpack_precision.o \
type_TridiagonalData.o \
type_Grid.o\
type_CenteredGrid.o\
type_StaggeredGrid.o\
genbun.o\
poistg.o\
type_TridiagonalSolver.f90
	$(FC) $(FFLAGS) -c type_TridiagonalSolver.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_FishpackSolver.o: \
fishpack_precision.o \
real_block_tridiagonal_linear_systems_solver.o \
complex_block_tridiagonal_linear_systems_solver.o\
cmgnbn.o\
genbun.o\
pois3d.o\
hw3crt.o\
poistg.o\
sepeli.o\
sepx4.o\
centered_helmholtz_solvers.o \
staggered_helmholtz_solvers.o \
type_FishpackSolver.f90
	$(FC) $(FFLAGS) -c type_FishpackSolver.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

fishpack_library.o: \
fishpack_precision.o \
real_block_tridiagonal_linear_systems_solver.o \
complex_block_tridiagonal_linear_systems_solver.o\
cmgnbn.o\
genbun.o\
pois3d.o\
hw3crt.o\
poistg.o\
sepeli.o\
sepx4.o\
centered_helmholtz_solvers.o \
staggered_helmholtz_solvers.o \
type_FishpackWorkspace.o\
type_FishpackSolver.o\
type_HelmholtzData.o\
type_HelmholtzSolver.o\
type_PoissonSolver.o \
type_TridiagonalData.o\
type_TridiagonalSolver.o\
type_FishpackGrid.o \
type_CenteredGrid.o\
type_StaggeredGrid.o\
type_MacGrid.o\
fishpack_library.f90
	$(FC) $(FFLAGS) -c fishpack_library.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

main.o: \
fishpack_library.o \
main.f90
	$(FC) $(FFLAGS) -c main.f90 $(LINKS) \
	-o $(OBJ_DIRECTORY)/$@

exec: $(OBJ)
	$(FC) $(FFLAGS) -o $(EXECUTABLE) $(OBJS) $(LINKS)

lib: 
	$(AR) rc lib$(LIB_NAME).a $(OBJ_DIRECTORY)/*.o
	mv lib$(LIB_NAME).a $(LIB_DIRECTORY)

all:
	make clean
	make exec; $(EXECUTABLE)
	make lib

clean:
	rm -f $(OBJS) $(LIB) $(LIB_DIRECTORY)/*.mod $(LIB_DIRECTORY)/*.smod
	rm -f $(EXECUTABLE)
