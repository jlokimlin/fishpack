
include ../make.inc

LINKS = \
-I$(OBJ_DIRECTORY) \
-I$(LIB_DIRECTORY)

EXECUTABLE = ./build_lib.exe

fishpack_precision.o: \
fishpack_precision.f90
	$(FC) $(FFLAGS) -c fishpack_precision.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_FishpackWorkspace.o: \
fishpack_precision.o \
type_FishpackWorkspace.f90
	$(FC) $(FFLAGS) -c type_FishpackWorkspace.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_GeneralizedCyclicReductionUtility.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
type_GeneralizedCyclicReductionUtility.f90
	$(FC) $(FFLAGS) -c type_GeneralizedCyclicReductionUtility.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

real_block_tridiagonal_linear_systems_solver.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
type_GeneralizedCyclicReductionUtility.o \
real_block_tridiagonal_linear_systems_solver.f90
	$(FC) $(FFLAGS) -c real_block_tridiagonal_linear_systems_solver.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

complex_block_tridiagonal_linear_systems_solver.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
type_GeneralizedCyclicReductionUtility.o \
complex_block_tridiagonal_linear_systems_solver.f90
	$(FC) $(FFLAGS) -c complex_block_tridiagonal_linear_systems_solver.f90 $(LINKS) \
	-o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

complex_linear_systems_solver.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
complex_linear_systems_solver.f90
	$(FC) $(FFLAGS) -c complex_linear_systems_solver.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_CyclicReductionUtility.o: \
fishpack_precision.o \
type_CyclicReductionUtility.f90
	$(FC) $(FFLAGS) -c type_CyclicReductionUtility.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_PeriodicFastFourierTransform.o: \
fishpack_precision.o \
type_PeriodicFastFourierTransform.f90
	$(FC) $(FFLAGS) -c type_PeriodicFastFourierTransform.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_CenteredCyclicReductionUtility.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
type_CyclicReductionUtility.o \
type_CenteredCyclicReductionUtility.f90
	$(FC) $(FFLAGS) -c type_CenteredCyclicReductionUtility.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_StaggeredCyclicReductionUtility.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
type_CyclicReductionUtility.o \
type_StaggeredCyclicReductionUtility.f90
	$(FC) $(FFLAGS) -c type_StaggeredCyclicReductionUtility.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

three_dimensional_solvers.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
type_PeriodicFastFourierTransform.o \
three_dimensional_solvers.f90
	$(FC) $(FFLAGS) -c three_dimensional_solvers.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

general_linear_systems_solver_3d.o: \
three_dimensional_solvers.o \
general_linear_systems_solver_3d.f90
	$(FC) $(FFLAGS) -c general_linear_systems_solver_3d.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

centered_cartesian_helmholtz_solver_3d.o: \
three_dimensional_solvers.o \
centered_cartesian_helmholtz_solver_3d.f90
	$(FC) $(FFLAGS) -c centered_cartesian_helmholtz_solver_3d.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

staggered_helmholtz_solvers.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
type_CenteredCyclicReductionUtility.o \
type_StaggeredCyclicReductionUtility.o \
type_GeneralizedCyclicReductionUtility.o \
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

staggered_real_linear_systems_solver.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
type_StaggeredCyclicReductionUtility.o \
staggered_real_linear_systems_solver.f90
	$(FC) $(FFLAGS) -c staggered_real_linear_systems_solver.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)
	 
centered_helmholtz_solvers.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
type_CenteredCyclicReductionUtility.o \
type_GeneralizedCyclicReductionUtility.o \
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

centered_real_linear_systems_solver.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
type_CenteredCyclicReductionUtility.o \
centered_real_linear_systems_solver.f90
	$(FC) $(FFLAGS) -c centered_real_linear_systems_solver.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_SepAux.o: \
fishpack_precision.o \
type_SepAux.f90
	$(FC) $(FFLAGS) -c type_SepAux.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

sepeli.o: \
fishpack_precision.o \
type_FishpackWorkspace.o\
type_GeneralizedCyclicReductionUtility.o \
type_SepAux.o \
sepeli.f90
	$(FC) $(FFLAGS) -c sepeli.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

sepx4.o: \
fishpack_precision.o \
type_FishpackWorkspace.o\
centered_real_linear_systems_solver.o \
type_SepAux.o \
sepx4.f90
	$(FC) $(FFLAGS) -c sepx4.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

fishpack.o: \
fishpack_precision.o \
real_block_tridiagonal_linear_systems_solver.o \
complex_block_tridiagonal_linear_systems_solver.o\
complex_linear_systems_solver.o\
centered_real_linear_systems_solver.o\
three_dimensional_solvers.o\
staggered_real_linear_systems_solver.o\
sepeli.o\
sepx4.o\
centered_helmholtz_solvers.o \
staggered_helmholtz_solvers.o \
type_FishpackWorkspace.o\
type_PeriodicFastFourierTransform.o\
fishpack.f90
	$(FC) $(FFLAGS) -c fishpack.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

main.o: \
fishpack.o \
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
