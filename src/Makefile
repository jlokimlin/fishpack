
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

type_ComfAux.o: \
fishpack_precision.o \
type_ComfAux.f90
	$(FC) $(FFLAGS) -c type_ComfAux.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

blktri.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
type_ComfAux.o \
blktri.f90
	$(FC) $(FFLAGS) -c blktri.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

cblktri.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
type_ComfAux.o \
cblktri.f90
	$(FC) $(FFLAGS) -c cblktri.f90 $(LINKS) \
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

helmholtz_solvers.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
genbun.o \
poistg.o \
helmholtz_solvers.f90
	$(FC) $(FFLAGS) -c helmholtz_solvers.f90 $(LINKS) \
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

hstcrt.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
genbun.o \
poistg.o \
hstcrt.f90
	$(FC) $(FFLAGS) -c hstcrt.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

hstcsp.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
blktri.o\
hstcsp.f90
	$(FC) $(FFLAGS) -c hstcsp.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

hstcyl.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
genbun.o \
poistg.o \
hstcrt.f90
	$(FC) $(FFLAGS) -c hstcyl.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

hstplr.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
genbun.o \
poistg.o \
hstcrt.f90
	$(FC) $(FFLAGS) -c hstplr.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

hstssp.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
genbun.o \
poistg.o \
hstssp.f90
	$(FC) $(FFLAGS) -c hstssp.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

hwscrt.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
genbun.o \
poistg.o \
hwscrt.f90
	$(FC) $(FFLAGS) -c hwscrt.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

hwscsp.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
blktri.o \
hwscsp.f90
	$(FC) $(FFLAGS) -c hwscsp.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

hwscyl.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
genbun.o \
hwscyl.f90
	$(FC) $(FFLAGS) -c hwscyl.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

centered_polar_solver.o: \
helmholtz_solvers.o \
centered_polar_solver.f90
	$(FC) $(FFLAGS) -c centered_polar_solver.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

hwsssp.o: \
fishpack_precision.o \
type_FishpackWorkspace.o \
genbun.o \
hwsssp.f90
	$(FC) $(FFLAGS) -c hwsssp.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_SepAux.o: \
fishpack_precision.o \
type_SepAux.f90
	$(FC) $(FFLAGS) -c type_SepAux.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

sepeli.o: \
fishpack_precision.o \
type_FishpackWorkspace.o\
blktri.o \
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
hwscrt.o \
hstcrt.o \
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
blktri.o\
cblktri.o\
cmgnbn.o\
genbun.o\
hstcrt.o\
hstcsp.o\
hstcyl.o\
hstplr.o\
hstssp.o\
hw3crt.o\
hwscrt.o\
hwscsp.o\
hwscyl.o\
helmholtz_solvers.o\
hwsssp.o\
pois3d.o\
poistg.o\
sepeli.o\
sepx4.o\
type_FishpackSolver.f90
	$(FC) $(FFLAGS) -c type_FishpackSolver.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

fishpack_library.o: \
fishpack_precision.o \
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
	rm -f $(OBJS) $(LIB) $(LIB_DIRECTORY)/*.mod
	rm -f $(EXECUTABLE)
