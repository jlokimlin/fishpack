
include ../make.inc

LINKS = \
-I$(OBJ_DIRECTORY) \
-I$(LIB_DIRECTORY)

EXECUTABLE = ./build_lib.exe

type_FishpackWorkspace.o: \
type_FishpackWorkspace.f90
	$(FC) $(FFLAGS) -c type_FishpackWorkspace.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

comf.o: \
comf.f90
	$(FC) $(FFLAGS) -c comf.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

blktri.o: \
type_FishpackWorkspace.o \
comf.o \
blktri.f90
	$(FC) $(FFLAGS) -c blktri.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

cblktri.o: \
type_FishpackWorkspace.o \
comf.o \
cblktri.f90
	$(FC) $(FFLAGS) -c cblktri.f90 $(LINKS) \
	-o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

cmgnbn.o: \
type_FishpackWorkspace.o \
cmgnbn.f90
	$(FC) $(FFLAGS) -c cmgnbn.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

gnbnaux.o: \
gnbnaux.f90
	$(FC) $(FFLAGS) -c gnbnaux.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_FFTpack.o: \
type_FFTpack.f90
	$(FC) $(FFLAGS) -c type_FFTpack.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)
	 
genbun.o: \
type_FishpackWorkspace.o \
gnbnaux.o \
genbun.f90
	$(FC) $(FFLAGS) -c genbun.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

poistg.o: \
type_FishpackWorkspace.o \
gnbnaux.o \
poistg.f90
	$(FC) $(FFLAGS) -c poistg.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

pois3d.o: \
type_FishpackWorkspace.o \
type_FFTpack.o \
pois3d.f90
	$(FC) $(FFLAGS) -c pois3d.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

hw3crt.o: \
type_FishpackWorkspace.o \
pois3d.o \
hw3crt.f90
	$(FC) $(FFLAGS) -c hw3crt.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)
	 
hstcrt.o: \
type_FishpackWorkspace.o \
genbun.o \
poistg.o \
hstcrt.f90
	$(FC) $(FFLAGS) -c hstcrt.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

hstcsp.o: \
type_FishpackWorkspace.o \
blktri.o\
hstcsp.f90
	$(FC) $(FFLAGS) -c hstcsp.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

hstcyl.o: \
type_FishpackWorkspace.o \
genbun.o \
poistg.o \
hstcrt.f90
	$(FC) $(FFLAGS) -c hstcyl.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

hstplr.o: \
type_FishpackWorkspace.o \
genbun.o \
poistg.o \
hstcrt.f90
	$(FC) $(FFLAGS) -c hstplr.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

hstssp.o: \
type_FishpackWorkspace.o \
genbun.o \
poistg.o \
hstssp.f90
	$(FC) $(FFLAGS) -c hstssp.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

hwscrt.o: \
type_FishpackWorkspace.o \
genbun.o \
poistg.o \
hwscrt.f90
	$(FC) $(FFLAGS) -c hwscrt.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

hwscsp.o: \
type_FishpackWorkspace.o \
blktri.o \
hwscsp.f90
	$(FC) $(FFLAGS) -c hwscsp.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

hwscyl.o: \
type_FishpackWorkspace.o \
genbun.o \
hwscyl.f90
	$(FC) $(FFLAGS) -c hwscyl.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

hwsplr.o: \
type_FishpackWorkspace.o \
genbun.o \
hwsplr.f90
	$(FC) $(FFLAGS) -c hwsplr.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

hwsssp.o: \
type_FishpackWorkspace.o \
genbun.o \
hwsssp.f90
	$(FC) $(FFLAGS) -c hwsssp.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

sepaux.o: \
sepaux.f90
	$(FC) $(FFLAGS) -c sepaux.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

sepeli.o: \
type_FishpackWorkspace.o\
blktri.o \
sepaux.o \
sepeli.f90
	$(FC) $(FFLAGS) -c sepeli.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

sepx4.o: \
type_FishpackWorkspace.o\
genbun.o \
sepaux.o \
sepx4.f90
	$(FC) $(FFLAGS) -c sepx4.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_RectangularDomain.o: \
type_RectangularDomain.f90
	$(FC) $(FFLAGS) -c type_RectangularDomain.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_Grid.o: \
type_RectangularDomain.o \
type_Grid.f90
	$(FC) $(FFLAGS) -c type_Grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_CenteredGrid.o: \
type_Grid.o \
type_CenteredGrid.f90
	$(FC) $(FFLAGS) -c type_CenteredGrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_StaggeredGrid.o: \
type_Grid.o \
type_StaggeredGrid.f90
	$(FC) $(FFLAGS) -c type_StaggeredGrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_MacGrid.o:\
type_Grid.o \
type_MacGrid.f90
	$(FC) $(FFLAGS) -c type_MacGrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_HelmholtzData.o: \
type_RectangularDomain.o \
type_Grid.o\
type_HelmholtzData.f90
	$(FC) $(FFLAGS) -c type_HelmholtzData.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_HelmholtzSolver.o: \
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
type_HelmholtzSolver.o \
type_HelmholtzSolver.f90
	$(FC) $(FFLAGS) -c type_PoissonSolver.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_TridiagonalData.o: \
type_Grid.o\
type_TridiagonalData.f90
	$(FC) $(FFLAGS) -c type_TridiagonalData.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_TridiagonalSolver.o: \
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
hwsplr.o\
hwsssp.o\
pois3d.o\
poistg.o\
sepeli.o\
sepx4.o\
type_FishpackSolver.f90
	$(FC) $(FFLAGS) -c type_FishpackSolver.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)
	 
fishpack_library.o: \
type_FishpackWorkspace.o\
type_FishpackSolver.o\
type_HelmholtzData.o\
type_HelmholtzSolver.o\
type_PoissonSolver.o \
type_TridiagonalData.o\
type_TridiagonalSolver.o\
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