
include ../make.inc

LINKS = \
-I$(MOD_DIRECTORY)

fishpack_library.o: \
type_FishpackSolver.o \
type_FishpackWorkspace.o \
type_PoissonSolver.o \
type_TridiagonalSolver.o \
type_CenteredGrid.o \
type_StaggeredGrid.o \
type_MacGrid.o

type_FishpackSolver.o: \
cblktri.o\
blktri.o\
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
sepx4.o

type_TridiagonalSolver.o: \
type_TridiagonalData.o \
type_CenteredGrid.o\
type_StaggeredGrid.o\
genbun.o\
poistg.o

type_TridiagonalData.o: type_Grid.o

type_PoissonSolver.o: type_HelmholtzSolver.o

type_HelmholtzSolver.o: \
hwscrt.o \
hstcrt.o \
type_StaggeredGrid.o\
type_CenteredGrid.o \
type_HelmholtzData.o

type_HelmholtzData.o: type_Grid.o

type_MacGrid.o: type_Grid.o

type_CenteredGrid.o: type_Grid.o

type_StaggeredGrid.o: type_Grid.o 

type_Grid.o: type_RectangularDomain.o

hstcrt.o: genbun.o poistg.o

hstcsp.o: blktri.o

hstcyl.o: genbun.o poistg.o

hstplr.o: genbun.o poistg.o

hstssp.o: genbun.o poistg.o

hwscrt.o: genbun.o poistg.o

hwscsp.o: blktri.o

hwscyl.o: genbun.o

hwsplr.o: genbun.o

hwsssp.o: genbun.o

sepeli.o: blktri.o sepaux.o

sepx4.o: genbun.o sepaux.o

blktri.o: type_FishpackWorkspace.o comf.o

cblktri.o: type_FishpackWorkspace.o comf.o

genbun.o: type_FishpackWorkspace.o gnbnaux.o

poistg.o: type_FishpackWorkspace.o gnbnaux.o

hw3crt.o: pois3d.o 

pois3d.o: type_FishpackWorkspace.o type_FFTpack.o

cmgnbn.o: type_FishpackWorkspace.o 

$(LIB): $(OBJS)
	$(AR) -rv $@ $?

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< $(LINKS) -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)
	 
clean:
	rm -f $(OBJS) $(LIB) $(MOD_DIRECTORY)/*.mod
