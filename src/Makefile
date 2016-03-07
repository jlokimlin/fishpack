
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

#cblktri.o: \
#type_FishpackWorkspace.o \
#comf.o \
#cblktri.f90
#	$(FC) $(FFLAGS) -c cblktri.f90 $(LINKS) \
#	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

gnbnaux.o: \
gnbnaux.f90
	$(FC) $(FFLAGS) -c gnbnaux.f90 $(LINKS) \
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

main.o: \
type_FishpackWorkspace.o \
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