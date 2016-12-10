
include ../make.inc

LINKS = \
-L$(LIB_DIRECTORY) -l$(LIB_NAME) \
-I$(MOD_DIRECTORY)

SRC=$(subst .o,.f90, $(OBJ))

EXES=$(subst .o,.exe, $(OBJ))

run: 
	make clean;
	make $(EXES)

$(EXES): $(SRC)

%.exe : %.f90
	rm -f $@
	$(FC) $(FFLAGS) $< -o $@ $(LINKS)
	./$@

clean:
	rm -f $(OBJS) $(EXES)
