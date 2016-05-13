
include make.inc

all: clean lib testlib

lib: 
	mkdir -p ./lib
	mkdir -p ./objs
	( cd ./src; $(MAKE) clean; $(MAKE) all )

testlib:
	( cd ./test; $(MAKE) clean; $(MAKE) run )

install:
	cp ./lib/lib$(LIB_NAME).a $(EXTERNAL_LIBRARY_PATH)
	cp -r ../$(LIB_NAME) $(BIN_PATH)

clean: 
	( cd ./src; $(MAKE) clean; cd ../test; $(MAKE) clean )

.PHONY: all lib testlib install