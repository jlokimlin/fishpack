
include make.inc

EX_DIR = examples

all: clean lib examples

lib: 
	mkdir -p ./lib
	mkdir -p ./objs
	( cd ./src; $(MAKE) clean; $(MAKE) all )

examples:
	( cd ./$(EX_DIR); $(MAKE) clean; $(MAKE) run )

install:
	cp ./lib/lib$(LIB_NAME).a $(EXTERNAL_LIBRARY_PATH)
	cp -r ../$(LIB_NAME) $(BIN_PATH)

clean: 
	( cd ./src; $(MAKE) clean; cd ../$(EX_DIR); $(MAKE) clean )

.PHONY: all lib examples install
