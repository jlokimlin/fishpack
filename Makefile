
include make.inc

EX_DIR = procedural_examples

all: clean lib procedural_examples

lib: 
	mkdir -p ./lib
	mkdir -p ./objs
	( cd ./src; $(MAKE) clean; $(MAKE) all )

procedural_examples:
	( cd ./$(EX_DIR); $(MAKE) clean; $(MAKE) run )

install:
	cp ./lib/lib$(LIB_NAME).a $(EXTERNAL_LIBRARY_PATH)
	cp -r ../$(LIB_NAME) $(BIN_PATH)

clean: 
	( cd ./src; $(MAKE) clean; cd ../$(EX_DIR); $(MAKE) clean )

.PHONY: all lib procedural_examples install
