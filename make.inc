MAKE = /usr/bin/make
AR = /usr/bin/ar
#EXTERNAL_LIBRARY_PATH = /usr/local/lib
EXTERNAL_LIBRARY_PATH = /usr/local/lib64
BIN_PATH = /usr/local/bin

OBJ := $(patsubst %.f90, %.o, $(wildcard *.f90))
OBJ_DIRECTORY = ../objs
OBJS = $(addprefix $(OBJ_DIRECTORY)/, $(OBJ))

LIB_DIRECTORY = ../lib
LIB_NAME = fishpack
LIB = $(LIB_DIRECTORY)/lib$(LIB_NAME).a
MOD_DIRECTORY = $(LIB_DIRECTORY)

FC = gfortran

FFLAGS = \
-O3 \
$(LANGUAGE_OPTIONS) \
$(ERROR_WARNING_OPTIONS)

LANGUAGE_OPTIONS = \
-std=f2008ts \
-fimplicit-none

ERROR_WARNING_OPTIONS = \
-Wall
