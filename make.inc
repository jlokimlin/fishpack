
LIB = ../lib/libmodern_fishpack.a

LIB_NAME = modern_fishpack

AR = /usr/bin/ar

MAKE = /usr/bin/make

#EXTERNAL_LIBRARY_PATH = /usr/local/lib
EXTERNAL_LIBRARY_PATH = /usr/local/lib64

BIN_PATH = /usr/local/bin

OBJ := $(patsubst %.f90, %.o, $(wildcard *.f90))

OBJ_DIRECTORY = ../objs

LIB_DIRECTORY = ../lib

OBJS = $(addprefix $(OBJ_DIRECTORY)/, $(OBJ))

FC = gfortran

FFLAGS = \
-O3 \
$(ENFORCE_64_BIT_PRECISION) \
$(LANGUAGE_OPTIONS)\
$(ERROR_WARNING_OPTIONS)

ENFORCE_64_BIT_PRECISION = \
-fdefault-real-8 \
-fdefault-double-8

LANGUAGE_OPTIONS = \
-fimplicit-none

ERROR_WARNING_OPTIONS = \
-Wall

CODE_GENERATION_OPTIONS = \
-fcheck=all
