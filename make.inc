
LIB_NAME = modern_fishpack

LIB = ../lib/lib$(LIB_NAME).a

AR = /usr/bin/ar

MAKE = /usr/bin/make

#EXTERNAL_LIBRARY_PATH = /usr/local/lib
EXTERNAL_LIBRARY_PATH = /usr/local/lib64

BIN_PATH = /usr/local/bin

OBJ := $(patsubst %.f90, %.o, $(wildcard *.f90))

OBJ_DIRECTORY = ../objs

LIB_DIRECTORY = ../lib

MOD_DIRECTORY = $(LIB_DIRECTORY)

OBJS = $(addprefix $(OBJ_DIRECTORY)/, $(OBJ))

FC = gfortran

FCFLAGS = \
-O3\
$(ENFORCE_64_BIT_PRECISION) \
$(LANGUAGE_OPTIONS)
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

