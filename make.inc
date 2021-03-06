
INCLUDES = -I. -I/usr/local/include

LIB = ../lib/libfortran_fftw.a

LIB_NAME = fortran_fftw

AR = /usr/bin/ar

MAKE = /usr/bin/make

# Set directories

#EXTERNAL_LIBRARY_PATH = /usr/local/lib
EXTERNAL_LIBRARY_PATH = /usr/local/lib64

BIN_PATH = /usr/local/bin

OBJ := $(patsubst %.f90, %.o, $(wildcard *.f90))

OBJ_DIRECTORY = ../objs

LIB_DIRECTORY = ../lib

OBJS = $(addprefix $(OBJ_DIRECTORY)/, $(OBJ))

# Set compiler options

FC = gfortran

FFLAGS = \
-O3 \
$(LANGUAGE_OPTIONS) \
$(ERROR_WARNING_OPTIONS) \
$(INCLUDES)

LANGUAGE_OPTIONS = \
-fimplicit-none

ERROR_WARNING_OPTIONS = \
-Wall

CODE_GENERATION_OPTIONS = \
-fcheck=all
