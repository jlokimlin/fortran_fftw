
include ../make.inc

LINKS = \
-I/usr/local/bin/fftw3/api \
-L/usr/local/bin -lfftw3 \
-I$(OBJ_DIRECTORY) \
-I$(LIB_DIRECTORY)

EXECUTABLE = ./build_lib.exe

type_FortranFFTW.o: \
type_FortranFFTW.f90
	$(FC) $(FFLAGS) -c type_FortranFFTW.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

main.o: \
type_FortranFFTW.o \
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