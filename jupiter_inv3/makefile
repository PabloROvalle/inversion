FC= /usr/local/intel/bin/intel64/ifort
FCFLAGS = -traceback -ftrapuv -fp-stack-check -check all

F_SOURCES = ludcmp.f lubksb.f matinv.f cont_h2he.f
F90_SOURCES = info_content.f90 convolution.f90 read_file.f90 voigt.f90 gravity.f90

OBJ= $(F_SOURCES:.f=.o)
OBJ90= $(F90_SOURCES:.f90=.o)
LIBRARY = /obs/tfouchet/lib/jwst/jwst.a
LIBRARY90 = /obs/tfouchet/lib/jwst/jwst90.a
LIBCODE= $(LIBRARY)($(OBJ))
LIBCODE90= $(LIBRARY90)($(OBJ90))

# target definitions

bin: jupiter

jupiter: /obs/tfouchet/bin/jwst/jupiter

/obs/tfouchet/bin/jwst/jupiter: jupiter.f90 $(LIBRARY90) $(LIBRARY)  
	$(FC) $(FCFLAGS) $^ -o $@ 

lib $(LIBRARY): $(LIBCODE)

lib90 $(LIBRARY90): $(LIBCODE90)

$(LIBRARY)(%.o): %.f
	$(FC) $(FCFLAGS) -c $<
	ar rv $(LIBRARY) $*.o
	ranlib $(LIBRARY)
	-rm $*.o

$(LIBRARY90)(%.o): %.f90 declaration.mod 
	$(FC) $(FCFLAGS) -c $<
	ar rv $(LIBRARY90) $*.o
	ranlib $(LIBRARY90)
	-rm $*.o

declaration.mod: declaration.f90
	$(FC) $(FCFLAGS) -c $<
	/bin/rm declaration.o
