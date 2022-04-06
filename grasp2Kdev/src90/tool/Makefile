.SUFFIXES: .f90 .mod

BIN = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
SRCLIBDIR = ../lib
LAPACK90_DIR = $(GRASPLIB)
MODDIR = ${SRCLIBDIR}/libcommons
MODL9290 = ${SRCLIBDIR}/lib9290
F90 = $(FC)
F90_FLAGS = $(FC_FLAGS)
GRASPLIBS =-l9290 -lcommons90
LAPACK90_LIB = -llapack90

APP_LIBS = -L${GRASPLIB} ${GRASPLIBS}

UTIL = fical rcsfsplit rmixaccumulate rseqenergy  \
		 rseqhfs rseqtrans rtablevels rtabtransE1 uniquelbl

install: EXE 

EXE : $(BIN)/fical\
		$(BIN)/rcsfsplit\
		$(BIN)/rmixaccumulate\
		$(BIN)/rseqenergy\
		$(BIN)/rseqhfs\
		$(BIN)/rseqtrans\
		$(BIN)/rtablevels\
		$(BIN)/rtabtransE1\
		$(BIN)/runiquelbl\
   

OBJ : fical.o rcsfsplit.o rmixaccumulate.o rseqenergy.o  \
      rseqhfs.o rseqtrans.o rtablevels.o rtabtransE1.o 

$(BIN)/fical: fical.o
	$(FC) -o $(BIN)/fical fical.o 

$(BIN)/rcsfsplit: rcsfsplit.o
	$(FC) -o $(BIN)/rcsfsplit rcsfsplit.o 

$(BIN)/rmixaccumulate: rmixaccumulate.o
	$(FC) -o $(BIN)/rmixaccumulate rmixaccumulate.o 

$(BIN)/rseqenergy: rseqenergy.o
	$(FC) -o $(BIN)/rseqenergy rseqenergy.o 

$(BIN)/rseqtrans: rseqtrans.o
	$(FC) -o $(BIN)/rseqtrans rseqtrans.o 

$(BIN)/rseqhfs: rseqhfs.o
	$(FC) -o $(BIN)/rseqhfs rseqhfs.o 

$(BIN)/rtablevels: rtablevels.o
	$(FC) -o $(BIN)/rtablevels rtablevels.o 

$(BIN)/rtabtransE1: rtabtransE1.o
	$(FC) -o $(BIN)/rtabtransE1 rtabtransE1.o 

$(BIN)/runiquelbl : uniquelbl.o 
	   $(F90) -o $(BIN)/runiquelbl $(F90_FLAGS) uniquelbl.o $(APP_LIBS) \
			   -L ${LAPACK90_DIR} ${LAPACK90_LIB} $(APP_LIBS)

.f90.o:
	$(F90) -c $(F90_FLAGS) $< -I $(MODDIR) -I ${MODL9290} -I $(MODDIR) -o $@

clean:
	-rm -f *.o core 

