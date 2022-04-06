.SUFFIXES: .f90 .mod

LIBDIR = $(GRASP)/lib
BIN = $(GRASP)/bin
LIB = $(GRASPLIBS)
BD = $(LIBDIR)/graspBD.o
UTIL = rmixextract  rcsfblock rwfnmchfmcdf  \
       rlevels rtabtrans1 rwfnrotate \
       ris_lsj risosupport rwfnplot \
       rtabtrans2 rlevelseV rtabhfsiso \
		 wfnplot rwfnrelabel

install: EXE 
	cp rsave $(GRASP)/bin; chmod +x $(GRASP)/bin/rsave
	cp lscomp.pl $(GRASP)/bin; chmod +x $(GRASP)/bin/lscomp.pl

EXE : $(BIN)/rmixextract $(BIN)/rcsfblock  \
      $(BIN)/rwfnmchfmcdf  \
      $(BIN)/rlevels $(BIN)/rtabtrans1  \
      $(BIN)/rwfnrotate \
      $(BIN)/ris_lsj $(BIN)/risosupport \
      $(BIN)/rwfnplot \
		$(BIN)/rtabtrans2\
		$(BIN)/rlevelseV\
		$(BIN)/rtabhfsiso\
		$(BIN)/wfnplot\
		$(BIN)/rwfnrelabel\
   

OBJ : rmixextract.o rcsfblock.o  rwfnmchfmcdf.o \
      rlevels.o rtabtrans1.o format_mix.o rlevelsj.o \
      rwfnrotate.o ris_lsj.o risosupport.o \
      rwfnplot.o \
		rtabtrans2.o rlevelseV.o rtabhfsiso \\
		wfnplot rwfnrelabel
      


$(BIN)/format_mix  : format_mix.o
	$(FC) -o $(BIN)/format_mix  $(FC_LD) -L$(LIBDIR) \
                 format_mix.o $(BD) $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rmixextract: rmixextract.o
	$(FC) -o $(BIN)/rmixextract $(FC_LD) -L$(LIBDIR)  rmixextract.o \
                 $(BD) $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rcsfblock: rcsfblock.o
	$(FC) -o $(BIN)/rcsfblock $(FC_LD)  -L$(LIBDIR) rcsfblock.o \
                 $(BD) $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rwfnmchfmcdf: rwfnmchfmcdf.o
	$(FC) -o $(BIN)/rwfnmchfmcdf $(FC_LD) -L$(LIBDIR) rwfnmchfmcdf.o \
                 $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rlevels: rlevels.o
	$(FC) -o $(BIN)/rlevels $(FC_LD) -L$(LIBDIR) rlevels.o $(LIB) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rtabtrans1: rtabtrans1.o
	$(FC) -o $(BIN)/rtabtrans1 $(FC_LD) -L$(LIBDIR) rtabtrans1.o $(LIB) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rlevelsj: rlevelsj.o
	$(FC) -o $(BIN)/rlevelsj $(FC_LD) -L$(LIBDIR) rlevelsj.o $(LIB) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rwfnrotate: rwfnrotate.o
	$(FC) -o $(BIN)/rwfnrotate  $(FC_LD) -L$(LIBDIR) -I $(GRASP_INCLUDES) \
                rwfnrotate.o $(BD) $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/ris_lsj: ris_lsj.o
	$(FC) -o $(BIN)/ris_lsj  $(FC_LD) -L$(LIBDIR) -I $(GRASP_INCLUDES) \
                ris_lsj.o $(BD) $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/risosupport: risosupport.o
	$(FC) -o $(BIN)/risosupport  $(FC_LD) -L$(LIBDIR) -I $(GRASP_INCLUDES) \
                risosupport.o $(BD) $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rwfnplot: rwfnplot.o
	$(FC) -o $(BIN)/rwfnplot rwfnplot.o 

$(BIN)/rtabtrans2: rtabtrans2.o
	$(FC) -o $(BIN)/rtabtrans2 rtabtrans2.o 

$(BIN)/rlevelseV: rlevelseV.o
	$(FC) -o $(BIN)/rlevelseV $(FC_LD) -L$(LIBDIR) rlevelseV.o $(LIB) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rtabhfsiso: rtabhfsiso.o
	$(FC) -o $(BIN)/rtabhfsiso rtabhfsiso.o 

$(BIN)/wfnplot: wfnplot.o
	$(FC) -o $(BIN)/wfnplot wfnplot.o 

$(BIN)/rwfnrelabel: rwfnrelabel.o
	$(FC) -o $(BIN)/rwfnrelabel  $(FC_LD) -L$(LIBDIR) -I $(GRASP_INCLUDES) \
                rwfnrelabel.o $(BD) $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

.f90.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

