************************************************************************
*                                                                      *
      SUBROUTINE SETHAM (myid, nprocs, jblock, ELSTO,ICSTRT, nelmntt
     &                   , atwinv,slf_en)
*                                                                      *
*   Sets up the Hamiltonian matrix and determines the average energy.  *
*
*   Serial I/O moved out; able to run on single/multiple processors
*   For this purpose a new common /setham_to_genmat2/ is created
*                                                                      *
*   Call(s) to: [LIB92]: ALCBUF, CONVRT, DALLOC, ICHOP, RKCO, TNSRJJ.  *
*               [RCI92]: BRINT, IABINTC, KEINT, RKINTC, VINT, VPINT.   *
*                                                                      *
*   Written by Farid A Parpia             Last revision: 30 Oct 1992   *
*   Block version by Xinghong He          Last revision: 15 Jun 1998   *
*   Shift diagonal elements by Per J                      March 2007   *
*                                                                      *
************************************************************************
* 
      USE PCFI_PT_MOD          ! ASIMINA
      
      IMPLICIT REAL*8          (A-H, O-Z)

      integer*8 nelmnt,nelmntt,nelmnttmp

      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      POINTER (PINDT1,INDT1DUMMY)
      POINTER (PINDT2,INDT2DUMMY)
      POINTER (PINDT3,INDT3DUMMY)
      POINTER (PINDT4,INDT4DUMMY)
      POINTER (PINDT5,INDT5DUMMY)
      POINTER (PINDT6,INDT6DUMMY)
      POINTER (PVALT1,VALT1DUMMY)
      POINTER (PVALT2,VALT2DUMMY)
      POINTER (PVALT3,VALT3DUMMY)
      POINTER (PVALT4,VALT4DUMMY)
      POINTER (PVALT5,VALT5DUMMY)
      POINTER (PVALT6,VALT6DUMMY)
      POINTER (PCOEIL,COEILDUMMY)
      POINTER (PCOEVL,COEVLDUMMY)
      POINTER (PCTEIL,CTEILDUMMY)
      POINTER (PCTEVL,CTEVLDUMMY)
      POINTER (PNEVAL,EVALDUMMY)
      POINTER (PINDKE,INDKEDUMMY)
      POINTER (PVALKE,VALKEDUMMY)
      POINTER (PIENDC,ENDCDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      POINTER (PNTJQS,JQSDUMMY)
      POINTER (PNJCUP,JCUPDUMMY)
      POINTER (PVINIL,VINILDUMMY)
      POINTER (PVINVL,VINVLDUMMY)
      POINTER (PINDVP,INDVPDUMMY)
      POINTER (PVALVP,VALVPDUMMY)
      POINTER (PNTRPF,RPFDUMMY)
      POINTER (PNTRQF,RQFDUMMY)
      POINTER (PNIVEC,NVECMXDUMMY)

      POINTER (PCTEVLRK,VALTEIRK(1))                                  
      POINTER (PCTEILRK, INDTEIRK(1))

      LOGICAL FIRST,FRSTCO,FRSTCT,FRSTKI,FRSTVI,FRSTVP,
     :        LDBPA,lshift,
     :        LFORDR,LTRANS,LVP,LSE,LNMS,LSMS,YES,GETYN
      CHARACTER*2 NH
*
      EXTERNAL BREID,CORD
*
      DIMENSION TSHELL(NNNW),SLF_EN(1)

      POINTER (PNTEMT,EMT(1))
      POINTER (PNIROW,IROW(1))
*
!----------------------------------------------------------------------------
!ASIMINA see subroutines talk --> alcbuf      
      POINTER (PLABEL,LABEL(6,1))
      POINTER (PCOEFF,COEFF(1))
!----------------------------------------------------------------------------
*
      COMMON/BCORE/ICORE(NNNW)
     :      /BILST/PINDT1,PINDT2,PINDT3,PINDT4,PINDT5,PINDT6,
     :             PVALT1,PVALT2,PVALT3,PVALT4,PVALT5,PVALT6,
     :             NDTPA(6),NTPI(6),FIRST(6)
     :      /BUFFER/NBDIM,PLABEL,PCOEFF,NVCOEF
     :      /COEILS/NDCOEA,NCOEI,PCOEIL,PCOEVL,FRSTCO
     :      /CTEILS/NDCTEA,NCTEI,PCTEIL,PCTEVL,FRSTCT
     :      /DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /DEBUGA/LDBPA(5)
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /EIGVAL/EAV,PNEVAL
     :      /FOPARM/ICCUT
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /HMAT/PNTEMT,PIENDC,PNIROW,NELMNT
     :      /KEILST/NDKEA,NKEI,PINDKE,PVALKE,FRSTKI
     :      /NCDIST/ZDIST(NNNP)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /STAT/PNTJQS,PNJCUP
     :      /STOR/KEEP(2,2)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /VINLST/NDVIN,NVINTI,PVINIL,PVINVL,FRSTVI
     :      /VPILST/NDVPA,NVPI,PINDVP,PVALVP,FRSTVP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /CTEILSRK/PCTEILRK,PCTEVLRK
     :      /BLIM/IPRERUN,NCSFPRE,COEFFCUT1,COEFFCUT2
     :      /WHERE/IMCDF
     :      /DEF10/AUCM,AUEV,CCMS,FASI,FBSI

* Per common for shift

      COMMON/hamshiftj/nshiftj(100),nasfshift(100,100),
     :       asfenergy(100,100),lshift

* Per end

*     ...For pre-run
      POINTER (PNEVEC1,EVEC1(1))
      COMMON/EIGVEC1/PNEVEC1

*
CGG      PARAMETER (KEYORB = 121)
      PARAMETER (KEY = KEYORB)
*
*   Matrix elements smaller than CUTOFF are not accumulated
*
      PARAMETER (CUTOFF = 1.0D-20)

      COMMON/setham_to_genmat2/CUTOFFtmp,
     &  NCOEItmp, NCOECtmp, NCTEItmp, NCTECtmp, NTPItmp(6), NMCBPtmp, 
     &  NCOREtmp, NVPItmp, NKEItmp, NVINTItmp, NELMNTtmp, ncftmp
*
!-----------------------------------------------------------------------
      PRINT *, 'Calling SETHAM ...'

      nelmnt = nelmntt
      
      !Per Bug fix 30/4-2013
      ATWINV = 1.D0/EMN
      !Per Bug fix 30/4-2013

      IF (IPRERUN .EQ. 2) THEN
         DO IPI = 1,NVEC
            DO IPJ = 1,NCF
               WRITE (*,*) IPI,IPJ,EVEC1(IPJ+(IPI-1)*NCF)
            ENDDO
         ENDDO
      ENDIF

*
*   Allocate storage to arrays in COMMON/BUFFER/; these are
*   used for the Coulomb and transverse two-electron integrals
*
      CALL ALCBUF (1)

*     ...Locals
      CALL alloc (pntemt, ncf, 8)
      CALL alloc (pnirow, ncf, 4)
*
      INC1 = 1
      INC2 = 1
*
*   Initialisations for contributions from the Dirac-Coulomb
*   operator
*
      KT  = 0
      IPT = 1
*
      INCOR = 1

      NCOEC = 0
*
      !FRSTCT = .TRUE.
      !NCTEI  = 0
      NCTEC   = 0

      IF (LTRANS) THEN

*        ...Initialisations for transverse interaction correction
         DO 2 I = 1, NW
            ICORE(I) = 0
            DO J = 1, NCF
               IF (ICHOP (I,J) .LE. 0) GOTO 2
            ENDDO
            ICORE(I) = 1
    2    CONTINUE

* ASIMINA
* Previously placed in the AUXBLK.F subroutine.            
         DO I = 1, 6
            FIRST(I) = .TRUE.
            NTPI(I)  = 0
         ENDDO

         NMCBP = 0
         NCORE = 0
      ENDIF

* ASIMINA 
* Load radial wavefunctions for the zero-order space and the first PCF
      PRINT *, 'Calling SETRWFA...'
      CALL SETRWFA (trim(PCFINAME(1))//'.w')      
*            
* Calculate all the needed one-electron and Rk  integrals
*   
      PRINT *, 'Calling GENINTIAB...'
      CALL GENINTIAB (myid, nprocs, ndum)
      
      PRINT *, 'Calling GENINTRK...'
      CALL GENINTRK (myid, nprocs, ndum, j2max)
*
*   If transverse interaction comput Breit integrals of type 1 and 2
*
      IF (LTRANS) THEN
         PRINT *, 'Calling GENINTBREIT1...'
         CALL GENINTBREIT1 (myid, nprocs, ndum, j2max)
         PRINT *, 'Calling GENINTBREIT2...'
         CALL GENINTBREIT2 (myid, nprocs, ndum, j2max)
      END IF
      
! ASIMINA counter initialization
      ip = 2
      
!ASIMINA  Constructing the upper triangle of the Hamiltonian matrix - 
!     due to symmetry the lower triangle will be the same
!--------------------------------------------------------------------------
! Loop over rows of the Hamiltonian matrix - distributed
!--------------------------------------------------------------------------
      icstrt = 1                       
      DO 10 ic = icstrt, ncf, nprocs
!ASIMINA-------------------------------------------------------------------
!-------------------------------------------------------------------------- 
! LOOP for reading the wavefunctions that correspond to first-order spaces
         if (ic .EQ. ICCUTBLK2(jblock,ip)) then

            CALL dalloc (PNTRPF) ! lodrwf or lodres
            CALL dalloc (PNTRQF) ! lodrwf or lodres
*     
*   Deallocate storage for the integral lists from the
*   Dirac-Coulomb operator; the storage was allocated
*   in the GENINTIAB and GENINTRK subroutines 
*            
            IF (NCOEI .GT. 0) THEN    
               CALL DALLOC (PCOEIL)
               CALL DALLOC (PCOEVL)
            ENDIF
            
            CALL DALLOC (PCTEVLRK) 
            CALL DALLOC (PCTEILRK)             
*
*   Deallocate storage for the integral lists from the
*   transverse photon interaction operator; this storage was
*   allocated in genintbreit1, genintbreit2 and brint3,...brint6.
*
            IF (LTRANS) THEN
               IF (NTPI(1) .GT. 0) THEN
                  CALL DALLOC (PINDT1)
                  CALL DALLOC (PVALT1)
               ENDIF
               IF (NTPI(2) .GT. 0) THEN
                  CALL DALLOC (PINDT2)
                  CALL DALLOC (PVALT2)
               ENDIF
               IF (NTPI(3) .GT. 0) THEN
                  CALL DALLOC (PINDT3)
                  CALL DALLOC (PVALT3)
               ENDIF               
               IF (NTPI(4) .GT. 0) THEN
                  CALL DALLOC (PINDT4)
                  CALL DALLOC (PVALT4)
               ENDIF
               IF (NTPI(5) .GT. 0) THEN
                  CALL DALLOC (PINDT5)
                  CALL DALLOC (PVALT5)
               ENDIF    
               IF (NTPI(6) .GT. 0) THEN
                  CALL DALLOC (PINDT6)
                  CALL DALLOC (PVALT6)
               ENDIF  
            ENDIF

            DO I = 3, 6
               FIRST(I) = .TRUE.
               NTPI(I)  = 0
            ENDDO            
            
* ASIMINA     
* Load the radial wavefunctions for the 'ip' PCF
*            
            PRINT *, 'Calling SETRWFA...'
            CALL SETRWFA (trim(PCFINAME(ip))//'.w')
*            
* Calculate all the needed one-electron and Rk  integrals
*            
            PRINT *, 'Calling GENINTIAB...'
            CALL GENINTIAB (myid, nprocs, ndum)
            
            PRINT *, 'Calling GENINTRK...'
            CALL GENINTRK (myid, nprocs, ndum, j2max)
*     
* If transverse interaction comput Breit integrals of type 1 and 2
*     
            IF (LTRANS) THEN
               PRINT *, 'Calling GENINTBREIT1...'
               CALL GENINTBREIT1 (myid, nprocs, ndum, j2max)
               PRINT *, 'Calling GENINTBREIT2...'
               CALL GENINTBREIT2 (myid, nprocs, ndum, j2max)
            END IF
*            
* ASIMINA
* Counter ip for indicating the first-order space index            
            ip = ip + 1
         end if         
!------------------------------------------------------------------------
         
         NELC = 0             ! counter - Number of non-zeros of this row
         
!     IF (LFORDR .AND. (IC .GT. ICCUT)) THEN
!     irstart = IC
!     ELSE
!     irstart = 1
!     ENDIF
         
!     Loop over columns of the current row
         
         irstart = 1
         DO 85 IR = irstart, IC
            
!ASIMINA lfordr is always true for perturbative approach            
            IF (LFORDR .AND. (IR .GT. ICCUTBLK2(jblock,1))) THEN 
!ASIMINA------All CSFs beyond iccut(1) are treated perturbatively---
!-------------------------------------------------------------------
               IF (IR.NE.IC) CYCLE
            END IF             
!     PER
            ELEMNT = 0.D0       ! accumulates various contributions to H 
*     
*     Generate the integral list for the matrix element of the
*     one-body operators
*     
!     ASIMINA   This is NOT USED anymore
            IF (IPRERUN .EQ. 1) THEN
               INC1 = 0
               INC2 = 0
               IF (IC.LE.NCSFPRE .OR. IC.EQ.IR) THEN 
                  INC1 = 1
               ENDIF
            ENDIF
            
            IF (IPRERUN .EQ. 2) THEN
*     
*     Diagonal elements are always included
*     Off diagonal elements are included only if the      
*     products of the weights from the prerun are larger
*     than the cutoff.
*     
               IF (IC .EQ. IR) THEN
                  INC1 = 1
                  INC2 = 1
               ELSE
                  INC1 = 0
                  INC2 = 0
               ENDIF
               DO IPI = 1,NVEC
                  PRECOEFF = 
     :                 DABS(EVEC1(IC+(IPI-1)*NCF)*EVEC1(IR+(IPI-1)*NCF))
                  IF (PRECOEFF .GT. COEFFCUT1) INC1 = 1 
                  IF (PRECOEFF .GT. COEFFCUT2) INC2 = 1 
               ENDDO
            ENDIF

! ---------------------------------------------------------------------
!     ASIMINA Begin to accumulate all contributions HERE
!     ...INC1.EQ.1 ------------>
            IF (INC1 .EQ. 1) THEN !inc1 is always 1 without PRE-RUN   
               
!     ASIMINA finds column and row and gives a and b orbitals and the
!     angular coefficient
               CALL ONESCALAR(IC,IR, IA,IB,TSHELL)
*     
*     Accumulate the contribution from the ONE-BODY operators:
*     kinetic energy, electron-nucleus interaction; update the
*     angular integral counter
*     
               IF (IA .NE. 0) THEN
                  IF (IA .EQ. IB) THEN
                     DO IA = 1,NW
!ASIMINA TCOEFF:angular part
                        TCOEFF = DBLE(TSHELL(IA))
                        IF (DABS (TCOEFF) .GT. CUTOFF) THEN
                           NCOEC = NCOEC + 1
!ASIMINA kinetic energy plus nuclear potential contribution
                           CALL IABINTC (IA, IA, TEGRAL)
                           ELEMNT = ELEMNT + TEGRAL*TCOEFF
                           IF (LNMS) THEN
                              CALL KEINT (IA,IA,TEGRAL)
                              ELEMNT = ELEMNT + TEGRAL*ATWINV*TCOEFF
                           ENDIF
                           IF (LVP) THEN
!ASIMINA vacuum polarization
                              CALL VPINT (IA, IA, TEGRAL)
                              ELEMNT = ELEMNT + TEGRAL*TCOEFF
                           ENDIF
                        ENDIF
                     ENDDO
                  ELSE
                     TCOEFF = DBLE(TSHELL(1))
                     IF (DABS (TCOEFF) .GT. CUTOFF) THEN
                        NCOEC = NCOEC + 1
                        CALL IABINTC (IA, IB, TEGRAL)
                        ELEMNT = ELEMNT + TEGRAL*TCOEFF                 
                        IF (LNMS) THEN
                           CALL KEINT (IA, IB, TEGRAL)
                           ELEMNT = ELEMNT + TEGRAL*ATWINV*TCOEFF
                        ENDIF
                        IF (LVP) THEN
                           CALL VPINT (IA, IB, TEGRAL)
                           ELEMNT = ELEMNT + TEGRAL*TCOEFF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
*     
               IBUG1 = 0
*     
*     Accumulate the contributions from the TWO-ELECTRON
*     Coulomb operator and the mass polarisation; the latter
*     is computed first because the orbital indices may be
*     permuted by RKINTC
*     
               NVCOEF = 0
*     
               CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
*     
               DO 7 I = 1, NVCOEF
                  VCOEFF = COEFF(I)
                  IF (DABS (VCOEFF) .GT. CUTOFF) THEN
                     NCTEC = NCTEC + 1
                     IF (LSMS) THEN
                        IF (LABEL(5,I) .EQ. 1) THEN
                           CALL VINT (LABEL(1,I), LABEL(3,I), TGRL1)
                           CALL VINT (LABEL(2,I), LABEL(4,I), TGRL2)
                           ELEMNT = ELEMNT - TGRL1*TGRL2*ATWINV*VCOEFF
                        ENDIF
                     ENDIF
                     CALL RKINTC (LABEL(1,I), LABEL(2,I),
     :                    LABEL(3,I), LABEL(4,I),
     :                    LABEL(5,I), TEGRAL)
                     ELEMNT = ELEMNT + TEGRAL*VCOEFF
                  ENDIF
 7             CONTINUE
               
*     
               IBUG1 = 0
               
            ENDIF               !inc1 is always 1 without PRE-RUN
!     ...INC1.EQ.1 <------------
************************************************************************
!     ...LTRANS .AND. (INC2.EQ.1) ------------>
            IF (LTRANS .AND. (INC2.EQ.1)) THEN
!     IF (INC2 .EQ. 1) THEN  !inc2 is always 1 without PRE-RUN
*     
*     Accumulate the contribution from the two-electron
*     transverse interaction operator
*     
               NVCOEF = 0
*               
               CALL RKCO_GG (IC, IR, BREID, 1, 2)               
*     
               DO 8 I = 1, NVCOEF
                  IF (DABS (COEFF(I)) .GT. CUTOFF) THEN
                     NMCBP = NMCBP + 1
                     ITYPE = ABS (LABEL(6,I))
                     IF (ITYPE .EQ. 1) THEN                 
                        CALL BRINT1 (LABEL(1,I), LABEL(2,I),
     :                       LABEL(3,I), LABEL(4,I),
     :                       LABEL(5,I), TEGRAL)
                     ELSEIF (ITYPE .EQ. 2) THEN
                        CALL BRINT2 (LABEL(1,I), LABEL(2,I), 
     :                       LABEL(3,I), LABEL(4,I),
     :                       LABEL(5,I), TEGRAL)
                     ELSEIF (ITYPE .EQ. 3) THEN
                        CALL BRINT3 (LABEL(1,I), LABEL(2,I),
     :                       LABEL(3,I), LABEL(4,I),
     :                       LABEL(5,I), TEGRAL)
                     ELSEIF (ITYPE .EQ. 4) THEN
                        CALL BRINT4 (LABEL(1,I), LABEL(2,I),
     :                       LABEL(3,I), LABEL(4,I),
     :                       LABEL(5,I), TEGRAL)
                     ELSEIF (ITYPE .EQ. 5) THEN
                        CALL BRINT5 (LABEL(1,I), LABEL(2,I),
     :                       LABEL(3,I), LABEL(4,I),
     :                       LABEL(5,I), TEGRAL)
                     ELSEIF (ITYPE .EQ. 6) THEN                       
                        CALL BRINT6 (LABEL(1,I), LABEL(2,I),
     :                       LABEL(3,I), LABEL(4,I),
     :                       LABEL(5,I), TEGRAL)
                     ENDIF 
                     CONTR = COEFF(I)*TEGRAL
                     IF (LABEL(6,I) .GT. 0) THEN
                        ELEMNT = ELEMNT + CONTR
                        
                     ELSE
!     ...It comes here only when ic=ir=1
!     clue: rkco<-breid<-talk<-label(6,i)
                        NCORE = NCORE + 1
                        ELSTO = ELSTO + CONTR
                     ENDIF
                  ENDIF
 8             CONTINUE
*     
               IBUG1 = 0
*     
!     ...ELSTO is a constant over all diagonals, thus its
!     contribution to the total energy can be added later
!     IF (IR .EQ. IC) ELEMNT = ELEMNT + ELSTO
*     
!     ENDIF   !inc2 is always 1 without PRE-RUN
            ENDIF
!     ...LTRANS .AND. (INC2.EQ.1) <------------
************************************************************************
!     
!     Store this element if it is diagonal or its value is greater than 
!     CUTOFF
!     
            IF ( (IR .EQ. IC) .OR. (DABS (ELEMNT) .GT. CUTOFF) ) THEN
               NELC       = NELC + 1
               EMT(NELC)  = ELEMNT
               IROW(NELC) = IR
            ENDIF
            
*     Per for shift
            
            IF (nshiftj(jblock).GT.0) THEN
               DO ishift = 1,nshiftj(jblock)
                  IF ((IR .EQ. IC) .AND.
     &                 (nasfshift(jblock,ishift).EQ. IR)) THEN
                     write(*,*)
                     write(*,*) 'Diagonalelement shifted for',IR
                     write(*,*) 'Energy shift in Hartree ',
     &                    asfenergy(jblock,ishift)
                     write(*,*) 'Energy shift in cm-1 ',
C     &                        2*109737.7*asfenergy(jblock,ishift)
     &                    aucm*asfenergy(jblock,ishift)
                     write(*,*)
                     EMT(NELC) = EMT(NELC) + asfenergy(jblock,ishift)
                  END IF
               END DO
            ENDIF
            
*     Per end
            
*     
 85      CONTINUE            
c     zou
c     print *, ic,SLF_EN(IC)
         IF(LSE) EMT(NELC) = EMT(NELC) + SLF_EN(IC)
c     zou
*     
*     This column is done; write it to disk
*     
         WRITE (imcdf) NELC, ELSTO, (EMT(IR), IR = 1, NELC),
     :        (IROW(IR), IR = 1, NELC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     This EAV (and the above EMT) does not have ELSTO.
         EAV = EAV + EMT(NELC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*     
!-----------------------------------------------------------------------
CFF   IF (MOD (IC, 20) .EQ. 0 .OR.
         IF (MOD (IC, 100) .EQ. 0 .OR.
     &        IC .LT. nprocs*2 .OR. IC .GT. (NCF-nprocs*2)) THEN
            PRINT *, 'Row ', IC, ': ', NELC, ' nonzero elements;'
     &           , '  block = ', jblock
         ENDIF
*     
*     Update the counter for the total number of elements
*     
         NELMNT = NELMNT + NELC         
*     
 10   CONTINUE
************************************************************************
*     
*     Deallocate storage for the arrays in /BUFFER/
*     
      CALL ALCBUF (3)
      
*     ...Locals
      CALL DALLOC (PNTEMT)
      CALL DALLOC (PNIROW)
      
!     Fill the common block /setham_to_genmat2/ for use in genmat2
      
      CUTOFFtmp = CUTOFF
      NCOEItmp = NCOEI
      NCOECtmp = NCOEC
      NCTEItmp = NCTEI
      NCTECtmp = NCTEC
      NTPItmp = NTPI
      NMCBPtmp = NMCBP
      NCOREtmp = NCORE
      NVPItmp = NVPI
      NKEItmp = NKEI
      NVINTItmp = NVINTI
      NELMNTtmp = NELMNT
      NCFtmp = NCF

      RETURN
      END
