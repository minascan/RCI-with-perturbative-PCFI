************************************************************************
*                                                                      *
      SUBROUTINE SMSNEW(DOIT,VINT,VINT2)
*                                                                      *
*   This routine controls the main sequence of routine calls for the   *
*   calculation  of the  sms parameter, the electron density at the    *
*   origin.
*                                                                      *
*   Call(s) to: [LIB92]: ALCBUF, CONVRT, GETYN                         *
*                        ITJPO,                                        *
*               [SMS92]: RINTISO, RINTDENS, VINTI                      *
*               [LIBRANG]: RKCO_GG                                     *
*                                                                      *
*   Written by Per Jonsson                                             *
*                                                                      *
*                                         Last revision: 10 Nov 1995   *
*                                                                      *
*   Modified by C. Naz\'e  Mai. 2012                                   *
*   Modified by J. Ekman   Mar. 2016                                   *      
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
*   Key is used to store the indices of a couple of orbitals
      PARAMETER (KEY = KEYORB)
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = 600)
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ,
Cww     :        PINDTE,PVALTE
      POINTER (PINDTE,INDTEDUMMY)
      POINTER (PVALTE,VALTEDUMMY)
      CHARACTER*11 CNUM
      CHARACTER*4 JLBL,LABJ,LABP
      CHARACTER*2 CK,NH
      LOGICAL GETYN,FIRSTT,LDBPA,VSH,NUCDE,SMSSH,YES
      LOGICAL LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
*
      DIMENSION VINT(NNNW,NNNW),VINT2(NNNW,NNNW)
*
      POINTER (PNSMS1,SMSC1(1))
      POINTER (PNSMS2,SMSC2(1))
**      POINTER (PNDENS1,DENS1(1))
**      POINTER (PNDENS2,DENS2(1))
**      POINTER (PNDENS3,DENS3(1))
**      POINTER (PNDENS4,DENS4(1))
**      POINTER (PNDENS5,DENS5(1))
**      POINTER (PNDENS6,DENS6(1))
      POINTER (PLABEL,LABEL(6,1))
      POINTER (PCOEFF,COEFF(1))
      POINTER (PNEVAL,EVAL(1))
      POINTER (PNEVEC,EVEC(1))
      POINTER (PNIVEC,IVEC(1))
      POINTER (PIATJP,IATJPO(1)),(PIASPA,IASPAR(1))
*
CGG      EXTERNAL COR,CORD
      EXTERNAL CORD
*
      COMMON/DEBUGA/LDBPA(5)
     :      /DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
**     :      /DEF1/ATW,IONCTY,NELEC,Z
**     :      /DEF3/EMPAM,RBCM
**     :      /DEF9/CVAC,PI
**     :      /DEF10/AUCM,AUEV,CCMS,FASI,FBSI
**     :      /DEF11/FMTOAU,B1
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /FOPARM/ICCUT
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /JLABL/JLBL(32),LABJ(32),LABP(2)
     :      /NPAR/PARM(2),NPARM
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /SYMA/PIATJP,PIASPA
     :      /ORB2/NCF,NW,PNTRIQ,NOB,NCFBL(100),NEVBL(100)            
     :      /TEILST/NDTEA,NTEI,PINDTE,PVALTE,FIRSTT
     :      /BUFFER/NBDIM,PLABEL,PCOEFF,NVCOEF
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /SMS1/PNSMS1,PNSMS2
*Cedric???,PNDENS1,PNDENS2,PNDENS3,PNDENS4,PNDENS5,PNDENS6


*
*   Matrix elements smaller than CUTOFF are not accumulated
*
      PARAMETER (CUTOFF = 1.0D-10)
*
      INCOR = 1
*
*   Allocate storage for the arrays in BUFFER
*     
      CALL ALCBUF (1)
*
*   Sweep through the Hamiltonian matrix to determine the
*   sms parameter
*
      DO II = 1,NOB
         DO 13 IC = NCFBL(II)+1,NCFBL(II+1)
*
*   Output IC on the screen to show how far the calculation has preceede
*
            CALL CONVRT (IC,CNUM,LCNUM)
            if (mod(IC,100).eq.0) then
               PRINT *, 'Column '//CNUM(1:LCNUM)//' complete;'
            end if
*
            ITJPOC = ITJPO (IC)
            DO 12 IR = IC,NCFBL(II+1)                   
*
*   Call the MCP package to generate V coefficients; ac and bd
*   are the density pairs
*
*   Initialize
*
*
*   Matrix elements are diagonal in J
*
               IF (ITJPO(IR) .EQ. ITJPOC) THEN
                  NVCOEF = 0
CGG            CALL RKCO (IC,IR,COR,CORD,INCOR)
                  CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
*
                  DO 11 I = 1,NVCOEF
                     VCOEFF = COEFF(I)
                     IF (ABS (VCOEFF) .GT. CUTOFF) THEN
                        IIA = LABEL(1,I)
                        IIB = LABEL(2,I)
                        IIC = LABEL(3,I)
                        IID = LABEL(4,I)
                        K   = LABEL(5,I)
*
*   Only K = 1   LABEL(5,I) .EQ. 1
*
                        IF (LABEL(5,I) .EQ. 1) THEN
                           IF (LDBPA(2)) THEN
                              WRITE (99,309) K,IC,IR,
     :                             NP(IIA),NH(IIA),NP(IIB),NH(IIB),
     :                             NP(IIC),NH(IIC),NP(IID),NH(IID),
     :                             VCOEFF
                           ENDIF
                           IF(DOIT.EQ.1) WRITE(51) IC,IR
*** Storage sequence
                           LAB  = ((IIA*KEY + IIC)*KEY+IIB)*KEY+IID
                           IF(DOIT.EQ.1) THEN
                              WRITE(51) VCOEFF,LAB
                           ENDIF
****
                           DO 10 J = NEVBL(II)+1,NEVBL(II+1) ! JE ADD   
                              LOC = (J-1)*NCF
                              CONTRIK1 = - EVEC(IC+LOC)*EVEC(IR+LOC)
     :                             * VCOEFF
     :                             * VINT (LABEL(1,I),LABEL(3,I))
     :                             * VINT (LABEL(2,I),LABEL(4,I))
                              CONTRI = - EVEC(IC+LOC)*EVEC(IR+LOC)
     :                             * VCOEFF
     :                             * ( VINT2(LABEL(1,I),LABEL(3,I))
     :                             * VINT (LABEL(2,I),LABEL(4,I)) 
     :                             + VINT2(LABEL(2,I),LABEL(4,I))
     :                             * VINT (LABEL(1,I),LABEL(3,I)) )/
     :                             2.0D 00
                              IF (IR.NE.IC) THEN 
                                 CONTRI = 2.0D 00 * CONTRI
                                 CONTRIK1 = 2.0D 00 * CONTRIK1
                              ENDIF
                              SMSC1(J) = SMSC1(J) + CONTRIK1
                              SMSC2(J) = SMSC2(J) + CONTRI
 10                        CONTINUE
                        ENDIF
                     ENDIF
 11               CONTINUE
               ENDIF
 12         CONTINUE
 13      CONTINUE
      END DO
      IF(DOIT.EQ.1) WRITE(51) -1
*
* Empty the buffer and close file
      IF(DOIT.EQ.1) CLOSE(51)
*
*   Deallocate storage for the arrays in BUFFER
*
      CALL ALCBUF (3)
      RETURN
 309  FORMAT (' V^[(',1I2,')]_[',1I3,',',1I3,']',
     :     ' (',1I2,1A2,',',1I2,1A2,';',
     :     1I2,1A2,',',1I2,1A2,') = ',1PD19.12)
*     
      END
