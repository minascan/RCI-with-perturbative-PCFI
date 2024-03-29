************************************************************************
*                                                                      *
      SUBROUTINE RIS_CAL(NAME)
*                                                                      *
*   This routine controls the main sequence of routine calls for the   *
*   calculation  of the MS parameters, the electron density at the     *
*   origin and radial expectation values.                              *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC, GETYN                          *
*               [SMS92]: RINTDENS, VINTI                               *
*                                                                      *
*   Written by Per Jonsson                                             *
*                                                                      *
*                                         Last revision: 10 Nov 1995   *
*                                                                      *
*   Modified by C. Naz\'e  Feb. 2011                                   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = 600)
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ,
Cww     :        PINDTE,PVALTE
      POINTER (PNTRIQ,RIQDUMMY)
      POINTER (PINDTE,INDTEDUMMY)
      POINTER (PVALTE,VALTEDUMMY)
      CHARACTER*11 CNUM
      CHARACTER*4 JLBL,LABJ,LABP
      CHARACTER*2 CK,NH
      CHARACTER*24 NAME
      
      LOGICAL GETYN,FIRSTT,LDBPA,VSH,NUCDE,SMSSH,YES,AVAIL_TB,AVAIL_OB
      LOGICAL LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
*
      DIMENSION VINT(NNNW,NNNW),VINT2(NNNW,NNNW),TSHELL(NNNW)
      DIMENSION DINT1(NNNW,NNNW),DINT2(NNNW,NNNW),DINT3(NNNW,NNNW)
      DIMENSION DINT4(NNNW,NNNW),DINT5(NNNW,NNNW),DINT6(NNNW,NNNW)
      DIMENSION DINT7(NNNW,NNNW)
*
      POINTER (PNSMS1,SMSC1(1))
      POINTER (PNSMS2,SMSC2(1))
      POINTER (PNDENS1,DENS1(1))
      POINTER (PNDENS2,DENS2(1))
      POINTER (PNDENS3,DENS3(1))
      POINTER (PNDENS4,DENS4(1))
      POINTER (PNDENS5,DENS5(1))
      POINTER (PNDENS6,DENS6(1))
      POINTER (PNDENS7,DENS7(1))
      POINTER (PLABEL,LABEL(6,1))
      POINTER (PCOEFF,COEFF(1))
      POINTER (PNEVAL,EVAL(1))
      POINTER (PNEVEC,EVEC(1))
      POINTER (PNIVEC,IVEC(1))
      POINTER (PIATJP,IATJPO(1)),(PIASPA,IASPAR(1))
*
CGG      EXTERNAL COR,CORD
      COMMON/DEBUGA/LDBPA(5)
     :      /DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /DEF1/ATW,IONCTY,NELEC,Z
     :      /DEF3/EMPAM,RBCM
     :      /DEF9/CVAC,PI
     :      /DEF10/AUCM,AUEV,CCMS,FASI,FBSI
     :      /DEF11/FMTOAU,AUMAMU
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /FOPARM/ICCUT
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /JLABL/JLBL(32),LABJ(32),LABP(2)
     :      /NPAR/PARM(2),NPARM
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /SYMA/PIATJP,PIASPA
     :      /OPT6/NTC(10)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /TEILST/NDTEA,NTEI,PINDTE,PVALTE,FIRSTT
     :      /BUFFER/NBDIM,PLABEL,PCOEFF,NVCOEF
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /SMS1/PNSMS1,PNSMS2,PNDENS1,PNDENS2,PNDENS3,PNDENS4
     :          ,PNDENS5,PNDENS6,PNDENS7
     :      /TIME/TIMMAX,TIMECK,TIME1,TIME2

*
*   Matrix elements smaller than CUTOFF are not accumulated
*
      PARAMETER (CUTOFF = 1.0D-10)
*
*   Allocate storage for local arrays
*
      CALL ALLOC (PNSMS1,NVEC,8)
      CALL ALLOC (PNSMS2,NVEC,8)
      CALL ALLOC (PNDENS1,NVEC,8)
      CALL ALLOC (PNDENS2,NVEC,8)
      CALL ALLOC (PNDENS3,NVEC,8)
      CALL ALLOC (PNDENS4,NVEC,8)
      CALL ALLOC (PNDENS5,NVEC,8)
      CALL ALLOC (PNDENS6,NVEC,8)
      CALL ALLOC (PNDENS7,NVEC,8)

      CALL STARTTIME (ncount1, 'RIS')

*
*   Initialise
*
      DO 1 I = 1,NVEC
        SMSC1(I) = 0.0D 00
        SMSC2(I) = 0.0D 00
        DENS1(I) = 0.0D 00
        DENS2(I) = 0.0D 00
        DENS3(I) = 0.0D 00
        DENS4(I) = 0.0D 00
        DENS5(I) = 0.0D 00
        DENS6(I) = 0.0D 00
        DENS7(I) = 0.0D 00
    1 CONTINUE

*C*      VSH   = .TRUE.
*C*      SMSSH = .TRUE.
*
*   Calculate all integrals needed for the volume shift calc.
*
*C*      IF (VSH) THEN
      DO 5 I = 1,NW
        DO 4 J = 1,NW
          IF (NAK(I).EQ.NAK(J)) THEN
            DINT1(I,J) = RINTDENS(I,J)
CGG                  DINT2(I,J) = RINTI(I,J,1)
            CALL RINTI_NMS(I,J,DINT2(I,J),DINT7(I,J))
            DINT3(I,J) = RINT(I,J,1)
            DINT4(I,J) = RINT(I,J,2)
            DINT5(I,J) = RINT(I,J,-1)
            DINT6(I,J) = RINT(I,J,-2)
        ELSE
          DINT1(I,J) = 0.0D 00
          DINT2(I,J) = 0.0D 00
          DINT3(I,J) = 0.0D 00
          DINT4(I,J) = 0.0D 00
          DINT5(I,J) = 0.0D 00
          DINT6(I,J) = 0.0D 00
          DINT7(I,J) = 0.0D 00
        ENDIF
    4   CONTINUE
    5 CONTINUE
*C*      ENDIF
*
*   Calculate and save the Vinti integrals
*
*C*      IF (SMSSH) THEN
      DO 7 I = 1,NW
        DO 6 J = 1,NW
          IF (I.NE.J) THEN
            RCRE = CRE(NAK(I),1,NAK(J))
            IF (DABS(RCRE) .GT. CUTOFF) THEN
              VINT  (I,J) = VINTI(I,J)
              VINT2(I,J) = VINT(I,J) 
     :                   + RINT_SMS2(I,J)/RCRE
     :                   + RINT_SMS3(I,J)
            ELSE
             VINT (I,J) = 0.0D 00  
             VINT2(I,J) = 0.0D 00  
            ENDIF
          ELSE
            VINT (I,J) = 0.0D 00  
            VINT2(I,J) = 0.0D 00  
          ENDIF
    6   CONTINUE
    7 CONTINUE
*C*      ENDIF

*
*   See if the appropriate angular data is available. If so,
*   then read the angular files and perform the calculation.
*   If not, the user can choose to save them or not
*
      DOIT_OB = 0
      DOIT_TB = 0
      CALL ANGDATA(NAME,AVAIL_OB,1)
      CALL ANGDATA(NAME,AVAIL_TB,2) 
      IF ((.NOT. AVAIL_OB) .AND. (.NOT. AVAIL_TB)) THEN
        PRINT *,' Save ang. coefficients of one- and two-body op.?'
        YES = GETYN ()
        PRINT *
!C        IF(.NOT. YES) THEN
!C          PRINT *,' Save ang. coefficients of one-body op. ONLY?'
!C          YES = GETYN ()
!C          PRINT *
!C          IF(.NOT. YES) THEN
!C            PRINT *,' Save ang. coefficients of two-body op. ONLY?'
!C            YES = GETYN ()
!C            PRINT *
!C            IF(YES) DOIT_TB = 1
!C          ELSE
!C            DOIT_OB = 1
!C          ENDIF
!C        ELSE
       IF(YES) THEN
          DOIT_OB = 1
          DOIT_TB = 1
        ENDIF
      ELSEIF (.NOT. AVAIL_OB) THEN
        PRINT *,' Save ang. coefficients of one-body op. ?'
        YES = GETYN ()
        PRINT *
        IF(YES) DOIT_OB = 1
      ELSEIF (.NOT. AVAIL_TB) THEN
        PRINT *,' Save ang. coefficients of two-body op. ?'
        YES = GETYN ()
        PRINT *
        IF(YES) DOIT_TB = 1
      ENDIF
      IF (AVAIL_OB) THEN
        J = INDEX(NAME,' ')
        OPEN(UNIT=50,FILE = NAME(1:J-1)//'.IOB',STATUS='UNKNOWN'
     :      ,FORM='UNFORMATTED')
*    Read angular data from file and compute matrix elements
        CALL DENSREAD(DINT1,DINT2,DINT3,DINT4,DINT5,DINT6,DINT7)
      ELSE
*Check if the user wants to save one-body ang. coeff.
        IF(DOIT_OB .EQ.1) THEN
          J = INDEX(NAME,' ')
          OPEN(UNIT=50,FILE = NAME(1:J-1)//'.IOB',STATUS='UNKNOWN'
     :      ,FORM='UNFORMATTED')
        ENDIF
        CALL DENSNEW(DOIT_OB,DINT1,DINT2,DINT3,DINT4,DINT5,DINT6,DINT7)
      ENDIF
*
*
*
      IF (AVAIL_TB) THEN
        J = INDEX(NAME,' ')
        OPEN(UNIT=51,FILE = NAME(1:J-1)//'.ITB',STATUS='UNKNOWN'
     :    ,FORM='UNFORMATTED')
        CALL SMSREAD(VINT,VINT2)
      ELSE
* Check if the user wants to save two-body ang. coeff.
        IF(DOIT_TB.EQ.1) THEN
          J = INDEX(NAME,' ')
          OPEN(UNIT=51,FILE = NAME(1:J-1)//'.ITB',STATUS='UNKNOWN'
     :       ,FORM='UNFORMATTED')
        ENDIF
         CALL SMSNEW(DOIT_TB,VINT,VINT2) 
      ENDIF
*
*   Printouts
* 
C      WRITE (24,301) CUTOFF
      WRITE (24,308)
      DO 18 I = 1,NVEC
         WRITE(24,314)
         WRITE (24,313) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
     :                  DENS7(I),(DENS2(I)-DENS7(I)),DENS2(I)
         HzNMSu=DENS7(I)*AUMAMU*AUCM*CCMS
         if (dabs(HzNMSu) .LE. (10**8)) THEN
             WRITE (24,315) HzNMSu*(1.0D-06),
     :      (DENS2(I)-DENS7(I))*AUMAMU*AUCM*CCMS*(1.0D-06),
     :       DENS2(I)*AUMAMU*AUCM*CCMS*(1.0D-06)
         else
             WRITE (24,316) HzNMSu*(1.0D-09),
     :      (DENS2(I)-DENS7(I))*AUMAMU*AUCM*CCMS*(1.0D-09),
     :       DENS2(I)*AUMAMU*AUCM*CCMS*(1.0D-09)
         endif
   18 CONTINUE

      WRITE (24,302)
      DO 14 I = 1,NVEC
        WRITE(24,314)
        WRITE (24,313) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
     :                 SMSC1(I),SMSC2(I)-SMSC1(I),SMSC2(I)
         HzSMSu=SMSC1(I)*AUMAMU*AUCM*CCMS
         if (dabs(HzSMSu) .LE. (10**8)) THEN
           WRITE (24,315) HzSMSu*(1.0D-06),
     :      (SMSC2(I)-SMSC1(I))*AUMAMU*AUCM*CCMS*(1.0D-06),
     :      SMSC2(I)*AUMAMU*AUCM*CCMS*(1.0D-06)
         else
           WRITE (24,316) HzSMSu*(1.0D-09),
     :      (SMSC2(I)-SMSC1(I))*AUMAMU*AUCM*CCMS*(1.0D-09),
     :      SMSC2(I)*AUMAMU*AUCM*CCMS*(1.0D-09)
         endif
   14 CONTINUE
      WRITE (24,307)
      DO 16 I = 1,NVEC
        WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
     :                 DENS1(I)
   16 CONTINUE
!      WRITE (24,308)
!      DO 18 I = 1,NVEC
!         WRITE(24,314)
!         WRITE (24,313) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
!     :                  DENS7(I),(DENS2(I)-DENS7(I)),DENS2(I)
!         HzNMSu=DENS7(I)*AUMAMU*AUCM*CCMS
!         if (dabs(HzNMSu) .LE. (10**8)) THEN
!             WRITE (24,315) HzNMSu*(1.0D-06),
!     :      (DENS2(I)-DENS7(I))*AUMAMU*AUCM*CCMS*(1.0D-06),
!     :       DENS2(I)*AUMAMU*AUCM*CCMS*(1.0D-06)
!         else
!             WRITE (24,316) HzNMSu*(1.0D-09),
!     :      (DENS2(I)-DENS7(I))*AUMAMU*AUCM*CCMS*(1.0D-09),
!     :       DENS2(I)*AUMAMU*AUCM*CCMS*(1.0D-09)
!         endif
!   18 CONTINUE
      WRITE (24,309)
      DO 20 I = 1,NVEC
         WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
     :                  DENS3(I)
   20 CONTINUE
      WRITE (24,310)
      DO 22 I = 1,NVEC
         WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
     :                  DENS4(I)
   22 CONTINUE
      WRITE (24,311)
      DO 24 I = 1,NVEC
         WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
     :                  DENS5(I)
   24 CONTINUE
      WRITE (24,312)
      DO 26 I = 1,NVEC
         WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
     :                  DENS6(I)
   26 CONTINUE
*
*   Dealloc
*
      CALL DALLOC (PNSMS1)
      CALL DALLOC (PNSMS2)
      CALL DALLOC (PNDENS1)
      CALL DALLOC (PNDENS2)
      CALL DALLOC (PNDENS3)
      CALL DALLOC (PNDENS4)
      CALL DALLOC (PNDENS5)
      CALL DALLOC (PNDENS6)
      CALL DALLOC (PNDENS7)

      CALL STOPTIME (ncount1, 'RIS')
      RETURN
*
  301 FORMAT (//' CUTOFF set to ',1PD22.15)
  302 FORMAT (//' Level  J Parity  Specific mass shift parameter')
  303 FORMAT (1X,I3,5X,2A4,3X,D20.10)
  307 FORMAT (//' Electron density in atomic units'
     :        //' Level  J Parity',8X,'DENS (a.u.)'/)
  308 FORMAT (//' Level  J Parity  Normal mass shift parameter')
  309 FORMAT (//' Radial expectationvalue'
     :        //' Level  J Parity',8X,'<r> (a.u.)'/)
  310 FORMAT (//' Radial expectationvalue'
     :        //' Level  J Parity',8X,'<r2> (a.u.)'/)
  311 FORMAT (//' Radial expectationvalue'
     :        //' Level  J Parity',8X,'<r-1> (a.u.)'/)
  312 FORMAT (//' Radial expectationvalue'
     :        //' Level  J Parity',8X,'<r-2> (a.u.)'/)
  313 FORMAT (1X,I3,5X,2A4,3x,3D20.10,'  (a.u.)')
  314 FORMAT (/,29X,'<K^1>',13X,'<K^2+K^3>',9X,'<K^1+K^2+K^3>')
  315 FORMAT (20X,3D20.10 ,2X,'(MHz u)')
  316 FORMAT (20X,3D20.10,2X,'(GHz u)')
      END
