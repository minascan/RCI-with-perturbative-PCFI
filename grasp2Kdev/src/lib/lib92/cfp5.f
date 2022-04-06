************************************************************************
*                                                                      *
      SUBROUTINE CFP5 (NEL,IJD,IVD,IJP,IVP,COEFP)
*                                                                      *
*   Table look-up for fractional parentage coefficients  of equival-   *
*   ent electrons with j = 5/2. See listing of CFP for argument list.  *
*                                                                      *
*                                           Last update: 16 Oct 1994   *
*                                                                      *
*   These CFP are coordinate with theprogram Xlsj                      *
*   NIST 2005.11.15                              Gediminas Gaigalas    *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H,O-Z)
*
      DIMENSION IJ(3,3),IV(3,3),NUM3(3,3),NORM(3)
*
*   1.0 Tables of data
*
*   1.1 Floating point constants
*
      PARAMETER (D0 = 0.0D 00,
     :           D1 = 1.0D 00,
     :           D7 = 7.0D 00)
*
*   1.2 Table for CFPs.
*
      DATA IJ/5,0,5,0,4,3,0,8,9/,
     :     IV/1,0,1,8,2,3,8,2,3/,
C GG tikras     :   NUM3/-4,0,0,5,-5,3,9,2,-11/,
     :   NUM3/4,0,0,-5,-5,3,-9,2,-11/,
     :   NORM/18,7,14/
*
*   2.0 Locate entry in CFP table.
*
      IF (NEL .LE. 0) GOTO 12
      IF (NEL .GE. 4) GOTO 1
*
      N = NEL
      IJ1 = IJD
      IV1 = IVD
      IJ2 = IJP
      IV2 = IVP
      GOTO 2
    1 IF (NEL .GT. 6) GOTO 12
      N = 7-NEL
      IJ1 = IJP
      IV1 = IVP
      IJ2 = IJD
      IV2 = IVD
*
*   2.1 Find 'daughter' index.
*
    2 K = 0
    3 K = K+1
      IF (K .GT. 3) GOTO 12
      IF (IJ(N,K) .NE. IJ1) GOTO 3
      IF (IV(N,K) .NE. IV1) GOTO 3
      KD = K
*
*   2.2 Find 'parent' index.
*
      IF (N .NE. 1) GOTO 4
      IF (IV2 .NE. 0) GOTO 12
      IF (IJ2 .EQ. 0) GOTO 6
      GOTO 12
    4 K = 0
    5 K = K+1
      IF (K .GT. 3) GOTO 12
      IF (IJ(N-1,K) .NE. IJ2) GOTO 5
      IF (IV(N-1,K) .NE. IV2) GOTO 5
      KP = K
*
*   3.0 Compute coefficients.
*
*   3.1 Table look-up
*
      GOTO (6,6,7),N
    6 COEFP = D1
      GOTO 10
    7 CONTINUE
      COEFP = DBLE (NUM3(KD,KP))
      DENOM = DBLE (NORM(KD))
      IF (COEFP) 9,11,8
    8 CONTINUE
      COEFP = SQRT (COEFP/DENOM)
      GOTO 10
    9 CONTINUE
      COEFP = -SQRT (-COEFP/DENOM)
*
*   3.2 Insert additional factors for hole states
*
   10 IF (NEL .LE. 3) GOTO 11
      DNEL = NEL
      FACT = ((D7-DNEL)/DNEL)
     :     *(D1+DBLE (IJP))/(D1+DBLE (IJD))
      COEFP = COEFP* SQRT (FACT)
C GG tikras      IS = IABS ((IJD-IJP-IVD+IVP)/2)
      IS = IABS ((IJD-IJP+IVD-IVP)/2+3)
      IF (MOD (IS,2) .EQ. 0) GOTO 11
      COEFP = -COEFP
   11 RETURN
*
*   4.0 Fault mode section.
*
   12 COEFP = D0
      WRITE (*,300) NEL,IJD,IVD,IJP,IVP
      STOP
*
  300 FORMAT ('CFP5: Error in trying to compute a CFP',
     :        ' for a state with ',1I2,' electrons with j = 5/2;'
     :       /' IJD = ',1I2,', IVD = ',1I2,
     :       ', IJP = ',1I2,', IVP = ',1I2,'.')
*
      END
