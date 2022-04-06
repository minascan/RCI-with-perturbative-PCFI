*******************************************************************************
*                                                                             *
      SUBROUTINE EDENSITYFIT(XVEC,YVEC,Z,PAR,NRNUC,F,RHO,RES)
*                                                                             *
*   This routine fits polynomial b(r) = b0 + b2*r^2 + b4*r^4 + b6*r^6 to      *
*   electron density rho(r) within nuclear volume using a least squares       *
*   method and return electronic factors.                                     *
*                                                                             *
*   Written by Jorgen Ekman Mar. 2016                                         *
*                                                                             *
*******************************************************************************
*      
      IMPLICIT NONE
      include 'parameters.def'
      DOUBLE PRECISION :: M(3,3), RM(3), MI(3,3), PM(3), MDET
      DOUBLE PRECISION :: X(NNNP), Y(NNNP), C(3,3), B(3), CI(3,3), P(5)
      DOUBLE PRECISION :: XVEC(NNNP), YVEC(NNNP), RHO, Z
      DOUBLE PRECISION :: CDET, AU2FM, W(NNNP), NORM, RES, F(5)
      DOUBLE PRECISION :: PI, CONST, CONST2
      DOUBLE PRECISION :: FDSUM, PAR(2), DR2(4),PARF(2),PARF2(2),R2
      INTEGER :: I, NMIN, NMAX, NRNUC

      PARAMETER (PI = 3.14159265358979d0)
      PARAMETER (AU2FM = 52917.721067d0)         ! Bohr radius in fm (CODATA 2014)    
!      PARAMETER (CONST = 27.21138505d0*AU2FM*1000.0d0)  ! To get electronic factors in meV
!      PARAMETER (CONST2 = CONST*241.7989348)  ! To get electronic factors in GHz
      PARAMETER (CONST2 = 6579683.920711d0*AU2FM)  ! To get electronic factors in GHz (CODATA 2014)
      
!     COPY ARRAYS
      X(:) = XVEC(:)
      Y(:) = YVEC(:)
      
!     DETERMINE FIRST POINT WHERE GRID IS RELIABLE CALLED NMIN
!     SINCE FIRST FEW DATA POINTS IN DENSITY ARE NOT RELIABLE DUE TO DIVISION WITH SMALL
!     R^2 VALUES (STAGGERING IS SEEN), NMIN IS THE FIRST POINT TO BE USED IN THE SUBSEQUENT LEAST SQUARES FIT
      DO I=1,NNNP
         IF(X(I+1).GT.9.D-7.AND.X(I).LE.9.D-7) THEN
            NMIN = I
            exit
         END IF
      END DO

!     SETS NMAX IN FIT TO NR
      NMAX = NRNUC          


!     DETERMINE RHO:
!     BELOW PM(1) + PM(2)*X(I)**2 + PM(3)*X(I)**4 IS FITTED TROUGH DATA POINTS
!     (X(NMIN),Y(NMIN)), (X(NMIN+2),Y(NMIN+2)), (X(NMIN+4),Y(NMIN+4)).
!     WE HAVE: RHO(0) = PM(1)

      RM(1) = Y(NMIN)
      RM(2) = Y(NMIN+2)
      RM(3) = Y(NMIN+4)

      M(1,1) = 1.0d0
      M(1,2) = X(NMIN)**2.0d0
      M(1,3) = X(NMIN)**4.0d0
      M(2,1) = 1.0d0
      M(2,2) = X(NMIN+2)**2.0d0
      M(2,3) = X(NMIN+2)**4.0d0
      M(3,1) = 1.d0
      M(3,2) = X(NMIN+4)**2.0d0
      M(3,3) = X(NMIN+4)**4.0d0
      MDET = M(1,1)*M(2,2)*M(3,3)-M(1,1)*M(2,3)*M(3,2)-
     :     M(1,2)*M(2,1)*M(3,3)+M(1,2)*M(2,3)*M(3,1)+
     :     M(1,3)*M(2,1)*M(3,2)-M(1,3)*M(2,2)*M(3,1)

!     DETERMINE INVERSE OF C MATRIX
      MI(1,1) = (M(2,2)*M(3,3)-M(2,3)*M(3,2))/MDET
      MI(1,2) = (M(1,3)*M(3,2)-M(1,2)*M(3,3))/MDET
      MI(1,3) = (M(1,2)*M(2,3)-M(1,3)*M(2,2))/MDET

      MI(2,1) = (M(2,3)*M(3,1)-M(2,1)*M(3,3))/MDET
      MI(2,2) = (M(1,1)*M(3,3)-M(1,3)*M(3,1))/MDET
      MI(2,3) = (M(1,3)*M(2,1)-M(1,1)*M(2,3))/MDET

      MI(3,1) = (M(2,1)*M(3,2)-M(2,2)*M(3,1))/MDET
      MI(3,2) = (M(1,2)*M(3,1)-M(1,1)*M(3,2))/MDET
      MI(3,3) = (M(1,1)*M(2,2)-M(1,2)*M(2,1))/MDET

!     DETERMINE PARAMETERS FROM FIT
      PM(1) = MI(1,1)*RM(1)+MI(1,2)*RM(2)+MI(1,3)*RM(3)
      PM(2) = MI(2,1)*RM(1)+MI(2,2)*RM(2)+MI(2,3)*RM(3)
      PM(3) = MI(3,1)*RM(1)+MI(3,2)*RM(2)+MI(3,3)*RM(3)

!     Finallay RHO(0) in au^{-3} is determined
      RHO = PM(1)

!     START LEAST SQUARE FIT PROCEDURE FOR DATA POINTS 
!     (X(NMIN),Y(NMIN)), (X(NMIN+1),Y(NMIN+1)),...,(X(N),Y(N))
 
      NORM = -(Y(NMAX)-RHO)/AU2FM**3.0d0

      DO I=NMIN,NMAX
         X(I) = AU2FM*X(I)
         Y(I) = (Y(I)-RHO)/AU2FM**3.0d0
         Y(I) = Y(I)/NORM
         W(I) = X(I)    ! WEIGHTS SET TO R(I) TO COMPENSATE FOR EXPONENTIAL GRID DENSITY
      END DO

!     THE DATA POINTS, SUBTRACTED SO THAT Y(I) = Y(I)-RHO, ARE FITTED TO POLYNOMIAL:
!     P(2)*X(I)**2 + P(3)*X(I)**4 + P(4)*X(I)**6      
!     DETERMINE B_l AND C_{kl} MATRIX ELEMENTS
      B(:) = 0.0d0
      C(:,:) = 0.0d0
      DO I=NMIN,NMAX
         B(1) = B(1) + Y(I)*X(I)**2.0d0*W(I)
         B(2) = B(2) + Y(I)*X(I)**4.0d0*W(I)
         B(3) = B(3) + Y(I)*X(I)**6.0d0*W(I)
         C(1,1) = C(1,1) + X(I)**4.0d0*W(I)
         C(2,2) = C(2,2) + X(I)**8.0d0*W(I)
         C(3,3) = C(3,3) + X(I)**12.0d0*W(I)
         C(1,2) = C(1,2) + X(I)**6.0d0*W(I)
         C(1,3) = C(1,3) + X(I)**8.0d0*W(I)
         C(2,3) = C(2,3) + X(I)**10.0d0*W(I)
      END DO
      C(2,1) = C(1,2)
      C(3,1) = C(1,3)
      C(3,2) = C(2,3)
         
!     COMPUTE DETERMINANT
      CDET = C(1,1)*C(2,2)*C(3,3)-C(1,1)*C(2,3)*C(3,2)-
     :     C(1,2)*C(2,1)*C(3,3)+C(1,2)*C(2,3)*C(3,1)+
     :     C(1,3)*C(2,1)*C(3,2)-C(1,3)*C(2,2)*C(3,1)
      
!     COMPUTE INVERSE OF C MATRIX
      CI(1,1) = (C(2,2)*C(3,3)-C(2,3)*C(3,2))/CDET
      CI(1,2) = (C(1,3)*C(3,2)-C(1,2)*C(3,3))/CDET
      CI(1,3) = (C(1,2)*C(2,3)-C(1,3)*C(2,2))/CDET
      
      CI(2,1) = (C(2,3)*C(3,1)-C(2,1)*C(3,3))/CDET
      CI(2,2) = (C(1,1)*C(3,3)-C(1,3)*C(3,1))/CDET
      CI(2,3) = (C(1,3)*C(2,1)-C(1,1)*C(2,3))/CDET
      
      CI(3,1) = (C(2,1)*C(3,2)-C(2,2)*C(3,1))/CDET
      CI(3,2) = (C(1,2)*C(3,1)-C(1,1)*C(3,2))/CDET
      CI(3,3) = (C(1,1)*C(2,2)-C(1,2)*C(2,1))/CDET
      
!     DETERMINE FITTING PARAMETERS
      P(2) = CI(1,1)*B(1)+CI(1,2)*B(2)+CI(1,3)*B(3)
      P(3) = CI(2,1)*B(1)+CI(2,2)*B(2)+CI(2,3)*B(3)
      P(4) = CI(3,1)*B(1)+CI(3,2)*B(2)+CI(3,3)*B(3)
      
      P(1) = RHO/AU2FM**3.0d0
      P(2) = P(2)*NORM
      P(3) = P(3)*NORM
      P(4) = P(4)*NORM

!     RESULTING ELECTRONIC FACTORS F(N) IN GHZ fm^{-2N}       
      F(1) = 2.0d0*PI/3.0d0*Z*P(1)*CONST2
      F(2) = 2.0d0*PI/10.0d0*Z*P(2)*CONST2
      F(3) = 2.0d0*PI/21.0d0*Z*P(3)*CONST2
      F(4) = 2.0d0*PI/36.0d0*Z*P(4)*CONST2

!    DETERMINE AVERAGE POINT DISCREPANCY PARAMETER      
      RES = 0.0d0
      DO I=NMIN,NMAX
         Y(I) = Y(I)*AU2FM**3.0d0*NORM+RHO
         RES = RES + (Y(I)-RHO
     :        -AU2FM**3.0d0*(P(2)*X(I)**2.0d0+P(3)*X(I)**4.0d0
     :        +P(4)*X(I)**6.0d0))**2.0d0
      END DO
      
      RES = sqrt(RES/(NMAX-NMIN+1))
      RES = RES/RHO*1000.0d0    ! In per mille of RHO
         
!     BELOW WE APPROXIMATE THE IS ENERGY. DIFFERENCE IN RADIAL MOMENTS ARE CALUCLATED BETWEEN THE 
!     CURRENT DISTRIBTUION (PARF) AND WITH A DISTRIBUTION WITH A R^2 1 FM LARGER. THE DIFFERNCE 
!     IN RADIAL MOMENTS ARE THEN MULTIPLIED WITH THE ELECTRONIC FACTORS.

      PARF(1) = PAR(1)*AU2FM    ! 50% FALL OFF RADIUS C
      PARF(2) = PAR(2)*AU2FM    ! SKIN THICKNESS a. 10% TO 90% FALL OFF DISTANCE t = 4.0*ln(3)*a

      DR2(1) = 1.0d0
      R2 = 3.d0/5.d0*PARF(1)**2.d0 
     :     + 7.d0/5.d0*PI**2.d0*PARF(2)**2.d0 + DR2(1)
      PARF2(1) =  SQRT(5.d0/3.d0*(R2 
     :     - 7.d0/5.d0*(PI*PARF(2))**2.d0))
      DR2(2) = 3.d0/7.d0*(PARF2(1)**4.d0 - PARF(1)**4.d0)
     :     + 18.d0/7.d0*(PI*PARF(2))**2.d0*(PARF2(1)**2.d0 - 
     :     PARF(1)**2.d0)
      DR2(3) = 3.d0/9.d0*(PARF2(1)**6.d0 - PARF(1)**6.d0)
     :     + 11.d0/3.d0*(PI*PARF(2))**2.d0*(PARF2(1)**4.d0 - 
     :     PARF(1)**4.d0) + 239.d0/15.d0*(PI*PARF(2))**4.d0*
     :     (PARF2(1)**2.d0 - PARF(1)**2.d0)
      DR2(4) = 3.d0/11.d0*(PARF2(1)**8.d0 - PARF(1)**8.d0)
     :     + 52.d0/11.d0*(PI*PARF(2))**2.d0*(PARF2(1)**6.d0 - 
     :     PARF(1)**6.d0) + 410.d0/11.d0*(PI*PARF(2))**4.d0*
     :     (PARF2(1)**4.d0 - PARF(1)**4.d0)
     :     + 1636.d0/11.d0*(PI*PARF(2))**6.d0*
     :     (PARF2(1)**2.d0 - PARF(1)**2.d0)
      
      FDSUM = 0.0d0
      DO I=1,4 
         FDSUM = FDSUM + F(I)*DR2(I) 
      END DO
      F(5) = FDSUM              ! IS ENERGY IN GHZ / FM^2
      
      RETURN
      
      END SUBROUTINE EDENSITYFIT
      
