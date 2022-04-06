!*******************************************************************
!                                                                  *
      SUBROUTINE RMEAJJ(LL,IT,LQ,J,ITS,LQS,J1S,COEF)
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE CONS_C
      USE ribojj_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rmeajj11_I
      USE rmeajj9_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN)  :: LL, IT, LQ, J, ITS, LQS, J1S
      REAL(DOUBLE), INTENT(OUT) :: COEF
!      DIMENSION IDS3(1,2),IDV3(1,2),IDS5(3,3),IDV5(3,3),IDS7(6,8),
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I1, I2, IE1,  JT, JTS
      INTEGER, DIMENSION(1,2) :: IDS3, IDV3
      INTEGER, DIMENSION(3,3) :: IDS5, IDV5
      INTEGER, DIMENSION(6,8) :: IDS7, IDV7
!-----------------------------------------------
      DATA IDS3/12,20/
      DATA IDV3/1,1/
      DATA IDS5/-24,2*0,-30,120,-90,-54,-48,330/
      DATA IDV5/4*1,2*7,1,2*7/
      DATA IDS7/2*0,-40,3*0,-54,-33,-40,-65,30,0,88,-1,0,-1755,-40,  &
      0,198,-36,-72,4500,-234,-360,-78,156,0,108,-30,96,66,49,0,27,  &
      130,136,0,195,-104,-245,-624,1224,3*0,2720,-2448,-7752/
      DATA IDV7/6*1,7,2*1,7,2*1,7,2*1,77,11,1,7,11,1,77,2*11,35,5,1, &
      91,1,13,2*5,2*1,2*7,1,11,1,3*11,3*1,143,77,91/
!
      COEF=ZERO
      IF(LL.GT.37) RETURN
      IF(LL.LE.9) THEN
        IF(IT.LT.300) THEN
          IF(IMPTJJ(IT) .NE. IMPNJJ(ITS)) RETURN
          IF(IT .LT. ITS) THEN
            JT=IT
            JTS=ITS
          ELSE
            JT=ITS
            JTS=IT
          ENDIF
        ENDIF
      ENDIF
      IF(LL.EQ.1) THEN
        COEF=-2
      ELSEIF(LL.EQ.3)THEN
        I1=JT-2
        I2=JTS-3
        COEF=-DSQRT(DBLE(IDS3(I1,I2))/DBLE(IDV3(I1,I2)))
      ELSEIF(LL.EQ.5)THEN
        I1=JT-5
        I2=JTS-8
        IF(IDS5(I1,I2).GE.0) THEN
          COEF=DSQRT(DBLE(IDS5(I1,I2))/DBLE(IDV5(I1,I2)))
        ELSE
          COEF=-DSQRT(-DBLE(IDS5(I1,I2))/DBLE(IDV5(I1,I2)))
        ENDIF
      ELSEIF(LL.EQ.7)THEN
        I1=JT-11
        I2=JTS-17
        IF(IDS7(I1,I2).GE.0) THEN
          COEF=DSQRT(DBLE(IDS7(I1,I2))/DBLE(IDV7(I1,I2)))
        ELSE
          COEF=-DSQRT(-DBLE(IDS7(I1,I2))/DBLE(IDV7(I1,I2)))
        ENDIF
      ELSEIF(LL.EQ.9) THEN
        IF(IT.GT.300) THEN
          IF(ITS.LT.300) THEN
            WRITE(6,'(A)') ' ERROR IN  RMEAJJ '
            STOP
          ENDIF
          IF(LL.EQ.J1S) THEN
            CALL RMEAJJ11(IT,ITS,LL,COEF)
          ELSE
            CALL RMEAJJ11(ITS,IT,LL,COEF)
            IE1=LQ-LQS+J-J1S+LL-1
            IF((IE1/4)*4.NE.IE1)COEF=-COEF
          ENDIF
        ELSE
          CALL RMEAJJ9(IT,LQ,J,ITS,LQS,J1S,COEF)
          IF(IT .GT. ITS) THEN
            IF(MOD(LQ+J-LQS-J1S+LL-1,4) .NE. 0) COEF=-COEF
          ENDIF
          WRITE(0,'(A)') ' KLAIDA RMEAJJ SUB. '
          STOP
        END IF
      ELSE
        IF(LL.EQ.J1S) THEN
          CALL RMEAJJ11(IT,ITS,LL,COEF)
        ELSE
          CALL RMEAJJ11(ITS,IT,LL,COEF)
          IE1=LQ-LQS+J-J1S+LL-1
          IF((IE1/4)*4.NE.IE1)COEF=-COEF
        END IF
      ENDIF
      IF(LL.LT.9) THEN
        IF(IT .GT. ITS) THEN
          IF(MOD(LQ+J-LQS-J1S+LL-1,4) .NE. 0) COEF=-COEF
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE RMEAJJ
