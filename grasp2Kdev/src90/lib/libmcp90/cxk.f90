!***********************************************************************
!                                                                      *
      SUBROUTINE CXK (S,IS,KAPS,NU,K,IBR,IEX)
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      IMPLICIT NONE
!      DIMENSION IS(4),KAPS(4),S(12)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NU
      INTEGER :: K
      INTEGER, INTENT(IN) :: IBR
      INTEGER, INTENT(IN) :: IEX
      INTEGER, INTENT(IN) :: IS(4)
      INTEGER, INTENT(IN) :: KAPS(4)
      REAL(DOUBLE) , INTENT(INOUT) :: S(12)
!-----------------------------------------------
!
      STOP 'CXK: Error '
      END
