      MODULE coulom_I
      INTERFACE
!
      SUBROUTINE COULOM(J1,J2,J3,J4,L1,L2,L3,L4,K,AA)
      USE vast_kind_param, ONLY:  DOUBLE
      INTEGER, INTENT(IN) :: J1, J2, J3, J4, L1, L2, L3, L4, K
      REAL(DOUBLE),INTENT(OUT)  :: AA
      END SUBROUTINE
      END INTERFACE
      END MODULE
