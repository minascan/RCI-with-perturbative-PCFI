      MODULE pack_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  10:50:08   2/14/04  
      SUBROUTINE pack (IUNPKD, ISUBSH, IPACKD) 
      USE vast_kind_param, ONLY: BYTE 
      INTEGER, INTENT(IN) :: IUNPKD 
      INTEGER, INTENT(IN) :: ISUBSH 
      INTEGER(BYTE), DIMENSION(*), INTENT(INOUT) :: IPACKD 
!VAST.../IOUNIT/ ISTDE(IN)
!...This routine performs I/O.
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
