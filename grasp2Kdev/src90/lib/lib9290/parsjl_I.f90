      MODULE parsjl_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  10:50:09   2/14/04  
      SUBROUTINE parsjl (MODE, NCORE, RECORD, LOC, JX, NJX, IERR) 
      INTEGER, INTENT(IN) :: MODE 
      INTEGER, INTENT(IN) :: NCORE 
      CHARACTER (LEN = 256), INTENT(IN) :: RECORD 
      INTEGER, INTENT(IN) :: LOC 
      INTEGER, DIMENSION(*), INTENT(OUT) :: JX 
      INTEGER, INTENT(OUT) :: NJX 
      INTEGER, INTENT(OUT) :: IERR 
!VAST.../ORB2/ NW(IN)
!VAST...Calls: CONVRT
!...This routine performs I/O.
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
