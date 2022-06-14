************************************************************************
*                                                                      *
      SUBROUTINE IABINTC (IA,IB, TEGRAL)
*                                                                      *
*                                                                      *
*     This routine returns I(ab) integrals.                            *
*                                                                      *
*     Written by Asimina Papoulia                 February 2017        *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)

      include 'parameters.def'
C     PARAMETER (KEY or KEYORB ? = 121)
      
Cww        INTEGER PNTRIQ
      POINTER(PNTRIQ,RIQDUMMY(1))                 
*
      POINTER (PCOEIL,INDOEI(1))
      POINTER (PCOEVL,VALOEI(1))
*
      COMMON/COEILS/NDCOEA,NCOEI,PCOEIL,PCOEVL
     :     /ORB2/NCF,NW,PNTRIQ
      
!----------------------------------------------------------------------

!      WRITE(*,*) 'NCOEI ', NCOEI
      
      KEY = NW + 1
      
! Ensure that the indices are in `canonical' order

      IF (IA.GT.IB) THEN
        ISWAP = IB
        IB = IA
        IA = ISWAP
      ENDIF

! Compute the composite (packed) index      
      INDEX = IA*KEY+IB
      
      JL = 1
      JU = NCOEI

!      WRITE(*,*) 'INDOEI(JL) ', INDOEI(JL), 'INDOEI(JU) ', INDOEI(JU)
!--------------------------------------------------------------------
      IF (INDEX.LT.INDOEI(JL).OR.INDEX.GT.INDOEI(JU)) THEN
        WRITE(*,*) 'Something wrong in iabintc'
        STOP
      ENDIF
!--------------------------------------------------------------------
! The index is within the range of the indices stored; search
! for it in the list of indices

    1 IF (JU-JL .GT. 1) THEN
        JM = (JU+JL)/2
        IF (INDOEI(JM) .GT. INDEX) THEN
          JU = JM
        ELSE
          JL = JM
        ENDIF
        GOTO 1
      ENDIF

! The range is bracketed to the extent possible

      IF (INDEX .EQ. INDOEI(JU)) THEN
        LOC = JU
      ELSEIF (INDEX .EQ. INDOEI(JL)) THEN
        LOC = JL
      ELSE
        WRITE(*,*) IA,IB,INDEX
        WRITE(*,*) 'Iab Integral not found'
        STOP
      ENDIF
      
! Return the value of the integral from storage
      TEGRAL = VALOEI(LOC)
      
      RETURN
      END
