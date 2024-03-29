************************************************************************
*                                                                      *
      FUNCTION ITJPOR (ICSF)
*                                                                      *
*   ITJPOR is the value of 2J+1 for CSF number ICSF.                   *
*                                                                      *
*   Call(s) to: [LIB92]: IUNPCK.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
*                                                                      *
************************************************************************
*
! Packed array declared I*4, common/iounit/ added
! XHH 1997.02.12
Cww      INTEGER PNJQSR,PNTIQR
      POINTER (PNJQSR,JQSRDUMMY)
      POINTER (PNTIQR,IQRDUMMY)
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)
*
      INTEGER*4 JCUPAR
      POINTER (PJCUPR,JCUPAR(NNNWP,1))
*
      COMMON/ORB2R/NCFR,NWR,PNTIQR
     :      /STATR/PNJQSR,PJCUPR
      COMMON/iounit/istdi,istdo,istde
*
      IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCFR)) THEN
         ITJPOR = IUNPCK (JCUPAR(1,ICSF),NNNW)
         ITJPOR = ABS (ITJPOR)
      ELSE
         WRITE(istde,*) 'ITJPOR: Argument ICSF is out of range.'
         STOP
      ENDIF
*
      RETURN
      END
