************************************************************************
*
      SUBROUTINE setcsl (name, ncore, nblkin, idblk)
      IMPLICIT REAL*8          (A-H, O-Z)
*
*  A container which calls setcsll to open, read <name>.c file to get 
*     nblock, ncfblk(), idblk(), ncf (it is ncftot here). 
*  It then calls lib92/lodcsh to get
*     ncore, nelec, nw, np(), nak(), nkl(), nkj(), nh()
*  The file pointer points to the first CSL record after this routine.
*
*  Xinghong He 98-06-23
*
************************************************************************
*
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      CHARACTER*8 idblk(*), name*(*)

      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      INTEGER*4 IQAdum
      POINTER (PNTRIQ,IQAdum)
      CHARACTER*2 NH
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /ORB2/NCFtot,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
!-----------------------------------------------------------------------

c     CALL ALLOC   (pncfblk, nblkin+1, 4)
c     CALL SETCSLL (21, name, nblkin, nblock, ncfblk(1), ncftot, idblk)
c     CALL RALLOC  (pncfblk, nblkin+1, nblock+1, 4)

c     REWIND (21)
c     READ (21,*)

c     !..Load header of <name>.c file
c     CALL LODCSH (21, NCORE)
      open(unit=21,file=trim(name))
      call lodcslo(ncore)

      RETURN
      END
