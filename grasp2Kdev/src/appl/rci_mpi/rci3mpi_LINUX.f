************************************************************************
************************************************************************
************************************************************************
***                                                                  ***
***                                                                  ***
***             ******    *****   ****   *****    *****              ***
***             **   **  **   **   **   **   **  **   **             ***
***             **   **  **        **   **   **       **             ***
***             ******   **        **    *****       **              ***
***             **  **   **        **      **       **               ***
***             **   **  **   **   **     **      **                 ***
***             **   **   *****   ****   **      *******             ***
***                                                                  ***
***          Relativistic Configuration-Interaction Program          ***
***                                                                  ***
***   This program is a derivative of GRASP2 (F. A. Parpia, I. P.    ***
***   Grant, and C. F. Fischer, 1990).                               ***
***                                                                  ***
***                            GRASP92                               ***
***          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
***                                                                  ***
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
      PROGRAM RCI3MPI
      IMPLICIT REAL*8          (A-H, O-Z)
*
*   Entry routine for RCI92. Controls the entire computation.          *
*                                                                      *
*   Call(s) to: [LIB92]: SETMC, SETCON.                                *
*               [RCI92]: CHKPLT, MATRIX, SETCSL, SETDBG, SETMIXmpi     *
*                        SETRES, SETSUM, STRSUM.                       *
*               [NJGRAF]: FACTT.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 15 Oct 1992   *
*   MPI version by Xinghong He            Last revision: 28 May 1998   *
*   Modified for ifort -i8 by A. Kramida 22 Mar 2016                   *
*   Modified by C. Froese  Fischer (removed imethod)      April 2016   *
*                                                                      *
************************************************************************

      CHARACTER*128 NAME,NAMESAVE, isofile

      PARAMETER (nblk0 = 20)
      CHARACTER*8 idblk(nblk0)

      LOGICAL GETYN,YES,LFORDR,LTRANS,LVP,LSE,LNMS,LSMS

      COMMON/DEFAULT/NDEF
     :      /BLIM/IPRERUN,NCSFPRE,COEFFCUT1,COEFFCUT2
     :      /WHERE/IMCDF

      EXTERNAL CONSTS
      COMMON/CONS/ZERO,HALF,TENTH,ONE,TWO,THREE,TEN

! Different options set in getcid
      COMMON/DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS


      ! ...setdbg.f
      LOGICAL LDBPA, LDBPG, LDBPR
      COMMON/DEBUGA/LDBPA(5)
     :      /DEBUGG/LDBPG(5)
     :      /DEBUGR/LDBPR(30)

      ! ...setmc
      COMMON/DEF0/TENMAX,EXPMAX,EXPMIN,PRECIS

      ! ...setcon.f
      COMMON/DEF3/EMPAM,RBCM
      COMMON/DEF9/CVAC,PI
      COMMON/DEF10/AUCM,AUEV,CCMS,FASI,FBSI
      COMMON/DEF11/FMTOAU,AUMAMU

      ! ...cslhmpi.f
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      INTEGER*4 IQAdum
      POINTER (PNTRIQ,IQAdum)
      CHARACTER*2 NH
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /ORB2/NCFtot,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)

      ! ...Even tnsrjj.f used this
       COMMON/TERMS/NROWS,ITAB(31),JTAB(32),NTAB(327)

      ! ...genintrkwrap.f
      POINTER (PCTEVLRK,VALTEIRK(1))
      POINTER (PCTEILRK, INDTEIRK(1))
      COMMON/CTEILSRK/PCTEILRK,PCTEVLRK

      ! ...genintrk.f, setham/rkintc.f
      PARAMETER (KMAX = 20)
      COMMON/KKSTART/KSTART(0:KMAX)
! Memories allocated in setmixmpi/lodmixmpi/lodstate/items tree

      ! ...lib92/items
      POINTER (PCCMIN,ICCMIN(1))
      COMMON/DEF7/PCCMIN,NCMIN,NCMAX   ! NCMAX not used throughout

      ! ...lodmixmpi.f and cslhmpi.f
      POINTER (pncfblk, ncfblkdum  )
      COMMON/hblock/nblock, pncfblk

      ! ...lodmixmpi
      POINTER (pnevblk, nevblk(1))
      POINTER (pncmaxblk, ncmaxblk(1))
      COMMON/hblock2/pnevblk, pncmaxblk

      ! ...lodmixmpi
      POINTER (pidxblk, idxblk(1))
      COMMON/blkidx/pidxblk

! Memories allocated in setres/getcid

      POINTER (piccutblk, iccutblk(1))
      COMMON/iccu/piccutblk

      COMMON/iounit/istdi,istdo,istde

! mpi common is written here and kept unchanged through out.
      INCLUDE 'mpif.h'
      INTEGER myid, nprocs, ierr, lenhost
      COMMON /mpi/ myid, nprocs, ierr
      CHARACTER host*(MPI_MAX_PROCESSOR_NAME), idstring*3

! Things for timing
      INTEGER   ncount1

! cpath uses
      INTEGER       lenperm, lentmp
      CHARACTER*128 startdir, permdir, tmpdir
!-----------------------------------------------------------------------
cAK Handling the -i8 option of ifort and -fdefault-integer-8 option of gfortran
cff   Set default integer length to 4
      ISIZE = 4
      I_MPI = MPI_INTEGER
      if (ISIZE.EQ.8) I_MPI = MPI_INTEGER8
c

*=======================================================================
*  Start mpi --- get processor info: myid, nprocs, host name; and print
*=======================================================================

      CALL startmpi2 (myid, nprocs, host, lenhost, ncount1,
     &                     startdir, permdir, tmpdir, 'RCIMPI')
      WRITE (idstring, '(I3.3)') myid

      lenperm = LEN_TRIM (permdir)
      lentmp = LEN_TRIM (tmpdir)

*=======================================================================
*  Get name of the state (used in files like <name>.c, <name>.s)
*  info used only on node-0
*=======================================================================

      lenname = 0     ! Otherwise it might cause trouble for myid!=0
      IF (myid .EQ. 0) THEN
         DO
cAK            WRITE (istde,'(A)',ADVANCE='NO') 'Name of state: '  ! AK: This does not work correctly in ifort/openmpi
            WRITE (istde,*) 'Name of state: '
            READ (*,'(A)') NAME
            K = INDEX (NAME,' ')
            IF (K .GT. 1) EXIT
            WRITE (istde,*) 'Name may not start with a blank. redo...'
         ENDDO
         NAMESAVE = NAME

C         open(unit=734,file=trim(name)//'.clog',status='unknown')
C         write(734,'(a)') 'y            ! Default settings'
C         write(734,'(a)') trim(name)

         !...Form the full name of the files used on node-0

         lenname = LEN_TRIM (NAME)
         isofile = permdir(1:lenperm) // '/isodata '
         NAME    = permdir(1:lenperm) // '/' // NAME(1:lenname) // ' '
         lenname = lenperm + lenname + 1
         open(unit=734,file=trim(name)//'.clog',status='unknown')
         write(734,'(a)') trim(namesave)
      ENDIF

*=======================================================================
   99 CONTINUE

*=======================================================================
*  Get NDEF on node-0 and then send to all nodes
*=======================================================================

      IF (myid .EQ. 0) THEN
         WRITE (istde,*)
cAK         WRITE (istde,'(A)',ADVANCE='NO') 'Default settings?  (y/n) '  i
! AK:    This does not work correctly in ifort/openmpi
         WRITE (istde,*) 'Full interaction?  (y/n) '
         YES = GETYN ()
         IF (YES) THEN
            NDEF = 0
            write(734,'(A)') 'y            ! Full interaction'
         ELSE
            NDEF = 1
            write(734,'(A)') 'n            ! Full interaction'
         ENDIF
cAK         WRITE (istde,'(A)',ADVANCE='NO') 'imethod ?  (1-4) '  
! AK:    This does not work correctly in ifort/openmpi
!        WRITE (istde,*) 'imethod ?  (1-4) '
!        write (istde,*)
!        write (istde,'(A)') '  1.   LAPACK (dense, in memory)'
!        write (istde,'(A)') '  2.   Davidson (dense, in memory)'
!        write (istde,'(A)') '  3.   Davidson (sparse, in memory)'
!        write (istde,'(A)') '  4.   Davidson (sparse, on disk)'
!        READ (*,*) imethod
      ENDIF
      CALL MPI_Bcast (NDEF,1,I_MPI,0,MPI_COMM_WORLD,ierr)
!     CALL MPI_Bcast (imethod,1,I_MPI,0,MPI_COMM_WORLD,ierr)




      imcdf = 26	! Unit for rci.res file
      IPRERUN = 0
*=======================================================================

*=======================================================================
*
*  Checks and settings... Mostly done in backyard.
*
*    CHKPLT - compatibility of plant substitutions
*    SETDBG - Debug output control parameters - all set to false
*    SETMC - machine- and installation- dependent constants
*    SETCON - physical constants
*    SETSUM - open the summary file
*    cslhmpi - load header of the CSL file
*    SETRES - setup the .res  file, one for each node
*    SETMIXmpi - set up mixing coefficients file AND process user input
*    FACTT - table of logarithms of factorials setup
*    SETCSLmpi - open, check, load data from, and close the  .csl  file
*=======================================================================

      CALL CHKPLT ('RCIMPI')
      CALL SETDBG
      CALL SETMC
      CALL SETCON
      IF (myid .EQ. 0)  CALL SETSUM (NAME)
      CALL cslhmpi (name(1:lenname) // '.c', ncore, nblk0, idblk)
      CALL SETRES (isofile, name(1:lenname) // '.w', idblk)
      CALL SETMIXmpi (NAME, idblk)
      IF (myid .EQ. 0)  CALL STRSUM
      CALL FACTT

*=======================================================================
*   Calculate all the needed Rk integrals - genintrkwrap calls genintrk
*   to compute the integrals (each processor allocates the same amount
*   of memory, but only computes part of the integrals); and then do the
*   "MPI_Allreduce" stuff.
*=======================================================================

      CALL GENINTRKwrap (myid, nprocs, j2max)

*
*   If transverse interaction comput Breit integrals of type 1 and 2
*
      IF (LTRANS) THEN
         CALL GENINTBREIT1WRAP (myid, nprocs, j2max)
         CALL GENINTBREIT2WRAP (myid, nprocs, j2max)
      END IF

*=======================================================================
*   Proceed with the CI calculation
*=======================================================================

!cff  CALL MATRIX (ncore, (j2max), imethod)
      CALL MATRIX (ncore, (j2max))

      IF (IPRERUN .EQ. 1) THEN
         IPRERUN = 2
         GOTO 99
      ENDIF

*=======================================================================
*  Execution finished; Statistics output
*=======================================================================

      CALL stopmpi2 (myid, nprocs, host, lenhost,
     &                     ncount1, 'RCI3MPI')

      STOP
      END
