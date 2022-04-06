************************************************************************
*                                                                      *
      SUBROUTINE SETRES (isofile, rwffile, idblk)
      IMPLICIT REAL*8          (A-H, O-Z)
*                                                                      *
*   Open, check, load data from the  .res  file.                       *
*                                                                      *
*   Call(s) to: [LIB92]: GETYN, OPENFL.
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 06 Oct 1992   *
*   Modified by Xinghong                  Last revision: 23 Jun 1998   *
*   Modified for ifort -i8 by A. Kramida 22 Mar 2016            *
*                                                                      *
************************************************************************
*
      CHARACTER*(*) isofile, rwffile, idblk(*)*8
      CHARACTER(LEN=*), PARAMETER:: FORM = 'UNFORMATTED', 
     &                              STATUS = 'UNKNOWN',
     &                              RESTITLE = 'R92RES'

      CHARACTER(LEN=LEN(RESTITLE)) R92RES

      COMMON/DEFAULT/NDEF
     &      /WHERE/IMCDF

      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      POINTER (piccutblk, iccutblk(1))
      COMMON/iccu/piccutblk

      COMMON/iounit/istdi,istdo,istde

      CHARACTER DEFNAM*11, idstring*3
      LOGICAL   FOUND, GETYN, RESTRT

      INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------
cAK Handling the -i8 option of ifort and -fdefault-integer-8 option of gfortran
cff   Set default integer length to 4
      ISIZE = 4
      L_MPI = MPI_LOGICAL
      if (L_MPI.eq.8) L_MPI = MPI_LOGICAL8     !! This needs to be checked
c
!
! Compose the "rciXXX.res" file names (for each node)
!
      WRITE (idstring, '(I3.3)') myid
      DEFNAM = 'rci' // idstring // '.res'
!
! Ask if this is a restart
!
      IF (myid .EQ. 0) THEN
C         IF (NDEF .NE. 0) THEN
C            WRITE (istde,*) 'Restarting RCI92 ?'
C            RESTRT = GETYN()
C         ELSE
            RESTRT = .FALSE.
C         ENDIF
      ENDIF

      CALL MPI_Bcast (RESTRT, 1, L_MPI, 0, MPI_COMM_WORLD, ierr)
!
! Do some settings and checks
!
      IF (RESTRT) THEN
!         ...Restart, make sure file exist
         INQUIRE (FILE = DEFNAM, EXIST = FOUND)
         IF (.NOT. FOUND) THEN
            CALL stopmpi ('setres: .res does not exist', myid)
         ENDIF
      ENDIF
!
! Open the .res file
!
      CALL OPENFL (imcdf, defnam, FORM, STATUS, IERR)
      IF (IERR .NE. 0) THEN
         CALL stopmpi ('setres: Error openning .res file', myid)
      ENDIF
!
! If restart, load the contents. Otherwise generate them via getcid
!
! But first of all, iccutblk() is needed in both cases
!
      CALL alloc (piccutblk, nblock, ISIZE)

      IF (RESTRT) THEN
!        ...Check the signature of the file
         READ (imcdf, IOSTAT = IOS) R92RES
         IF ((IOS .NE. 0) .OR. (R92RES .NE. RESTITLE)) THEN
            CLOSE (imcdf)
            CALL stopmpi ('setres: Not RCI92 .res file', myid)
         ENDIF

!         ...Read and check restart information
         CALL LODRES

      ELSE

!         ...Write the file header
!         ...Generate the first part of the .res file
         WRITE (imcdf) RESTITLE
         CALL GETCID (isofile, rwffile, idblk)

      ENDIF

      RETURN
      END
