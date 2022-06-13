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
      PROGRAM RCI92
*                                                                      *
*   Entry routine for RCI92. Controls the entire computation.          *
*                                                                      *
*   Call(s) to: [LIB92]: SETMC, SETCON.                                *
*               [RCI92]: CHKPLT, MATRIX, SETCSL, SETDBG, SETMIX,       *
*                        SETRES, SETSUM, STRSUM.                       *
*               [NJGRAF]: FACTT.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 15 Oct 1992   *
*   Updated by Xinghong He                Last revision: 23 Jun 1998   *
*   Modified by Gediminas Gaigalas for new spin-angular integration.   *
*                                                                      *
************************************************************************
* 
      USE PCFI_PT_MOD          ! ASIMINA
      
      IMPLICIT REAL*8          (A-H, O-Z)

! cpath uses

!      CHARACTER*128 NAME, tmpdir, permdir, isofile
      CHARACTER*128 tmpdir, permdir, isofile
      
CGG      PARAMETER (nblk0 = 20)
      PARAMETER (nblk0 = 50)
      CHARACTER*8 idblk(nblk0)

      LOGICAL GETYN,YES,LFORDR,LTRANS,LVP,LSE,LNMS,LSMS

      COMMON/DEFAULT/NDEF
     :      /BLIM/IPRERUN,NCSFPRE,COEFFCUT1,COEFFCUT2
     :      /WHERE/IMCDF

      EXTERNAL CONSTS
      COMMON/CONS/ZERO,HALF,TENTH,ONE,TWO,THREE,TEN

! Different options set in getcid
      COMMON/DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS

! Memories allocated in setmix/lodmix/lodstate/items tree

      ! ...lib92/items
      POINTER (PCCMIN,ICCMIN(1))
      COMMON/DEF7/PCCMIN,NCMIN,NCMAX   ! NCMAX not used throughout

      ! ...lodmix
      POINTER (pncfblk, ncfblkdum  )
      COMMON/hblock/nblock, pncfblk

      ! ...lodmix
      POINTER (pnevblk, nevblk(1))
      POINTER (pncmaxblk, ncmaxblk(1))
      COMMON/hblock2/pnevblk, pncmaxblk

      ! ...lodmix
      POINTER (pidxblk, idxblk(1))
      COMMON/blkidx/pidxblk

! Memories allocated in setres/getcid

      POINTER (piccutblk, iccutblk(1))
      COMMON/iccu/piccutblk

      COMMON/iounit/istdi,istdo,istde

! Things for timing

      INTEGER ncount1, ncount2, ncount_rate, ncount_max

      CHARACTER chdate*8, chtime*10, chzone*5
               !ccyymmdd  hhmmss.sss  Shhmm
      INTEGER  nYMDUHMSM(8)
               !Year Month Day Universal Hour Minute Sesond Millisecond

      CHARACTER str*8, msg*128
!-----------------------------------------------------------------------
                
      imcdf = 26	! Unit for rci.res file
      IPRERUN = 0
      myid = 0
      nprocs = 1
      open(UNIT=31,STATUS="SCRATCH",FORM="FORMATTED")

      write(*,*)
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*) '    RCI combined with perturbative PCFI method    '
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*) 'This is the NEW configuration interaction program '
      write(*,*) '__________________________________________________'
      write(*,*) 'Input files:  isodata, name.c, and name.pt (the   '
      write(*,*) '  (latter .pt file defines the different .w files)'
      write(*,*) 'Output files: name.cm, name.csum, name.clog, and  '
      write(*,*) '                 rci.res (can be used for restart)'
      write(*,*) '__________________________________________________'
      write(*,*)

      CALL starttime (ncount1, 'RSCF2')
!
! Start timing
!
      CALL STARTTIME (ncount1, 'RCI')
CGG      CALL SYSTEM_CLOCK (ncount1, ncount_rate, ncount_max)

CGG      CALL DATE_AND_TIME (chdate, chtime, chzone, nYMDUHMSM)
CGG      msg = ' Date: ' // chdate //
CGG     &      ' Time: ' // chtime //
CGG     &      ' Zone: ' // chzone
CGG      PRINT *, msg
!
! Get NDEF 
!
CGG         WRITE (istde,*) 'RCI2: Execution begins ...'
CGG         WRITE (istde,*)
!
! Get name of the state (used in files like <name>.c, <name>.s)
!
         DO
            WRITE (istde,'(A)',ADVANCE='NO') ' File that contains'// 
     &' zero/first order CSFs: '
            READ (*,'(A)') NAME
            K = INDEX (NAME,' ')
            IF (K .GT. 1) EXIT
            WRITE (istde,*) 'Name may not start with a blank. redo...'
         ENDDO

! Now the name of the state is known, open the log file

         if (ndef.eq.0) then
            open(unit=734, file=trim(name)//'.clog',status='unknown')
            write(734,'(a)') trim(name)
         end if


!         ...Form the full name of the files used on node-0

         lenname = LEN_TRIM (NAME)
         
         isofile = 'isodata'
C         print *, 'isofile = ', isofile(1:LEN_TRIM (isofile))
C         print *, 'name = ', name(1:LEN_TRIM (name))


! ASIMINA ---------------------------------------------------------------
! For implementing the pertubative PCFI: always NOT FULL interaction
         
!         WRITE (istde,'(A)',ADVANCE='NO') ' Full interaction? '
!         YES = GETYN ()
!         IF (YES) THEN
!            NDEF = 0
!            write(734,'(A)') 'y            ! Full interaction'
!         ELSE
            !ASIMINA 
            NDEF = 1
!            write(734,'(A)') 'n            ! Full interaction'
!         ENDIF


! ASIMINA ---------------------------------------------------------------
! Subroutine that uses the output (.pt) file from rpcfcollect program
            WRITE(*,*)
            WRITE(*,*) 'Iccut input follows '
            CALL PCFIPTINP
!------------------------------------------------------------------------
 99   CONTINUE
!     
!     Check compatibility of plant substitutions. 
!     
C     PRINT *, 'Calling CHKPLT...'
      CALL CHKPLT ('RCI92')
      
!     
! In SETDBG of this version all control logicals are set to 
! false thus no debug output will be made
!
C      PRINT *, 'Calling SETDBG...'
      CALL SETDBG
!
! Perform machine- and installation-dependent setup
!
C      PRINT *, 'Calling SETMC...'
      CALL SETMC
!
! Set up the physical constants 
!
C      PRINT *, 'Calling SETCON...'
      CALL SETCON
!
! Open summary file
!
C      PRINT *, 'Calling SETSUM...'
      CALL SETSUM (NAME)

C      PRINT *, 'Calling setcsl...'
      CALL SETCSL (name(1:lenname) // '.c', ncore, nblk0, idblk)
!
!         Set up the  .res  file; determine if this is a restart.
!ASIMINA  .res file contains the computed matrix elements   - binary file  
!
C      PRINT *, 'Calling SETRES...'
      CALL SETRES (isofile)
      !CALL SETRES (isofile, name(1:lenname) // '.w', idblk)      !ASIMINA
*
*         Open the  .mix  file; determine the eigenpairs required
*
C      PRINT *, 'Calling SETMIX...'
      CALL SETMIX (NAME, idblk)
*
*        Append a summary of the inputs to the  .sum  file     
*
      PRINT *, 'Calling STRSUM...'
      CALL STRSUM
*
*   Set up the table of logarithms of factorials
*
      PRINT *, 'Calling FACTT...'
      CALL FACTT
*
*   Proceed with the CI calculation
*
      PRINT *, 'Calling MATRIX...'

      CALL MATRIX (ncore, (j2max))

      IF (IPRERUN .EQ. 1) THEN
         IPRERUN = 2
         GOTO 99
      ENDIF

      IF (myid .EQ. 0) THEN
         PRINT *
         PRINT *
         PRINT *, 'Finish time, Statistics'
         PRINT *
      ENDIF
      CALL STOPTIME (ncount1, 'RCI')
CGG      CALL SYSTEM_CLOCK (ncount2, ncount_rate, ncount_max)
CGG		ncount2 = ncount2 - ncount1
CGG		nseconds = ncount2 / ncount_rate
CGG		WRITE (str, '(I8)') nseconds
CGG		msg = str // ' seconds '
CGG		PRINT *, msg

CGG      CALL DATE_AND_TIME (chdate, chtime, chzone, nYMDUHMSM)

CGG      msg = ' Date: ' // chdate //
CGG     &      ' Time: ' // chtime //
CGG     &      ' Zone: ' // chzone
CGG      PRINT *, msg
*
*   Print completion message
*
CGG      PRINT *, 'RCI2: Execution complete.'
*
      STOP
      END
