************************************************************************
*                                                                      *
      SUBROUTINE GETCID (isofile)
      !SUBROUTINE GETCID (isofile, rwffile, idblk)      !ASIMINA
*                                                                      *
*   Interactively determines the data governing the CI problem.        *
*   iccut is replaced by an array iccutblk(1:nblock)
*                                                                      *
*   Call(s) to: [LIB92]: GETYN, NUCPOT, RADGRD, SETQIC.                *
*               [RCI92]: SETISO, SETRWF.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
*   Block version by Xinghong He          Last revision: 15 Jun 1998   *
*                                                                      *
************************************************************************
*     
      USE PCFI_PT_MOD          ! ASIMINA 
      
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*(*) isofile
      !CHARACTER*(*) isofile, rwffile, idblk(*)*8
      
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)

      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL GETYN,LFORDR,LTRANS,LVP,LSE,LNMS,LSMS,YES,lshift
*
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      ! set in setiso/lodiso
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /NPAR/PARM(2),NPARM
     :      /NSMDAT/SQN,DMOMNM,QMOMB

      COMMON/DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /DEF2/C
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /DEF9/CVAC,PI
     :      /DEFAULT/NDEF
     :      /FOPARM/ICCUT
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPOT/ZZ(NNNP),NNUC
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /WFAC/WFACT
     :      /BLIM/IPRERUN,NCSFPRE,COEFFCUT1,COEFFCUT2
     :      /WHERE/IMCDF
     :      /QEDCUT/NQEDCUT,NQEDMAX

      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      POINTER (pnevblk, nevblk(1))
      POINTER (pncmaxblk, ncmaxblk(1))
      COMMON/hblock2/pnevblk, pncmaxblk

      POINTER (piccutblk, iccutblk(1))
      COMMON/iccu/piccutblk

      COMMON/iounit/istdi,istdo,istde
      COMMON/hamshiftj/nshiftj(100),nasfshift(100,100),
     :       asfenergy(100,100),lshift


!------------------------------------------------------------------------
      
!
! Open, check, load data from, and close the  .iso  file
! Data loaded are: EMN,Z,/NPAR/,/NSMDAT/ where /.../ means whole 
! common block.
!
C      PRINT *, 'Calling SETISO ...'
      CALL SETISO (isofile)
!
! The speed of light, if non-default then spread from node-0
! Quantities to be obtained: C
!
C      IF (NDEF .NE. 0) THEN
C         WRITE (istde,*) 'Revise the physical speed of light (',CVAC,
C     &                     ' in a.u.) ?'
C         YES = GETYN ()
C         IF (YES) THEN
C            WRITE (istde,*) 'Enter the revised value:'
C            READ *,C
C         ELSE
C            C = CVAC
C         ENDIF
C      ELSE
C         C = CVAC
C      ENDIF

       C = CVAC
      
!
! Treat some CSF's as perturbation ? Broadcasting is used only in 
! non-default mode, as the above case for speed of light.
! Quantities to be obtained: LFORDR
!
!      IF (NDEF .NE. 0) THEN
C         WRITE (istde,*) 'Treat contributions of some CSFs'
C     &,                    ' as first-order perturbations?'
C         LFORDR = GETYN ()
         LFORDR = .TRUE.
!      ELSE
!         LFORDR = .FALSE.
!      ENDIF

! Get ICCUTBLK2() from the user-input

!      IF (.NOT. LFORDR) THEN
         !...Default first
!         DO i = 1, nblock
!            DO j = 1, npcfi
!               ICCUTBLK2(i,j) = ncfblk(i)
!            ENDDO
!         ENDDO
!      ENDIF
!ASIMINA -----------------------------------------------------------
! Now that the PCFIPTINP subroutine exists I comment the following lines    
!ASIMINA -----------------------------------------------------------
!For implementing the pertubative PCFI: always NOT FULL interaction
         
      ! Let master do the i/o, then broadcast
           ! WRITE (istde,*) 'There are ', nblock, 'blocks. They are:'
           ! WRITE (istde,*) '  block     J Parity     No of CSFs'
           ! DO i = 1, nblock
           !    WRITE (istde,*) i, idblk(i)(1:5), ncfblk(i)
           ! ENDDO

!ASIMINA------------------------------------------------------------
!Several iccut: (iccut(i)) with NCSFPCFI(i) values for the relative
!PCFINAME(i) files.            
        
           ! WRITE (istde,*)
           ! WRITE (istde,*) 'Enter iccut (larger than 1) for each block'
           ! DO jblock = 1, nblock              
           !    WRITE (istde,*) 'Block ', jblock, '   ncf = ',
!     &                        ncfblk(jblock)
!     &              , ' id = ', idblk(jblock)(1:5)
               
!               DO ip = 1, npcfi                            ! new loop
!                  WRITE(istde,'(A,I3)') ' Give iccut ', ip !!
! 123              READ (istdi,*) ntmp2(ip)
!                  IF (ntmp2(ip) .GE. 0 .AND. 
!     &                ntmp2(ip) .LE. ncfblk(jblock)) THEN
!                     ICCUTBLK2(jblock,ip) = ntmp2(ip)
!                  ELSE
!                     WRITE (istde,*) 'ICCUT out of range, re-enter:'
!                     GOTO 123
!                  ENDIF
!                  write(734,*) ntmp2(ip),'! ICCUT for block',jblock
!               ENDDO                                      ! end new loop
!               
!            ENDDO
!         ENDIF

******************************************************************
!
! Pre-run ?
!
!     IF (IPRERUN .EQ. 0) THEN

!        WRITE (istde,*) ' Prerun with limited interaction?'
!        YES = GETYN ()

!        IF (YES) THEN
!           IPRERUN = 1
!           LTRANS = .FALSE.
!           LVP = .FALSE.
!           LNMS = .FALSE.
!           LSMS = .FALSE.
!           LSE = .FALSE.

!           WRITE (istde,*)  ' Give CSL cut'
!           READ *, NCSFPRE
!           WRITE (istde,*)  ' Give coefficient cut for H_0'
!           READ *, COEFFCUT1
!           WRITE (istde,*) ' Give coefficient cut for the transvers'
!    &,                  ' interaction'
!           GOTO 99
!        ENDIF
!     ENDIF
******************************************************************
!
! Include transverse ?
! Quantities to be obtained: LTRANS, WFACT
!
      WRITE (istde,*) 'Include contribution of H (Transverse)?'
      LTRANS = GETYN ()
      WRITE (istde,*) 'Modify all transverse photon frequencies?'
      YES = GETYN ()
      IF (YES) THEN
         WRITE (istde,*) 'Enter the scale factor:'
         READ *, WFACT
      ELSE
         WFACT = 1.0D 00
      ENDIF
!
! Other interactions ? One logical for each case. Done altogether
!
      WRITE (istde,*) 'Include H (Vacuum Polarisation)?'
      LVP = GETYN ()

      WRITE (istde,*) 'Include H (Normal Mass Shift)?'
      LNMS = GETYN ()

      WRITE (istde,*) 'Include H (Specific Mass Shift)?'
      LSMS = GETYN ()

      WRITE (istde,*) 'Estimate self-energy?'
      LSE = GETYN ()
      IF (LSE.EQV..TRUE.) THEN
         NQEDCUT = 1
         WRITE (istde,*)
     :'Largest n quantum number for including self-energy for orbital'
         WRITE (istde,*) 'n should be less or equal 8'
         READ *, NQEDMAX
      ELSE
         NQEDCUT = 0
      END IF

      IF (NDEF .EQ. 0) THEN
         lshift = .FALSE.
      ELSE
C         WRITE (istde,*) 'Shift diagonal matrix elements?'
C         lshift = GETYN ()
         lshift = .FALSE.
      END IF

C      IF (NDEF.EQ.0) THEN
       IF (LTRANS) THEN
        WRITE(734,'(a)') 'y            ! Contribution of H Transverse?'
       ELSE
        WRITE(734,'(a)') 'n            ! Contribution of H Transverse?'
       END IF
       IF (YES) THEN
        WRITE(734,'(a)') 'y            ! Modify photon frequencies?'
        WRITE(734,*) WFACT,'! Scale factor'
       ELSE
        WRITE(734,'(a)') 'n            ! Modify photon frequencies?'
       END IF
       IF (LVP) THEN
        WRITE(734,'(a)') 'y            ! Vacuum polarization?'
       ELSE
        WRITE(734,'(a)') 'n            ! Vacuum polarization?'
       END IF
       IF (LNMS) THEN
        WRITE(734,'(a)') 'y            ! Normal mass shift?'
       ELSE
        WRITE(734,'(a)') 'n            ! Normal mass shift?'
       END IF
       IF (LSMS) THEN
        WRITE(734,'(a)') 'y            ! Specific mass shift?'
       ELSE
        WRITE(734,'(a)') 'n            ! Specific mass shift?'
       END IF
       IF (LSE) THEN
        WRITE(734,'(a)') 'y            ! Self energy?'
        WRITE(734,*) NQEDMAX, '! Max n for including self energy'
       ELSE
        WRITE(734,'(a)') 'n            ! Self energy?'
       END IF

C      END IF

!
! Parameters controlling the radial grid
!
!  99  IF (NPARM .EQ. 0) THEN
!         RNT = EXP (-65.D0/16.D0) / Z
!         H = 0.5D0**4
!         N = MIN (220,NNNP)
!      ELSE            ! Already defined by setiso
!         RNT = 2.D-6  ! Jon Grumer (Lund, 2013)     
!         H = 5.D-2    !
!         N = NNNP     !
!      ENDIF
!      HP = 0.D0       !
*
*     Print grid parameters to screen
*     
C      PRINT*
C      PRINT*, 'DEFAULT GRID PARAMETERS'
C      PRINT*, '----------------------------------------------------'
C      PRINT*, 'RNT =', RNT 
C      PRINT*, 'H   =', H 
C      PRINT*, 'HP  =', HP
C      PRINT*, 'N   =', N
C      PRINT*, '----------------------------------------------------'
C      PRINT*
*      
C      IF ( NDEF.NE.0) THEN
C         WRITE (istde,*) 'The default radial grid parameters'
C     &,                 ' for this case are:'
C         WRITE (istde,*) ' RNT = ',RNT,';'
C         WRITE (istde,*) ' H = ',H,';'
C         WRITE (istde,*) ' HP = ',HP,';'
C         WRITE (istde,*) ' N = ',N,';'
C         WRITE (istde,*) ' revise these values?'
C         YES = GETYN ()
C
C         IF (YES) THEN
C            WRITE (istde,*) 'Enter RNT:'
C            READ *, RNT
C            WRITE (istde,*) 'Enter H:'
C            READ *, H
C            WRITE (istde,*) 'Enter HP:'
C            READ *, HP
C            WRITE (istde,*) 'Enter N:'
C            READ *, N
C         ENDIF
C      ENDIF
*
*   ACCY is an estimate of the accuracy of the numerical procedures
*
      ACCY = H**6
*
*   Set up the coefficients for the numerical procedures
*
      CALL SETQIC
*
*   Generate the radial grid and all associated arrays
*
      CALL RADGRD
*
*   Generate $- r \times V_ (r)$
*
      CALL NUCPOT
*
!-----------------------------------------------------------------------
!ASIMINA This is moved to the SETHAM subroutine      
*   Load the radial wavefunctions
*
!      CALL SETRWFA (rwffile)
!-----------------------------------------------------------------------     
*
*   Write the basic parameters of the model electron cloud to the
*   .res  file; this is the second record on the file --- the
*   first is the header (see SUBROUTINE SETRES)
*
      WRITE (imcdf) NELEC, NCF, NW, nblock ! ncf is ncftot, from setcsll
*
*   Write the nuclear parameters and $- r \times V_ (r)$
*   to the  .res  file; these are the third, fourth, and fifth
*   records on the file
*
      WRITE (imcdf) Z, EMN
      WRITE (imcdf) NPARM,(PARM(I), I = 1, NPARM)
      WRITE (imcdf) N, (ZZ(I), I = 1, N), NNUC
*
*   Write the physical effects specification to the  .res  file.
*   iccutblk() is now an array of length nblock.
*   This is the sixth record on the file
*
!ASIMINA ----------------------------------------------------------------      
      WRITE (imcdf) C, LFORDR, (ICCUTBLK2(i,j), i = 1, nblock),  ! fix j
     &     LTRANS, WFACT, LVP, LNMS, LSMS
!------------------------------------------------------------------------      
*
*   Write the grid data to the  .res  file; this is the seventh
*   record on the file
*
      NP10 = N + 10
      WRITE (imcdf) RNT, H, HP, (R(I), I = 1, NP10),
     &           (RP(I), I = 1, NP10), (RPOR(I), I = 1, NP10) ! (imcdf)
*
*   Write out the interpolated radial wavefunctions; there are
*   2*NW such records; thus the total number of records written
*   at this point is 7+2*NW
*
      DO J = 1, NW
         WRITE (imcdf) E(J), GAMA(J), PZ(J), MF(J)
         WRITE (imcdf) (PF(I,J), I = 1, MF(J)), (QF(I,J), I = 1, MF(J))
      ENDDO
*
      RETURN
      END
