************************************************************************
*                                                                      *
      SUBROUTINE STRSUM
*                                                                      *
*   Generates the first part of  grasp92.sum  (on stream 24).          *
*                                                                      *
*   Call(s) to: [LIB92] CALEN, convrt2.                                *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 09 Dec 1992   *
*   Modified by Xinghong He               Last revision: 22 Dec 1997   *
*                                                                      *
************************************************************************
*
!     ASIMINA
      USE PCFI_PT_MOD
      
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
      CHARACTER*256 RECORD
      CHARACTER*26 CDATA
      CHARACTER*2 CLEVEL,NH
*
      POINTER (PNIVEC,IVEC(1))
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      COMMON/DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPAR/PARM(2),NPARM
     :      /NPOT/ZZ(NNNP),NNUC
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /WFAC/WFACT

      COMMON/iounit/istdi,istdo,istde

      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      POINTER (piccutblk, iccutblk(1))
      COMMON/iccu/piccutblk

*
*   Get the date and time of day; make this information the
*   header of the summary file
*
*
*   Write out the basic dimensions of the electron cloud
*
      WRITE (24,*)
      CALL convrt2 (NELEC,RECORD,LENTH, 'strsum.nelec')
      WRITE (24,*) 'There are '//RECORD(1:LENTH)
     :           //' electrons in the cloud'
      CALL convrt2 (NCF,RECORD,LENTH, 'strsum.ncf')
      WRITE (24,*) ' in '//RECORD(1:LENTH)
     :           //' relativistic CSFs'
      CALL convrt2 (NW,RECORD,LENTH, 'strsum.nw')
      WRITE (24,*) ' based on '//RECORD(1:LENTH)
     :           //' relativistic subshells.'
*
*   If the CSFs are not treated uniformly, write out an
*   informative message
*
!ASIMINA
      PRINT *, 'Hello 1'
      
      IF (LFORDR) THEN
         WRITE (24,*)
         do nb = 1, nblock
            do npp = 1, npcfi
               iccut = ICCUTBLK2(nb,npp)
!ASIMINA ----some priniting----
               PRINT *, 'nb ', nb, 'npp ', npp, 'iccut ', iccut
               
               CALL convrt2 (ICCUT,RECORD,LENTH, 'strsum.icccut')
               WRITE (24,*) ' CSFs 1--'//RECORD(1:LENTH)//' constitute'
     :              //' the zero-order space;  nb = ', nb,
     &              ' ncf = ', ncfblk(nb)
            enddo
         enddo
      ENDIF
!     ASIMINA
      PRINT *, 'Hello 2'
      
*
*   Write out the nuclear parameters
*
      WRITE (24,*)
      WRITE (24,300) Z
      IF (EMN .EQ. 0.0D 00) THEN
         WRITE (24,*) ' the nucleus is stationary;'
      ELSE
         WRITE (24,301) EMN
      ENDIF
      IF (NPARM .EQ. 2) THEN
         WRITE (24,*) ' Fermi nucleus:'
         WRITE (24,302) PARM(1),PARM(2)
         CALL convrt2 (NNUC,RECORD,LENTH, 'strsum.nnuc')
         WRITE (24,*) ' there are '//RECORD(1:LENTH)
     :              //' tabulation points in the nucleus.'
      ELSE
         WRITE (24,*) ' point nucleus.'
      ENDIF
!     ASIMINA
      PRINT *, 'Hello 3'
*
*   Write out the physical effects specifications
*
      WRITE (24,*)
      WRITE (24,303) C
*
      WRITE (24,*)
      IF (LTRANS .OR. LVP .OR. LNMS .OR. LSMS) THEN
         WRITE (24,*) 'To H (Dirac Coulomb) is added'
         IF (LTRANS) WRITE (24,304) WFACT
         IF (LVP) WRITE (24,*) ' H (Vacuum Polarisation);'
         IF (LNMS) WRITE (24,*) ' H (Normal Mass Shift);'
         IF (LSMS) WRITE (24,*) ' H (Specific Mass Shift);'
         WRITE (24,*) ' the total will be diagonalised.'
      ELSE
         WRITE (24,*) 'H (Dirac Coulomb) will be'
     :             //' diagonalised by itself.'
      ENDIF
*
      IF (LSE) THEN
         WRITE (24,*) 'Diagonal contributions from'
     :              //' H (Self Energy) will be estimated'
         WRITE (24,*) ' from a screened hydrogenic'
     :              //' approximation.'
      ENDIF
!     ASIMINA
      PRINT *, 'Hello 4'      
*
*   Write out the parameters of the radial grid
*
      WRITE (24,*)
      IF (HP .EQ. 0.0D 00) THEN
         WRITE (24,305) RNT,H,N
      ELSE
         WRITE (24,306) RNT,H,HP,N
      ENDIF
      WRITE (24,307) R(1),R(2),R(N)
!     ASIMINA
      PRINT *, 'Hello 5'     
*
*   Write out the orbital properties
*
!      WRITE (24,*)
!      WRITE (24,*) 'Subshell radial wavefunction summary:'
!      WRITE (24,*)
!      WRITE (24,308)
!      WRITE (24,*)
!      DO 1 I = 1,NW
!         WRITE (24,309) NP(I),NH(I),E(I),PZ(I),
!     :                  GAMA(I),PF(2,I),QF(2,I),MF(I)
!    1 CONTINUE
!     ASIMINA
      PRINT *, 'Hello 6'      
*
*   Write the list of eigenpair indices
*
      WRITE (24,*)
      IF (NVEC .EQ. 1) THEN
         CALL convrt2 (IVEC(1),RECORD,LENTH, 'strsum.ivec(1)')
         WRITE (24,*) 'Level '//RECORD(1:LENTH)
     :              //' will be computed.'
      ELSE
         CALL convrt2 (NVEC,RECORD,LENTH, 'strsum.nvec')
         WRITE (24,*) RECORD(1:LENTH)//' levels will be computed;'
         RECORD (1:20) = ' their indices are: '
         IEND = 20
         DO 2 I = 1,NVEC
            IBEG = IEND+1
            CALL convrt2 (IVEC(I),CLEVEL,LENTH, 'strsum.ivec(i)')
            IF (I .NE. NVEC) THEN
               IEND = IBEG+LENTH+1
               RECORD(IBEG:IEND) = CLEVEL(1:LENTH)//', '
            ELSE
               IEND = IBEG+LENTH
               RECORD(IBEG:IEND) = CLEVEL(1:LENTH)//'.'
            ENDIF
            IF (IEND .GE. 120) THEN
               WRITE (24,*) RECORD(1:IEND)
               RECORD(1:2) = '  '
               IEND = 2
            ENDIF
    2    CONTINUE
         IF (IEND .NE. 2) WRITE (24,*) RECORD(1:IEND)
      ENDIF
!     ASIMINA
      PRINT *, 'Hello 7'      
*
      RETURN
*
  300 FORMAT ('The atomic number is ',1F14.10,';')
  301 FORMAT (' the mass of the nucleus is ',1PD19.12,
     :        ' electron masses;')
  302 FORMAT ('  c =',1P,1D19.12,' Bohr radii,'
     :       /'  a =',   1D19.12,' Bohr radii;')
  303 FORMAT ('Speed of light = ',1PD19.12,' atomic units.')
  304 FORMAT ('  H (Transverse) --- factor multiplying the',
     :        ' photon frequency: ',1PD15.8,';')
  305 FORMAT ( 'Radial grid: R(I) = RNT*(exp((I-1)*H)-1),',
     :         ' I = 1, ..., N;'
     :       //' RNT  = ',1P,D19.12,' Bohr radii;'
     :        /' H    = ',   D19.12,' Bohr radii;'
     :        /' N    = ',1I4,';')
  306 FORMAT ( 'Radial grid: ln(R(I)/RNT+1)+(H/HP)*R(I) = (I-1)*H,',
     :         ' I = 1, ..., N;'
     :       //' RNT  = ',1P,D19.12,' Bohr radii;'
     :        /' H    = ',   D19.12,' Bohr radii;'
     :        /' HP   = ',   D19.12,' Bohr radii;'
     :        /' N    = ',1I4,';')
  307 FORMAT ( ' R(1) = ',1P,1D19.12,' Bohr radii;'
     :        /' R(2) = ',   1D19.12,' Bohr radii;'
     :        /' R(N) = ',   1D19.12,' Bohr radii.')
  308 FORMAT ('Subshell',6X,'e',13X,'p0',5X,
     :        'gamma',5X,'P(2)',7X,'Q(2)',4X,'MTP')
  309 FORMAT (1X,I2,A2,1X,1P,D17.10,1P,D11.3,0P,F6.2,1P,2(D11.3),I5)
*
      END
