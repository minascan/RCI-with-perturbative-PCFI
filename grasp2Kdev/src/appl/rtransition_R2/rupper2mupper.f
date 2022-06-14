*                                                                      *
      SUBROUTINE RUPPER2MUPPER
*                                                                      *
*   Input the upper integration bound for transition integrals and     *
*   determine the corresponding index number                           *
*                                                                      *
*   Written by Per Jonsson                Last revision: May 2019      *
*                                                                      *
************************************************************************
*                                                                                                                                            
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)

      include 'parameters.def'

      CHARACTER*2 NHII,NHFF
*                                                                                                                                            
      POINTER (PNTRPFII,PFII(NNNP,1)),(PNTRQFII,QFII(NNNP,1))
      POINTER (PNTRPFFF,PFFF(NNNP,1)),(PNTRQFFF,QFFF(NNNP,1))
*                                                                                                                                            
*  Initial state common                                                                                                                      
*                                                                                                                                            
      COMMON/ORB2II/NCFII,NWII
     :      /ORB4II/NPII(NNNW),NAKII(NNNW)
     :      /ORB10II/NHII(NNNW)
     :      /WAVEII/PZII(NNNW),PNTRPFII,PNTRQFII,MFII(NNNW)
*                                                                                                                                            
*  Final state common                                                                                                                        
*                                                                                                                                            
      COMMON/ORB2FF/NCFFF,NWFF
     :      /ORB4FF/NPFF(NNNW),NAKFF(NNNW)
     :      /ORB10FF/NHFF(NNNW)
     :      /WAVEFF/PZFF(NNNW),PNTRPFFF,PNTRQFFF,MFFF(NNNW)

      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /RUPPER/RUPPER0,RUPPER1,MUPPER
*                                                                                                                                            
      MTPMAX = 0
      DO I=1,NWII
        IF (MFII(I).GT.MTPMAX) THEN
           MTPMAX = MFII(I)
        END IF
      END DO
      DO I=1,NWFF
        IF (MFFF(I).GT.MTPMAX) THEN
           MTPMAX = MFFF(I)
        END IF
      END DO
*                                                                                                                                            
      WRITE(*,*) ' Max. extension of radial orbitals',R(MTPMAX),MTPMAX
      WRITE(*,*)
      WRITE(*,*) ' Give upper integration bound for transition integral'
      READ(*,*) RUPPER0
      IF (RUPPER0.GT.R(MTPMAX)) THEN
         WRITE(*,*) ' Upper bound should be lower than',R(MTPMAX)
         STOP
      END IF

      DO I = 1,NNN1
         IF (R(I).GT.RUPPER0) THEN
            MUPPER = I
            RUPPER1 = R(I)
            EXIT
         END IF
      END DO
      WRITE(*,*) 'RUPPER,MUPPER',RUPPER1,MUPPER

*     Now determine which orbitals for which the upper integration bound should be changed

      WRITE(*,*) 'Initial state orbitals'
      DO I=1,NWII
         WRITE(*,'(I2,A)') NPII(I),NHII(I)
         WRITE(*,*) 'Change yes(1), no(0)'
         READ(*,*) NUPPER
         IF (NUPPER.EQ.1) MFII(I) = MUPPER
      END DO
      WRITE(*,*) 'Final state orbitals'
      DO I=1,NWFF
         WRITE(*,'(I2,A)') NPFF(I),NHFF(I)
         WRITE(*,*) 'Change yes(1), no(0)'
         READ(*,*) NUPPER
         IF (NUPPER.EQ.1) MFFF(I) = MUPPER
      END DO

      RETURN
      END