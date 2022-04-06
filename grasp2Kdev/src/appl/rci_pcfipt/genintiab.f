************************************************************************
*                                                                      *
      SUBROUTINE GENINTIAB (myid, nprocs, N)
*                                                                      *
*  Input:                                                              *
*     myid, nprocs                                                     *
*  Output:                                                             * 
*   N - Number of integrals                                            * 
*                                                                      *
*     Generate the list of Iab Integrals that could arise from a set   *
*     of orbitals.                                                     *
*                                                                      *
*     Call(s) to: [LIB92]: ALLOC, RINTI                                *
*                                                                      *
*     Written by Asimina Papoulia                                      *
*                                                                      *
************************************************************************
      
!      USE test_mod
      IMPLICIT REAL*8          (A-H, O-Z)

      include 'parameters.def'
CGG   PARAMETER (NNNW = 120)      
CGG   PARAMETER (KEYORB = 121)

      LOGICAL GEN
      POINTER(PNTRIQ,RIQDUMMY(1))
      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
! The following were changed accordingly to the IABINT subroutine
     :      /COEILS/NDCOEA,NCOEI,PCOEIL,PCOEVL   
! Pointer arrays for the I ab integrals 
      POINTER (PCOEIL,INDOEI(1))
      POINTER (PCOEVL,VALOEI(1))

!-----------------------------------------------------------------------
      WRITE (*,*) 'Hello, Integrals I are calculated here'
      WRITE(*,*) 'Hello, N_ASIMINA value ', n_ASIMINA
      WRITE(*,*) 'NCOEI ', NCOEI
      KEY = NW + 1
      WRITE (*,*) 'NW ', NW
      DO i = 1, NW
         WRITE(*,*) 'NKJ ', NKJ(i)
      ENDDO
!      NCOEI = 0
     
! Make the I integrals
! When GEN is false, sweep through to find dimension
      GEN = .FALSE.      
      
  999 N = 0     
      DO IA = 1, NW
         DO IB = IA,NW
            IF (NAK(IA) .EQ. NAK(IB)) THEN
               N = N + 1
               NCOEI = N
            ENDIF            
            IF (GEN .AND. (MOD(N,nprocs) .EQ. myid)) THEN               
               IF (NAK(IA) .EQ. NAK(IB)) THEN               
                  INDOEI(N) = IA*KEY+IB
                  VALOEI(N) = RINTI (IA,IB,0)
                  WRITE (*,*) 'IA ' , IA, 'IB ', IB, 'N ', N,
     :              'INDOEI ', INDOEI(N), 'VALOEI ', VALOEI(N)
               ENDIF
            ENDIF           
         ENDDO
      ENDDO
            
! Allocate memory for integral book keeping      
      IF (.NOT. GEN) THEN
         CALL ALLOC (PCOEIL,N,4)
         CALL ALLOC (PCOEVL,N,8)

! Initialization is necessary in the mpi version
         DO i = 1, N
            INDOEI(i) = 0
            VALOEI(i) = 0.d0
         ENDDO

         IF (myid .EQ. 0) 
     &      PRINT *, 'Computing ',N,' I integrals'

         GEN = .TRUE.
         GOTO 999
      ENDIF
      
! now we go back to the loop and actually calculate the integrals 
      RETURN
      END
      
