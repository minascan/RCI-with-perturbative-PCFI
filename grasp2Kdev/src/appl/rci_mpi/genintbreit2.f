************************************************************************
*                                                                      *
      SUBROUTINE genintbreit2 (myid, nprocs, NB, j2max)
*                                                                      *
*  Input:
*     myid, nprocs
*  Output:
*   N - Number of integrals
*   j2max - max of 2*J
*
*       Generate the list of Breit Integrals of type 2                 *
*       that could arise from a set of orbitals                        *
*       of orbitals.                                                   *
*                                                                      *
*       This routine is similar to genintrk                            *
*                                                                      *
*     Written by Per Jonsson            October 2014                   *
*     Modified for ifort -i8 by A. Kramida 22 Mar 2016            *
*                                                                      *
************************************************************************

      IMPLICIT REAL*8          (A-H, O-Z)

      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      PARAMETER (KMAX = 20)

      LOGICAL GEN,TRIANGBREIT2
      POINTER(PNTRIQ,RIQDUMMY(1))
      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
      COMMON/BILST/PINDT1,PINDT2,PINDT3,PINDT4,PINDT5,PINDT6,
     :             PVALT1,PVALT2,PVALT3,PVALT4,PVALT5,PVALT6,
     :             NDTPA(6),NTPI(6),FIRST(6)
     :      /KKSTARTBREIT/KSTARTBREIT1(0:KMAX),KSTARTBREIT2(0:KMAX)

      POINTER (PINDT1,INDTP1(1))
      POINTER (PVALT1,VALTP1(1))
      POINTER (PINDT2,INDTP2(1))
      POINTER (PVALT2,VALTP2(1))

      POINTER (PINDT3,INDTP3DUMMY(1))
      POINTER (PVALT3,VALTP3DUMMY(1))
      POINTER (PINDT4,INDTP4DUMMY(1))
      POINTER (PVALT4,VALTP4DUMMY(1))
      POINTER (PINDT5,INDTP5DUMMY(1))
      POINTER (PVALT5,VALTP5DUMMY(1))
      POINTER (PINDT6,INDTP6DUMMY(1))
      POINTER (PVALT6,VALTP6DUMMY(1))


!-----------------------------------------------------------------------
cAK Handling the -i8 option of ifort and -fdefault-integer-8 option of gfortran
cff   Setting the default integer  length to 4
      ISIZE = 4
c
      KEY = NW + 1
      KSTARTBREIT2(0) = 1
*
*   Find 2*JMAX; it should not be greater than PARAMETER KMAX
*
      J2MAX = NKJ(1)
      DO I = 2, NW
         IF (NKJ(I) .GT. J2MAX) J2MAX = NKJ(I)
      ENDDO

      IF (J2MAX .GT. KMAX) THEN
         STOP 'genintbreit2: KMAX too small'
      ENDIF
*
*   Make the breit integrals. Note that there are no symmetry relations that
*   can be utilized. See paper by Grasnt.
*   When GEN is false, sweep through to find dimension
*
      GEN = .FALSE.

  999 NB = 0
      DO K = 0, J2MAX
         DO IA = 1, NW
            DO IB = 1, NW
               DO IC = 1, NW
                  DO ID = 1, NW
                     IF (TRIANGBREIT2(IA,IB,IC,ID,K)) THEN
                        NB = NB + 1
                        IF (GEN .AND.(MOD(NB,nprocs) .EQ. myid)) THEN
                           INDTP2(NB) = ((IA*KEY+IB)*KEY+IC)*KEY+ID
                           VALTP2(NB) = BRINTF (2,IA,IB,IC,ID,K)
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         KSTARTBREIT2(K+1) = NB + 1
      ENDDO
*
*     Allocate memory for integral book keeping
*
      IF (.NOT. GEN) THEN
         CALL ALLOC (PINDT2,NB,ISIZE)
         CALL ALLOC (PVALT2,NB,8)

! Initialization is necessary in the mpi version

         DO i = 1, NB
            INDTP2(i) = 0
            VALTP2(i) = 0.d0
         ENDDO

         IF (myid .EQ. 0)
     &      PRINT *, 'Computing',NB,' Breit integrals of type 2'

         GEN = .TRUE.
         GOTO 999
      ENDIF

      RETURN
      END
