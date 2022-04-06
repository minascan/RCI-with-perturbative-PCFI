************************************************************************
*                                                                      *
      FUNCTION RINTDENS (I,J)
*                                                                      *
*   The value of RINTDENS is an approximation to:                      *
*                                                                      *
*                                                                      *
*      (4pi^)-1  r^-2 |P (r)*P (r) + Q (r)*Q (r) | r -> 0              *
*                     I     J       I     J                            *
*                                                                      *
*   Call(s) to: [SMS92]:  POLINT                                       *
*                                                                      *
*   Written by Per Jonsson             Last revision: 24 Dec 1992      *
*   Rewritten by Jorgen Ekman                         17 Jul 2016      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = 600)
CGG      PARAMETER (NNNW = 120)
*
      DIMENSION XA(3),YA(3)
      POINTER (PNTRPF,PF(NNNP,1)),(PNTRQF,QF(NNNP,1))
*
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :     /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :     /DEF9/CVAC,PI
      
      NMIN = 10
      
      II = 1
      DO 1 L = NMIN,NMIN+4,2
         XA(II) = R(L)
         YA(II) = (PF(L,I)*PF(L,J)+QF(L,I)*QF(L,J))/
     :        (4.0D 00*PI*R(L)*R(L))
         II = II + 1
    1 CONTINUE
      
      CALL POLINT(XA,YA,DENS)
      RINTDENS = DENS
*
      RETURN
      END
