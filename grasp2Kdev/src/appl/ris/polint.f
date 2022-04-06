*******************************************************************************
*                                                                             *
      SUBROUTINE POLINT(XA,YA,DENS)
*                                                                             *
*   This routine uses interpolating polynomial d(r) = d0 + d2*r^2 + d4*r^4    *
*   to extrapolate electron density at r = 0.                                 *
*                                                                             *
*   Written by Jorgen Ekman Jul. 2016                                         *
*                                                                             *
*******************************************************************************
*      
      IMPLICIT NONE
      include 'parameters.def'
      DOUBLE PRECISION :: M(3,3), RM(3), MI(3,3), PM(3), MDET
      DOUBLE PRECISION :: XA(3), YA(3), DENS
      
      RM(1) = YA(1)
      RM(2) = YA(2)
      RM(3) = YA(3)

      M(1,1) = 1.0d0
      M(1,2) = XA(1)**2.0d0
      M(1,3) = XA(1)**4.0d0
      M(2,1) = 1.0d0
      M(2,2) = XA(2)**2.0d0
      M(2,3) = XA(2)**4.0d0
      M(3,1) = 1.d0
      M(3,2) = XA(3)**2.0d0
      M(3,3) = XA(3)**4.0d0
      MDET = M(1,1)*M(2,2)*M(3,3)-M(1,1)*M(2,3)*M(3,2)-
     :     M(1,2)*M(2,1)*M(3,3)+M(1,2)*M(2,3)*M(3,1)+
     :     M(1,3)*M(2,1)*M(3,2)-M(1,3)*M(2,2)*M(3,1)

!     DETERMINE INVERSE OF MATRIX
      MI(1,1) = (M(2,2)*M(3,3)-M(2,3)*M(3,2))/MDET
      MI(1,2) = (M(1,3)*M(3,2)-M(1,2)*M(3,3))/MDET
      MI(1,3) = (M(1,2)*M(2,3)-M(1,3)*M(2,2))/MDET

      MI(2,1) = (M(2,3)*M(3,1)-M(2,1)*M(3,3))/MDET
      MI(2,2) = (M(1,1)*M(3,3)-M(1,3)*M(3,1))/MDET
      MI(2,3) = (M(1,3)*M(2,1)-M(1,1)*M(2,3))/MDET

      MI(3,1) = (M(2,1)*M(3,2)-M(2,2)*M(3,1))/MDET
      MI(3,2) = (M(1,2)*M(3,1)-M(1,1)*M(3,2))/MDET
      MI(3,3) = (M(1,1)*M(2,2)-M(1,2)*M(2,1))/MDET

!     DETERMINE PARAMETERS FROM FIT
      PM(1) = MI(1,1)*RM(1)+MI(1,2)*RM(2)+MI(1,3)*RM(3)
      PM(2) = MI(2,1)*RM(1)+MI(2,2)*RM(2)+MI(2,3)*RM(3)
      PM(3) = MI(3,1)*RM(1)+MI(3,2)*RM(2)+MI(3,3)*RM(3)

!     Finallay RHO(0) in au^{-3} is determined
      DENS = PM(1)

      RETURN
      
      END SUBROUTINE POLINT
      
