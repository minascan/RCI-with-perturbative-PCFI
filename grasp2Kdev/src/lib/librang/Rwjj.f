      SUBROUTINE RWJJ(J,J1,J2,K1,K2,COEF)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/RIBOJJ/IMPTJJ(63),IMGTJJ(63),IMPNJJ(63),IMGNJJ(63)
	IF(J.EQ.1) THEN
	  CALL RMEW1JJ(J1,J2,K1,K2,COEF)
	ELSEIF(J.EQ.3) THEN
          CALL RMEW3JJ(J1,J2,K1,K2,COEF)
	ELSEIF(J.EQ.5) THEN
          CALL RMEW5JJ(J1,J2,K1,K2,COEF)
	ELSEIF(J.EQ.7) THEN
          CALL RMEW7JJ(J1,J2,K1,K2,COEF)
CGG        ELSEIF(J.EQ.9) THEN
CGG          CALL SUWJJ(K1,K2,J,J1,J2,COEF)
	ELSE
	  WRITE(0,'(A,I5)') ' KLAIDA SUB. RWJJ J=',J
	  STOP
	ENDIF
      RETURN
      END
