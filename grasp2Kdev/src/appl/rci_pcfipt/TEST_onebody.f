*   Accumulate the contribution from the one-body operators:
*   kinetic energy, electron-nucleus interaction; update the
*     angular integral counter


         NVCOEF = 0
*
         CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
*
         DO 7 I = 1, NVCOEF
            VCOEFF = COEFF(I)
            IF (DABS (VCOEFF) .GT. CUTOFF) THEN
               NCTEC = NCTEC + 1
               IF (LSMS) THEN
                  IF (LABEL(5,I) .EQ. 1) THEN
                     CALL VINT (LABEL(1,I), LABEL(3,I), TGRL1)
                     CALL VINT (LABEL(2,I), LABEL(4,I), TGRL2)
                     ELEMNT = ELEMNT - TGRL1*TGRL2*ATWINV*VCOEFF
                  ENDIF
               ENDIF
               CALL RKINTC (LABEL(1,I), LABEL(2,I),
     :                      LABEL(3,I), LABEL(4,I),
     :                      LABEL(5,I), TEGRAL)
               ELEMNT = ELEMNT + TEGRAL*VCOEFF
            ENDIF
    7    CONTINUE

!ASIMINA
         write(*,*) 'ic ', ic, 'ir ', ir, 'ELEMNT ', ELEMNT
*
         IBUG1 = 0

      ENDIF                     !inc1 is always 1 without PRE-RUN



      IF (IA .NE. 0) THEN
         DO 13 I = 1, NTCOEF
            TCOEFF = DBLE(TSHELL(IA))
            IF (DABS (TCOEFF) .GT. CUTOFF) THEN
               NCOEC = NCOEC + 1
               IF (LNMS) THEN
                  CALL KEINT (IA, IB, TEGRAL)
                  !------------------------
                  ELEMNT = ELEMNT + TEGRAL*ATWINV*TCOEFF
               ENDIF
               IF (LVP) THEN
                  CALL VPINT (IA, IB, TEGRAL)
                  !------------------------
                  ELEMNT = ELEMNT + TEGRAL*TCOEFF
               ENDIF
               CALL IABINTC (LABEL(1,I), LABEL(2,I)
     :              , TEGRAL)
              ELEMNT = ELEMNT + TEGRAL*TCOEFF 
           ENDIF
 13     CONTINUE
      ENDIF
      

        
      
