************************************************************************
*                                                                      *
      PROGRAM RIS
*                                                                      *
*   Entry routine for RIS. Controls the entire computation.            *
*                                                                      *
*   Call(s) to: [LIB92]: GETMIX, SETCSL, SETMC, SETCON.                *
*               [SMS92]: CHKPLT, GETSMD, SETDBG, SETSUM                *
*                        STRSUM, SETDENS.                              *
*               [NJGRAF]: FACTT.                                       *
*                                                                      *
*   Written by Per Jonsson                Last revision: 17 Jan 1996   *
*   Modify  by Gediminas Gaigalas                        26 Oct 2009   *
*   Modified by J. Ekman                                 28 Mar 2016   *
*                                                                      *
************************************************************************
*
      DOUBLE PRECISION DR2
      LOGICAL GETYN, YES
      CHARACTER*24 NAME
      COMMON/DEFAULT/NDEF

      PRINT *
      PRINT *, 'RIS: Execution begins ...'
      CALL STARTTIME (ncount1, 'RIS')
      PRINT *
      PRINT *, 'Default settings?'
      YES = GETYN ()
      PRINT *
      IF (YES) THEN
         NDEF = 0
      ELSE
         NDEF = 1
      ENDIF

    9 PRINT *, 'Name of state'
      READ(*,'(A)') NAME
      K=INDEX(NAME,' ')
      IF (K.EQ.1) THEN
         PRINT *, 'Names may not start with a blank'
         GOTO 9
      ENDIF
      PRINT *

      PRINT *, 'Mixing coefficients from a CI calc.?'
      YES = GETYN ()
      IF (YES) THEN
         NCI = 0
      ELSE
         NCI = 1
      ENDIF
      PRINT *

*
*   Check compatibility of plant substitutions
*
      CALL CHKPLT
*
*   Determine if there is to be any debug printout; this will be
*   made on the  .dbg  file
*
      CALL SETDBG
*
*   Perform machine- and installation-dependent setup
*
      CALL SETMC
*
*   Set up the physical constants
*
      CALL SETCON
*
*   Open the  .i  file
*
      CALL SETSUM(NAME,NCI)
*
*   Open, check, load data from, and close, the  .csl  file
*
      CALL SETCSLA(NAME,ncore_not_used)
*
*   Get the remaining information
*
      CALL GETSMD(NAME)
*
*   Get the eigenvectors
*
*      PRINT *, 'Block format?'
*      YES = GETYN ()
*      PRINT *

*      IF (YES) THEN
         CALL GETMIXBLOCK(NAME,NCI)
*      ELSE
*         IF (NCI.EQ.0) THEN
*            CALL GETMIXC(NAME)
*         ELSE
*            CALL GETMIXA(NAME)
*         ENDIF
*      ENDIF
*
*   Append a summary of the inputs to the  .sum  file
*
*      CALL STRSUM
*
*   Set up the table of logarithms of factorials
*
      CALL FACTT
*
*   Proceed with the RIS calculation
*
      CALL RIS_CAL(NAME)
*
*   Print completion message
*
      PRINT *
      CALL STOPTIME (ncount1, 'RIS') 
*      PRINT *, 'RIS: Execution complete.'
*
*      STOP
      END
