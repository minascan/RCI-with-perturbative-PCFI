!***********************************************************************
!                                                                      *
      PROGRAM uniquelbl
!-----------------------------------------------
!                                                                      *
!   Program defining unique labels for energy levels                   * 
!                                                                      *
!   Written by  G. Gaigalas                         Vilnius, May 2016  *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!     The maximum number of levels in the list
      INTEGER, PARAMETER :: Lev_No  = 100
!     The maximum number of mixing coefficients of ASF expension
      INTEGER, PARAMETER :: Vec_No  = 100
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER, DIMENSION(Vec_No,Lev_No)    :: COUPLING*256
      REAL(DOUBLE), DIMENSION(Vec_No,Lev_No) :: COMP, MIX
      CHARACTER, DIMENSION(Lev_No)           :: Str_No*3,Str_J*5,Str_P*1
      CHARACTER, DIMENSION(Lev_No)           :: OPT_COUPLING*256
      CHARACTER, DIMENSION(Vec_No)           :: tmp_COUPLING*256
      REAL(DOUBLE), DIMENSION(Lev_No)        :: ENERGY, PRO
      REAL(DOUBLE), DIMENSION(Vec_No)        :: tmp_COMP, tmp_MIX
      INTEGER, DIMENSION(Lev_No)             :: MAS_MAX, ICOUNT, IPRGG
      CHARACTER    :: RECORD*53, FILNAM*256
      CHARACTER    :: FORM*11, STATUS*3
      CHARACTER    :: Str*7
      REAL(DOUBLE) :: MAX_COMP
      INTEGER      :: IOS, timecount, IERR, INUM, NUM_Lev, I1,I2,I3,I4
      INTEGER      :: IPR, NUM_OPT, Lev_OPT, New_J
      INTEGER*4    :: last, length
!-----------------------------------------------
      call starttime (timecount, 'uniquelbl')
!
      print *, ""
      print *, "uniquelbl: Defines  unique labels  for energy levels."
      print *, "           (C) Copyright by G. Gaigalas and P. Rynkun"
      print *, "           (Fortran 95 version)      Vilnius  (2016)."
      print *, "           Input files: **.lbl"
      print *, "           Output file: *.uni.*.lbl, *.uni.*.sum"
      print *, ""
!
!    Loading inpit files
!
      print *,  "Give the name of the input file"
      read(*,'(a)') FILNAM
      write (6, *) 'Loading the file ',trim(FILNAM)//'.lbl'
      last = index(trim(FILNAM),'.',.true.)
      length = len(trim(FILNAM))
      form = 'FORMATTED'
      status = 'OLD'
      call OPENFL (21,trim(FILNAM)//'.lbl',FORM,STATUS,IERR)
      if (IERR == 1) then
         write (6, *) 'Error when opening ',trim(FILNAM)//'.lbl'
         stop
      end if
!
      read (21, '(1A53)', IOSTAT=IOS) RECORD
      if (IOS/=0) then
         write (6, *) 'Not a ',trim(FILNAM)//'.lbl',' File;'
         close(21)
      endif
!
      status = 'UNKNOWN'
      call OPENFL (22,                                                 &
                    FILNAM(1:last)//'uni'//FILNAM(last:length)//'.sum',&
                                                       FORM,STATUS,IERR)
      call OPENFL (23,                                                 &
                    FILNAM(1:last)//'uni'//FILNAM(last:length)//'.lbl',&
                                                       FORM,STATUS,IERR)
      write(23, '(1A53)') RECORD
!
      print*,""
      DO WHILE (New_J < 1)
        NUM_Lev = 0
        New_J = 1
        DO
          if(NUM_Lev > Lev_No) then
            print*, "Please extand the arrays. Now we have Lev_No=",   &
                                                                  Lev_No
            stop
          end if
          NUM_Lev = NUM_Lev + 1
          read(21,'(A3,A5,5X,A1,8X,F16.9,5X,F7.3)')Str_No(NUM_Lev),    &
              Str_J(NUM_Lev),Str_P(NUM_Lev),Energy(NUM_Lev),PRO(NUM_Lev)
          INUM = 0
          DO
            INUM = INUM + 1
            if(INUM > Vec_No) then
              print*, "Please extand the arrays. Now we have Vec_No=", &
                                                                  Vec_No
              stop
              end if
            read(21,'(A7,F12.8,3X,F11.8,3X,A)',IOSTAT=IOS)Str,         &
             MIX(INUM,NUM_Lev),COMP(INUM,NUM_Lev),COUPLING(INUM,NUM_Lev)
            if (IOS==-1) then
               MAS_MAX(NUM_Lev) = INUM - 1
               go to 1
            endif
            if(Str /= '       ') then
               backspace(21)
               MAS_MAX(NUM_Lev) = INUM - 1
             exit
            else if                                                    &
              (MIX(INUM,NUM_Lev)==0.00.and.COMP(INUM,NUM_Lev)==0.00)then
               MAS_MAX(NUM_Lev) = INUM - 1
               New_J = 0
               go to 1
            end if
          END DO
        END DO
   1    CONTINUE
!
        NUM_OPT = 0
        ICOUNT = 1
        IPRGG = 1
        write(22,'(A)')                                                &
                    "          Composition  Serial No.         Coupling"
        write(22,'(A)') "                       of compos."
        write(22,'(A5,A)') " J = ",trim(Str_J(NUM_Lev))
        write(22,'(A)')                                                &
                    "--------------------------------------------------"
        DO I1 = 1, NUM_Lev
   2      MAX_COMP = 0.0
          DO I2 = 1, NUM_Lev
            if(ICOUNT(I2) == 0) CYCLE
            if(COMP(1,I2) > MAX_COMP) then
               Lev_OPT = I2
               MAX_COMP = COMP(1,I2)
            end if
          END DO
          NUM_OPT = NUM_OPT + 1
          ICOUNT(Lev_OPT) = 0
          IPR = 1
          if(NUM_OPT == 1) then
            IPR = 1
          else
            I3 = 1
            DO WHILE (I3 < NUM_OPT)
             if(trim(OPT_COUPLING(I3))==trim(COUPLING(IPR,Lev_OPT)))then
               IPR = IPR + 1
               if(IPR > Vec_No) then
                 print*,"Please extand the arrays. Now we have Vec_No="&
                                                                 ,Vec_No
                 stop
               end if
               I3 = 1
             else
               I3 = I3 + 1
             end if
            END DO
          end if
          IPRGG(Lev_OPT) = IPR + IPRGG(Lev_OPT) - 1
          if(IPRGG(Lev_OPT) >= MAS_MAX(Lev_OPT)) then
             print*,                                                   &
             "The program is not able perform the identification for", &
             " level = ",Lev_OPT
             stop
          end if
          if(IPR == 1) then
            OPT_COUPLING(NUM_OPT) = COUPLING(IPR,Lev_OPT)
            write(22,'(A,I4,2X,F12.9,I5,3X,A,A,I4)')"Pos",Lev_OPT,     &
            COMP(IPR,Lev_OPT),IPRGG(Lev_OPT),trim(OPT_COUPLING(NUM_OPT))
          end if
          if(IPR > 1) then
            tmp_MIX = 0
            tmp_COMP = 0
            tmp_COUPLING = ""
            I4 = MAS_MAX(Lev_OPT)
            tmp_MIX(1:1) = MIX(IPR:IPR,Lev_OPT)
            tmp_MIX(I4-IPR+2:I4) = MIX(1:IPR-1,Lev_OPT)
            tmp_MIX(2:I4-IPR+1) = MIX(IPR+1:I4,Lev_OPT)
            MIX(1:I4,Lev_OPT) = tmp_MIX(1:I4)
!
            tmp_COMP(1:1) = COMP(IPR:IPR,Lev_OPT)
            tmp_COMP(I4-IPR+2:I4) = COMP(1:IPR-1,Lev_OPT)
            tmp_COMP(2:I4-IPR+1) = COMP(IPR+1:I4,Lev_OPT)
            COMP(1:I4,Lev_OPT) = tmp_COMP(1:I4)
!
            tmp_COUPLING(1:1) = COUPLING(IPR:IPR,Lev_OPT)
            tmp_COUPLING(I4-IPR+2:I4) = COUPLING(1:IPR-1,Lev_OPT)
            tmp_COUPLING(2:I4-IPR+1) = COUPLING(IPR+1:I4,Lev_OPT)
            COUPLING(1:I4,Lev_OPT) = tmp_COUPLING(1:I4)
!
            ICOUNT(Lev_OPT) = 1
            NUM_OPT = NUM_OPT - 1
            go to 2
          end if
        END DO
        write(22,'(A)')                                                &
                    "--------------------------------------------------"
        write(22,'(A)')""
!
        DO I1 = 1, NUM_Lev
          write(23,'(A3,A5,5X,A1,8X,F16.9,5X,F7.3,A)')                 &
          Str_No(I1),Str_J(I1),Str_P(I1),Energy(I1),PRO(I1),"%"
          DO I2 = 1, MAS_MAX(I1)
            write(23,'(A7,F12.8,3X,F11.8,3X,A)')                       &
            Str,MIX(I2,I1),COMP(I2,I1),trim(COUPLING(I2,I1))
          END DO
        END DO
        if(New_J < 1) write(23,'(A)') " "
!
      end do
      call stoptime (timecount, 'uniquelbl')
      STOP  
      END PROGRAM uniquelbl
