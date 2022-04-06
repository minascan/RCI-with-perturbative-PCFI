      MODULE PCFI_PT_MOD
      INTEGER       :: N_ASIMINA, NPCFI, L
      CHARACTER*32  :: PCFINAME(10)
      INTEGER, DIMENSION(10,10) :: ICCUTBLK2
      INTEGER, DIMENSION(10) :: NTMP2
      CHARACTER*32  :: COLLECTOUT

   !   ICCUT(1) = NCSFPCFI(1) 
      
   !   DO I = 2, NPCFI
         
   !      ICCUT(I) = ICCUT(I-1) + NCSFPCFI(I)
         
   !   ENDDO
      
      END MODULE PCFI_PT_MOD
      
      
