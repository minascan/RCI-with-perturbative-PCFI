************************************************************************
*                                                                      *
      subroutine pcfiptinp 
*                                                                      *
*     It takes as input the output file from the rpcfcollect program   *
*     and gives as output the values of the:  npcfi, pcfiname(npcfi)   *
*     and iccutblk2(nblock,npcfi) variables                            *
*                                                                      * 
*     Written by Asimina Papoulia                           May 2017   * 
*                                                                      *
************************************************************************

      use pcfi_pt_mod
      
      implicit none
   
      character(200)  :: string, dummy
      character(22)   :: dummy22
      integer         :: count, dash, found, line
      integer         :: j, i
      integer         :: nblock ! as it is in the getcid subroutine ?
      integer         :: csf(10,10)
      
!     Initialize
      
!     Counter for the total number of occupied rows in the input file
      count = 0
!     Counter for the number of the "dashed" lines in the input file
      dash = 0
      
      open(20, file=trim(name)//'.pt', status='old',
     &  form='formatted') 
      
      do
         read(20, '(a)', end=99) string
         found = index(string, '-----')

         if (found.eq.1) then
            dash = dash + 1
            if (dash.eq.1) then
               line = count
            endif
         endif
         count = count + 1
      enddo
 99   continue
      rewind(20)

! Get some preliminary data
!-----------------------------------------------------------------------
!     Here we take the integer value for the npcfi variable
      npcfi = line - 3
!-----------------------------------------------------------------------
      nblock = dash/2         ! nblock is later again determined 
                              ! in the setcsll (lib/lib92)
      
! Some printing follows      
!      write(*,'(a36,i3)') ' The number of lines in the file is ', count
!      write(*,'(a28,i3)') ' The partitions end in line ', line
!      write(*,'(a35,i3)') ' The total number of dash-lines is ', dash

      write(*,'(a29,i3,a29,i3)'), ' The number of partitions is ', npcfi  
     :, ' and the number of blocks is ', nblock
            
!     Now the actual reading of the data 
      do j =1, 3
         read(20,'(a)') dummy
!         print*, trim(dummy)
      enddo
!-----------------------------------------------------------------------      
!     Here we take the character values for the NAMEPCFI(npcfi) variables
      write(*,'(a)') ' The files that contain the partitions are: '
      do j =1,npcfi
         read(20,'(a)') pcfiname(j)
         write(*,*) trim(pcfiname(j))
      enddo
!-----------------------------------------------------------------------      
      do i = 1, nblock         
         do j =1, 3      
            read(20, '(a)') dummy
!            print*, trim(dummy)
         enddo
         
         do j =1, npcfi+1
            read(20, '(a22,i12)') dummy22,csf(i,j)
!            write(*,'(a22,i12)') dummy22, csf(i,j)
         enddo         
      enddo      
!     The reading of the file has finished
      close(20)
      
!     Test that the values I need are saved
!      do i=1,nblock
!         print*, (csf(i,j), j=1,npcfi+1)
!      enddo
      
!     Now I actually get the information I need
      
      do i =1,nblock
         iccutblk2(i,1) = csf(i,1) + 1
      enddo

      do i =1,nblock
         do j =2,npcfi
            iccutblk2(i,j) = csf(i,j) + iccutblk2(i,j-1)
         enddo
      enddo
      
      do i=1,nblock
         write(*,'(a35,i3,a6)') ' The values of the iccut for block '
     :, i, ' are: '
         write(*,*) (iccutblk2(i,j), j=1,npcfi)
      enddo
            
      end 
