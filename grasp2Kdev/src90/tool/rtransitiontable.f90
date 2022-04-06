program rtransitiontable

! This program makes ASCII and LaTeX tables over transition data

! Per Jonssson, Malmo University, August 2014

implicit none
integer :: h,i,j,k,l, nfile, ntransition, ncase, maxlengthascii1, maxlengthlatex1
integer :: maxlengthascii2, maxlengthlatex2, nt, nformat, nskip, lastpos
character(len=1) :: char1, char2, char3
character(len=2) :: j1,j2
character(len=100) :: filename
character(len=200) :: line1, line2, line3, line4, line5, line6, line7, linedummy
character(len=180) :: labelstring1, labelstring2, dummystring
character(len=180) :: latexstring1(10000), latexstring2(10000),asciistring1(10000),asciistring2(10000)
!character(len=11) :: energystring(10000),sstring(10000), gfstring(10000), Astring(10000)
character(len=9) :: sstring(10000), gfstring(10000), Astring(10000)
character(len=1) :: extra1(10000), extra2(10000)
character(len=8) :: energystring(1000)
character(len=13) :: wavelengthstring(10000)
!character(len=9) :: dTstring(10000)
character(len=7) :: dTstring(10000)

open(unit=19,file='transitiontablelatex.tex',status='unknown')
open(unit=20,file='transitiontableascii.txt',status='unknown')

write(*,*) 
write(*,*) ' RTRANSITIONTABLE'
write(*,*) ' Makes LaTeX tables of transition data from transition files'
write(*,*) ' name1.name2.ct.lsj '
write(*,*) ' Input file: name1.name2.ct.lsj'
write(*,*) ' Output file: transitiontablelatex.tex'
write(*,*) 

write(*,*) ' Specify tabel format '
write(*,*) ' (1). Lower & Upper & Energy diff. & wavelength & S & gf & A & dT '   
write(*,*) ' (2). Lower & Upper & Energy diff. & wavelength & gf & A & dT '   
write(*,*) ' (3). Lower & Upper & Energy diff. & wavelength & gf & A '   
write(*,*) ' (4). Lower & Upper & Energy diff. & S & gf & A & dT '   
write(*,*) ' (5). Lower & Upper & Energy diff. & gf & A & dT '   
write(*,*) ' (6). Lower & Upper & Energy diff. & gf & A '   
read(*,*) nformat

write(*,*) ' Inspect the name1.name2.ct.lsj file and determine how many positions'
write(*,*) ' should be skipped in the string that determines the label '
write(*,*) ' e.g. if the string is 1s(2).2s_2S.2p(2)3P2_4P and 1s(2) is a core'
write(*,*) ' then you would like to skip 1s(2). i.e. 6 positions and determine'
write(*,*) ' the label from 2s_2S.2p(2)3P2_4P'
write(*,*) 
write(*,*) ' How many positions should be skipped?'
read(*,*) nskip

!write(*,*) 'Give the number of files'
!read(*,*) nfile

write(19,'(a)') '\documentclass[10pt]{article}'
write(19,'(a)') '\usepackage{longtable}'
write(19,'(a)') '\begin{document}'


if (nformat.eq.1) then
   write(19,'(a)') '\begin{longtable}{llrrrrrr}'
   write(19,'(a)') ' Lower state & Upper state & $\Delta E$ & $\lambda$ & $S$ & $gf$ & $A$ & $dT$ \\ \hline'
elseif (nformat.eq.2) then
   write(19,'(a)') '\begin{longtable}{llrrrrr}'
   write(19,'(a)') ' Lower state & Upper state & $\Delta E$ & $\lambda$ & $S$ & $gf$ & $A$ \\ \hline'
elseif (nformat.eq.3) then
   write(19,'(a)') '\begin{longtable}{llrrrr}'
   write(19,'(a)') ' Lower state & Upper state & $\Delta E$ & $\lambda$ & $gf$ & $A$ \\ \hline'
elseif (nformat.eq.4) then
   write(19,'(a)') '\begin{longtable}{llrrrrr}'
   write(19,'(a)') ' Lower state & Upper state & $\Delta E$ & $S$ & $gf$ & $A$ & $dT$ \\ \hline'
elseif (nformat.eq.5) then
   write(19,'(a)') '\begin{longtable}{llrrrr}'
   write(19,'(a)') ' Lower state & Upper state & $\Delta E$ & $gf$ & $A$ & $dT$ \\ \hline'
elseif (nformat.eq.6) then
   write(19,'(a)') '\begin{longtable}{llrrr}'
   write(19,'(a)') ' Lower state & Upper state & $\Delta E$ & $gf$ & $A$  \\ \hline'
end if

nfile = 1
do h = 1,nfile



   write(*,'(a,i2)') 'Name of file',h
   read(*,'(a)') filename

   open(unit=20+h,file=trim(filename),status='old')

! Start reading the file 

   do j = 1,3
      read(20+h,'(a)') linedummy
   end do


   nt = 0
   maxlengthascii1 = 0
   maxlengthlatex1 = 0
   maxlengthascii2 = 0
   maxlengthlatex2 = 0
   do
      read(20+h,'(a)',end=999) linedummy
      read(20+h,'(a)') linedummy
      read(20+h,'(a)') line1
      read(20+h,'(a)') line2
      read(20+h,'(a)') line3
      read(20+h,'(a)') line4
      read(20+h,'(a)') line5
      nt = nt + 1

      labelstring1 = line1(21+nskip:200)
      labelstring2 = line2(21+nskip:200)

      lastpos = len_trim(labelstring1)

      select case (labelstring1(lastpos:lastpos))
      case ('S','P','D','F','G','H','I','K','L','M','N')
         extra1(nt) = ' '
      case default
         extra1(nt) = labelstring1(lastpos:lastpos)
         labelstring1(lastpos:lastpos) = ' '
      end select

      lastpos = len_trim(labelstring2)

      select case (labelstring2(lastpos:lastpos))
      case ('S','P','D','F','G','H','I','K','L','M','N')
         extra2(nt) = ' '
      case default
         extra2(nt) = labelstring2(lastpos:lastpos)
         labelstring2(lastpos:lastpos) = ' '
      end select







      asciistring1(nt) = line1(21+nskip:200)
      asciistring2(nt) = line2(21+nskip:200)
      j1 = line1(3:4)
      j2 = line2(3:4)

!      energystring(nt) = line3(1:11)
!      wavelengthstring(nt) = line3(17:29)
!      sstring(nt) = line4(11:21)
!      gfstring(nt) = line4(31:41)
!      Astring(nt) = line4(52:62)
!      dTstring(nt) = line4(70:78)

      energystring(nt) = line3(1:8)
      wavelengthstring(nt) = line3(17:29)
      sstring(nt) = line4(11:15)//line4(18:21)
      gfstring(nt) = line4(31:35)//line4(38:41)
      Astring(nt) = line4(52:56)//line4(59:62)
      dTstring(nt) = line4(70:76)
!      write(*,*) labelstring1
!      write(*,*) labelstring2
  

! Convert quantum labels
! to LaTeX


      do i = 1,177

!  Replace (n) with ^n

         if ((labelstring1(i:i).eq.'(').and.(labelstring1(i+2:i+2).eq.')')) then
            labelstring1(i:i) = '^'
            labelstring1(i+2:i+2) = ' '
         end if
         if ((labelstring2(i:i).eq.'(').and.(labelstring2(i+2:i+2).eq.')')) then
            labelstring2(i:i) = '^'
            labelstring2(i+2:i+2) = ' '
         end if
      end do

      do i = 1,177

!  Replace . with \,

         if (labelstring1(i:i).eq.'.') then
            dummystring = labelstring1
            labelstring1(1:i-1) = dummystring(1:i-1)
            labelstring1(i:i) = '\'
            labelstring1(i+1:i+1) = ','
            labelstring1(i+2:145) = dummystring(i+1:143) 
         end if
         if (labelstring2(i:i).eq.'.') then
            dummystring = labelstring2
            labelstring2(1:i-1) = dummystring(1:i-1)
            labelstring2(i:i) = '\'
            labelstring2(i+1:i+1) = ','
            labelstring2(i+2:145) = dummystring(i+1:143) 
         end if
      end do

      do i = 1,177

!  Replace _ with ~

         if (labelstring1(i:i).eq.'_') labelstring1(i:i) = '~'
         if (labelstring2(i:i).eq.'_') labelstring2(i:i) = '~'
      end do
!      write(*,'(a)') trim(labelstring1)
!      write(*,'(a)') trim(labelstring2)

!  If integer1 and S, P, D, F, G, H, I, K, L, M, N and integer2 replace with (^integer1_integer2S), (^integer1_integer2P), etc

      do l = 1,15
         ncase = 0
         do i = 1,177
            do j = 48,57
               do k = 48,57 
                  char1 = labelstring1(i:i)
                  char2 = labelstring1(i+1:i+1)
                  char3 = labelstring1(i+2:i+2)
                  if ((ichar(char1).eq.j).and.(ichar(char3).eq.k).and.((char2.ne.'~').and.(char2.ne.' ').and.(char2.ne.'_'))) then
                     dummystring = labelstring1
                     labelstring1(1:i-1) = dummystring(1:i-1)
                     labelstring1(i:i+6) = '(^'//char1//'_'//char3//char2//')'
                     labelstring1(i+7:145) = dummystring(i+3:141)
                     ncase = ncase + 1
                  end if
               end do
            end do
            if (ncase.eq.1) exit
         end do

!         write(*,'(a)') trim(labelstring1)
      end do

      do l = 1,15
         ncase = 0
         do i = 1,177
            do j = 48,57
               do k = 48,57 
                  char1 = labelstring2(i:i)
                  char2 = labelstring2(i+1:i+1)
                  char3 = labelstring2(i+2:i+2)
                  if ((ichar(char1).eq.j).and.(ichar(char3).eq.k).and.((char2.ne.'~').and.(char2.ne.' ').and.(char2.ne.'_'))) then
                     dummystring = labelstring2
                     labelstring2(1:i-1) = dummystring(1:i-1)
                     labelstring2(i:i+6) = '(^'//char1//'_'//char3//char2//')'
                     labelstring2(i+7:145) = dummystring(i+3:141)
                     ncase = ncase + 1
                  end if
               end do
            end do
            if (ncase.eq.1) exit
         end do

!         write(*,'(a)') trim(labelstring2)
      end do

!  If integer1 and S, P, D, F, G, H, I, K, L, M, N and not integer2 replace with ^integer1S, ^integer1P, etc 

      do i = 1,177
!  
         if (labelstring1(i:i).eq.'~') then
            dummystring = labelstring1
            labelstring1(1:i) = dummystring(1:i)
            labelstring1(i+1:i+1) = '^'
            labelstring1(i+2:145) = dummystring(i+1:143)
         end if
         if (labelstring2(i:i).eq.'~') then
            dummystring = labelstring2
            labelstring2(1:i) = dummystring(1:i)
            labelstring2(i+1:i+1) = '^'
            labelstring2(i+2:145) = dummystring(i+1:143)
         end if
      end do

      if (j1.eq.' 0') labelstring1 = '$'//trim(labelstring1)//'_{0}$' 
      if (j2.eq.' 0') labelstring2 = '$'//trim(labelstring2)//'_{0}$' 
      if (j1.eq.' 1') labelstring1 = '$'//trim(labelstring1)//'_{1/2}$' 
      if (j2.eq.' 1') labelstring2 = '$'//trim(labelstring2)//'_{1/2}$' 
      if (j1.eq.' 2') labelstring1 = '$'//trim(labelstring1)//'_{1}$' 
      if (j2.eq.' 2') labelstring2 = '$'//trim(labelstring2)//'_{1}$' 
      if (j1.eq.' 3') labelstring1 = '$'//trim(labelstring1)//'_{3/2}$' 
      if (j2.eq.' 3') labelstring2 = '$'//trim(labelstring2)//'_{3/2}$' 
      if (j1.eq.' 4') labelstring1 = '$'//trim(labelstring1)//'_{2}$' 
      if (j2.eq.' 4') labelstring2 = '$'//trim(labelstring2)//'_{2}$' 
      if (j1.eq.' 5') labelstring1 = '$'//trim(labelstring1)//'_{5/2}$' 
      if (j2.eq.' 5') labelstring2 = '$'//trim(labelstring2)//'_{5/2}$' 
      if (j1.eq.' 6') labelstring1 = '$'//trim(labelstring1)//'_{3}$' 
      if (j2.eq.' 6') labelstring2 = '$'//trim(labelstring2)//'_{3}$' 
      if (j1.eq.' 7') labelstring1 = '$'//trim(labelstring1)//'_{7/2}$' 
      if (j2.eq.' 7') labelstring2 = '$'//trim(labelstring2)//'_{7/2}$' 
      if (j1.eq.' 8') labelstring1 = '$'//trim(labelstring1)//'_{4}$' 
      if (j2.eq.' 8') labelstring2 = '$'//trim(labelstring2)//'_{4}$' 
      if (j1.eq.' 9') labelstring1 = '$'//trim(labelstring1)//'_{9/2}$' 
      if (j2.eq.' 9') labelstring2 = '$'//trim(labelstring2)//'_{9/2}$' 
      if (j1.eq.'10') labelstring1 = '$'//trim(labelstring1)//'_{5}$' 
      if (j2.eq.'10') labelstring2 = '$'//trim(labelstring2)//'_{5}$' 
      if (j1.eq.'11') labelstring1 = '$'//trim(labelstring1)//'_{11/2}$' 
      if (j2.eq.'11') labelstring2 = '$'//trim(labelstring2)//'_{11/2}$' 
      if (j1.eq.'12') labelstring1 = '$'//trim(labelstring1)//'_{6}$' 
      if (j2.eq.'12') labelstring2 = '$'//trim(labelstring2)//'_{6}$' 
      if (j1.eq.'13') labelstring1 = '$'//trim(labelstring1)//'_{13/2}$' 
      if (j2.eq.'13') labelstring2 = '$'//trim(labelstring2)//'_{13/2}$' 
      if (j1.eq.'14') labelstring1 = '$'//trim(labelstring1)//'_{7}$' 
      if (j2.eq.'14') labelstring2 = '$'//trim(labelstring2)//'_{7}$' 
      if (j1.eq.'15') labelstring1 = '$'//trim(labelstring1)//'_{15/2}$' 
      if (j2.eq.'15') labelstring2 = '$'//trim(labelstring2)//'_{15/2}$' 
      if (j1.eq.'16') labelstring1 = '$'//trim(labelstring1)//'_{8}$' 
      if (j2.eq.'16') labelstring2 = '$'//trim(labelstring2)//'_{8}$' 
      if (j1.eq.'17') labelstring1 = '$'//trim(labelstring1)//'_{17/2}$' 
      if (j2.eq.'17') labelstring2 = '$'//trim(labelstring2)//'_{17/2}$' 
      if (j1.eq.'18') labelstring1 = '$'//trim(labelstring1)//'_{9}$' 
      if (j2.eq.'18') labelstring2 = '$'//trim(labelstring2)//'_{9}$' 
      if (j1.eq.'19') labelstring1 = '$'//trim(labelstring1)//'_{19/2}$' 
      if (j2.eq.'19') labelstring2 = '$'//trim(labelstring2)//'_{19/2}$' 

      if (len_trim(labelstring1).gt.maxlengthlatex1) maxlengthlatex1 = len_trim(labelstring1)   
      if (len_trim(labelstring2).gt.maxlengthlatex2) maxlengthlatex2 = len_trim(labelstring2)   

      latexstring1(nt) = labelstring1
      latexstring2(nt) = labelstring2

      if (j1.eq.' 1') asciistring1(nt) = trim(asciistring1(nt))//' 1/2' 
      if (j2.eq.' 1') asciistring2(nt) = trim(asciistring2(nt))//' 1/2' 

      if (len_trim(asciistring1(nt)).gt.maxlengthascii1) maxlengthascii1 = len_trim(asciistring1(nt))   
      if (len_trim(asciistring2(nt)).gt.maxlengthascii2) maxlengthascii2 = len_trim(asciistring2(nt))   

   end do
999 continue

   do i = 1,nt
      if (nformat.eq.1) then
      write(19,'(a)') latexstring1(i)(1:maxlengthlatex1)//'~'//extra1(i)//' & '//&
                     & latexstring2(i)(1:maxlengthlatex2)//'~'//extra2(i)//' & '//&
                     &energystring(i)//' & '//wavelengthstring(i)//' & '//sstring(i)//' & '//gfstring(i)//' & '//&
                     &Astring(i)//' & '//dTstring(i)//'\\'
      elseif (nformat.eq.2) then
      write(19,'(a)') latexstring1(i)(1:maxlengthlatex1)//'~'//extra1(i)//' & '//&
                     & latexstring2(i)(1:maxlengthlatex2)//'~'//extra2(i)//' & '//&
                     &energystring(i)//' & '//wavelengthstring(i)//' & '//gfstring(i)//' & '//&
                     &Astring(i)//' & '//dTstring(i)//'\\'
      elseif (nformat.eq.3) then
      write(19,'(a)') latexstring1(i)(1:maxlengthlatex1)//'~'//extra1(i)//' & '//&
                     & latexstring2(i)(1:maxlengthlatex2)//'~'//extra2(i)//' & '//&
                     &energystring(i)//' & '//wavelengthstring(i)//' & '//gfstring(i)//' & '//&
                     &Astring(i)//'\\'
      elseif (nformat.eq.4) then
      write(19,'(a)') latexstring1(i)(1:maxlengthlatex1)//'~'//extra1(i)//' & '//&
                     & latexstring2(i)(1:maxlengthlatex2)//'~'//extra2(i)//' & '//&
                     &energystring(i)//' & '//sstring(i)//' & '//gfstring(i)//' & '//&
                     &Astring(i)//' & '//dTstring(i)//'\\'
      elseif (nformat.eq.5) then
      write(19,'(a)') latexstring1(i)(1:maxlengthlatex1)//'~'//extra1(i)//' & '//&
                     & latexstring2(i)(1:maxlengthlatex2)//'~'//extra2(i)//' & '//&
                     &energystring(i)//' & '//gfstring(i)//' & '//&
                     &Astring(i)//' & '//dTstring(i)//'\\'
      elseif (nformat.eq.6) then
      write(19,'(a)') latexstring1(i)(1:maxlengthlatex1)//'~'//extra1(i)//' & '//&
                     & latexstring2(i)(1:maxlengthlatex2)//'~'//extra2(i)//' & '//&
                     &energystring(i)//' & '//gfstring(i)//' & '//&
                     &Astring(i)//'\\'
      end if
      if (nformat.eq.1) then
      write(20,'(a)') asciistring1(i)(1:maxlengthascii1)//'  '// asciistring2(i)(1:maxlengthascii2)//'  '//&
                     &energystring(i)//'  '//wavelengthstring(i)//'  '//sstring(i)//'  '//gfstring(i)//'  '//&
                     &Astring(i)//'  '//dTstring(i)
      elseif (nformat.eq.2) then
      write(20,'(a)') asciistring1(i)(1:maxlengthascii1)//'  '// asciistring2(i)(1:maxlengthascii2)//'  '//&
                     &energystring(i)//'  '//wavelengthstring(i)//'  '//gfstring(i)//'  '//&
                     &Astring(i)//'  '//dTstring(i)
      elseif (nformat.eq.3) then
      write(20,'(a)') asciistring1(i)(1:maxlengthascii1)//'  '// asciistring2(i)(1:maxlengthascii2)//'  '//&
                     &energystring(i)//'  '//wavelengthstring(i)//'  '//gfstring(i)//'  '//&
                     &Astring(i)
      elseif (nformat.eq.4) then
      write(20,'(a)') asciistring1(i)(1:maxlengthascii1)//' & '// asciistring2(i)(1:maxlengthascii2)//' & '//&
                     &energystring(i)//' & '//sstring(i)//' & '//gfstring(i)//' & '//&
                     &Astring(i)//' & '//dTstring(i)//'\\'
      elseif (nformat.eq.5) then
      write(20,'(a)') asciistring1(i)(1:maxlengthascii1)//' & '// asciistring2(i)(1:maxlengthascii2)//' & '//&
                     &energystring(i)//' & '//gfstring(i)//' & '//&
                     &Astring(i)//' & '//dTstring(i)//'\\'
      elseif (nformat.eq.6) then
      write(20,'(a)') asciistring1(i)(1:maxlengthascii1)//' & '// asciistring2(i)(1:maxlengthascii2)//' & '//&
                     &energystring(i)//' & '//gfstring(i)//' & '//&
                     &Astring(i)//'\\'
      end if
   end do

write(19,'(a)') '\hline\\'
write(19,'(a)') '\caption{Transition data from the file '//trim(filename)//'}'
write(19,'(a)') '\end{longtable}'
write(19,'(a)') '\end{document}'

end do

end program rtransitiontable


