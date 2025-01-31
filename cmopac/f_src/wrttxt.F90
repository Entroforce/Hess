      subroutine wrttxt(iprt)  
      use molkst_C, only : koment, title, refkey, keywrd, line
      implicit none
      integer , intent(in) :: iprt 
      integer :: i, j
      logical :: l_chains = .false., l_start = .false.
!-----------------------------------------------
!
! Is CHAINS or START_RES present in the data set?
!
      do i = 1, 6
        if (index(refkey(i), " NULL") /= 0) exit
        line = " "//trim(refkey(i))
        call upcase(line, len_trim(line))
        if (.not. l_chains) l_chains = (index(line, " CHAINS") /= 0)
        if (.not. l_start)  l_start  = (index(line, " START_RES") /= 0)
      end do
!
!  Is CHAINS present in the keyword?
!
      i = index(keywrd," CHAINS")
      if(i /= 0 .and. .not. l_chains) then
!
!  CHAINS is present in keyword, but was not present in the data-set,
!  so add CHAINS keyword to refkey(1)
!
        j = index(keywrd(i + 7:), ")") + i + 7
        refkey(1) = keywrd(i:j)//trim(refkey(1))
      end if
!
!  Is START_RES present in the keyword?
!
      i = index(keywrd," START_RES")
      if(i /= 0 .and. .not.l_start) then
!
!  START_RES is present in keyword, but was not present in the data-set,
!  so add START_RES keyword to refkey(1)
!
        j = index(keywrd(i + 10:), ")") + i + 10
        refkey(1) = keywrd(i:j)//trim(refkey(1))
      end if   
      do i = 1, 6
        if (index(refkey(i), " NULL") /= 0) exit
        write(iprt,'(a)')trim(refkey(i))
      end do
      if (index(koment,' NULL ') == 0) write (iprt, '(A)') trim(koment)  
      if (index(title,' NULL ') == 0) write (iprt, '(A)') trim(title) 
      return  
      end subroutine wrttxt 
