      subroutine datin(iw)  
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      USE parameters_C, only : partyp
      use Common_arrays_C, only : ijpars, parsij
      use molkst_C, only : keywrd, lpars, verson
      use chanel_C, only : iext
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  10:47:05  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use upcase_I 
      use reada_I 
      use update_I
      implicit none
      integer, intent(in) :: iw
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      integer :: i, j, l, iparam, k, ielmnt, jelmnt, nref, &
      loop, mpar

      real(double) :: param
      logical :: exists
      character (len=60), dimension (20) :: file
      character :: text*80, line*241
      character, dimension(107) :: elemnt*2 
      save elemnt 
!----------------------------------------------- 
      data (elemnt(i),i=1,107)/ 'H ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ', 'O '&
        , 'F ', 'NE', 'NA', 'MG', 'AL', 'SI', 'P ', 'S ', 'CL', 'AR', 'K ', &
        'CA', 'SC', 'TI', 'V ', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'GA'&
        , 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y ', 'ZR', 'NB', 'MO', &
        'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE', 'I ', 'XE'&
        , 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', &
        'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF', 'TA', 'W ', 'RE', 'OS', 'IR'&
        , 'PT', 'AU', 'HG', 'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', &
        'AC', 'TH', 'PA', 'U ', 'NP', 'PU', 'AM', 'CM', 'BK', 'MI', 'XX', 'FM'&
        , 'MD', 'CB', '++', '+', '--', '-', 'TV'/  
    if (.not. allocated(ijpars))  allocate(ijpars(5,5000), parsij(5000))
    i = Index(keywrd, "EXTERNAL=") + Index(keywrd, "PARAMS=")
    i = i + Index (keywrd(i:), "=") ! i = start of parameter file name list
    j = Index (keywrd(i:), " ") + i ! j = end of parameter file name.
!
! In between "i" and "j" are the names of the parameter file name, separated by ";"
!
    nref = 0
    do l = 1,10
      k = Index(keywrd(i:j),";")
      if (k /= 0) then
        nref = nref + 1
        file(nref) = keywrd(i:i+k-2)
        i = i + k 
      end if
      if (k == 0) exit
    end do
    nref = nref + 1
    file(nref) = keywrd(i:j-1)
  !
  !   Read in parameters from a previous run - these will overwrite
  !   the default values of the parameters.
  !
    lpars = 0
    do 10 loop = 1, nref
      mpar = 0
      !call add_path(file(loop))
      inquire (file=trim(file(loop)), exist = exists)
      if (.not. exists) then
        if (index(keywrd,' 0SCF') == 0) call mopend("EXTERNAL file: '"//trim(file(loop))//"' does not exist!")
        line = trim(file(loop))
        do i = len_trim(line), 1, -1
          if (line(i:i) == "\" .or. line(i:i) == "/") then
            inquire (file=line(:i), exist = exists)
            if (exists) then
              write(iw,"(10x,a)")" (but folder: '"//line(:i)//"' does exist.)"
            else
              if (verson(7:7) == "W") write(iw,"(10x,a)")" (Note: the folder: '"//line(:i)//"' also does not exist.)"
            end if
            exit
          end if
        end do
        return
      end if
      open (unit=iext, form="FORMATTED", status="OLD", file=trim(file(loop)), action="READ", iostat = i)
      if (i /= 0) then
        if (lpars > 0) exit
        if (loop == 1) then
          write(line,'(a)')" EXTERNAL file """//trim(file(loop))//""" could not be opened"
          write(iw,'(/,a)')trim(line)
          if (index(keywrd,' 0SCF') + index(keywrd, " RESEQ") == 0 ) then 
            call mopend(trim(line))
            inquire (file=trim(file(loop)), exist = exists)
            if (exists) then
              write(line,'(a)')" (The EXTERNAL file exists, but could not be read)"
              write(iw,'(a)')trim(line)
              call mopend(trim(line))
            else
              write(line,'(a)')" (The EXTERNAL file does not exist)"
              write(iw,'(a)')trim(line)
              call mopend(trim(line))
            end if
          end if
        end if
      end if
      rewind (iext,err = 10)
      do
        read (iext, "(A60)", err=11, end=11) text
        call upcase (text, 80)
        if (Index (text, "END") /= 0 .or. text == " ")  exit
        if (text(1:1) == "*") cycle
        if (index(text, "PAR") /= 0) then
          i = index(text, "PAR") + 3
          j = ichar(text(i:i)) - ichar("0")
          i = i + 1
          line = text (1:80)
          if (text(i:i) /= " ") j = j*10 + ichar(text(i:i)) - ichar("0")     
          if (j < 13) then       
            k = mod(j - 1,3) + 1
            j = (j - 1)/3 + 1
            write(text,"('FN',2i1,' XX ',a)") k, j, trim(line(i + 1:))
          else
            j = j - 12
            k = mod(j - 1,3) + 1
            j = (j - 1)/3 + 1
             write(text,"('FN',2i1,' MI ',a)") k, j, trim(line(i + 1:))
           end if
        end if
        line = text (1:80)
        do
          if (line(1:1) /= " ") exit
          line = line(2:)
        end do
!
! Clean up line - delete anything after the third set of spaces
!
        do
          j = 0
          do i = 1, len_trim(line)
            if (line(i:i + 1) == "  ") then
              j = 1
              line(i:) = line(i + 1:)
            end if
          end do
          if (j == 0) exit
        end do
        i = Index (line, " ")
        i = index(line(i + 1:)," ") + i + 1
        i = index(line(i + 1:)," ") + i 
        text = line(1:i)
!
!  Force in spaces needed for parsing
!
        i = index(text," ")
        text(i:) = "  "//text(i:)
        i = index(text(i + 3:)," ") + i + 3
        text(i:) = "  "//text(i:)
        line = text
        do j = 1, 37
          if (Index(text, partyp(j)) /= 0) go to 1000
        end do
        write(iw,"(3a)")" EXTERNAL parameter type: '",trim(text),"' unrecognized"
        close(iext)
        goto 99 
1000    iparam = j
        jelmnt=0
        if (iparam == 34 .or. iparam == 35) then
!
!  This is a di-atomic parameter - read in the other element number
!
          i = Index(text, partyp(j))+5
          do j = 1, 99
            if (Index (" "//text(i:i+2), " "//elemnt(j)) /= 0) exit
          end do
          jelmnt = j
          end if
        i = Index (line, " ")
        text = line(i:)
        line = text
        do i = 1, 88
          if (line(1:1) /= " ") exit
          text = line(2:60)
          line = text
        end do
        text = " " // line (1:79)
        do j = 1, 100
          if (Index (text, " "//elemnt(j)) /= 0) go to 1100
        end do
        write(iw,"(3a)")" EXTERNAL element type: '",trim(text),"' unrecognized"
        close(iext)
        goto 99 
1100    param = reada (text, Index (text, " "//elemnt(j)))
        if (j > jelmnt) then
          ielmnt = j+200*jelmnt
        else
          ielmnt = jelmnt+200*j
        end if       
        do i = 1, lpars
          if (ijpars(1, i) == ielmnt .and. ijpars(2, i) == iparam) go to 1200
        end do
        lpars = lpars + 1
        i = lpars
1200    ijpars(1, i) = ielmnt
        ijpars(2, i) = iparam
        if (Abs(param) < 1.d-7) param = 1.d-7  ! Don't allow new parameters to be exactly zero
        parsij(i) = param
        text = " "
        mpar = 1
      end do
11  continue
    if(mpar == 1) then
      write(iw,'(3a)')" Parameters read in from file: """, file(loop)(:len_trim(file(loop))),""""
    else
      write(iw,'(3a)')" No parameters read in from """,file(loop)(:len_trim(file(loop))),""""
      if (index(keywrd,' 0SCF') + index(keywrd, " RESEQ") == 0 ) then 
        call mopend("No parameters read in from '"//file(loop)(:len_trim(file(loop)))//"'")
      end if
    end if
10  continue
    close (iext, status="KEEP")
    call write_params(iw)
    do i = 1, lpars
      call update(ijpars(2, i), ijpars(1, i), parsij(i), 0.d0)
    end do
    close(iext)
!99  deallocate(ijpars, parsij)
 99   return
    end subroutine datin
!
!
!
    subroutine write_params(iw)
      USE parameters_C, only : partyp
      use Common_arrays_C, only : ijpars, parsij
      use molkst_C, only : lpars, line
      implicit none
      integer, intent(in) :: iw
      logical :: lold = .true.
!
!  Local
!
      integer :: i, j, k, l, il, iu, ii, jj, iparam, ielmnt, jelmnt, j1
      character :: elemnt(107)*2, elemnt2*2
      save elemnt, lold
!----------------------------------------------- 
      data (elemnt(i),i=1,107)/ 'H ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ', 'O '&
        , 'F ', 'NE', 'NA', 'MG', 'AL', 'SI', 'P ', 'S ', 'CL', 'AR', 'K ', &
        'CA', 'SC', 'TI', 'V ', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'GA'&
        , 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y ', 'ZR', 'NB', 'MO', &
        'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE', 'I ', 'XE'&
        , 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', &
        'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF', 'TA', 'W ', 'RE', 'OS', 'IR'&
        , 'PT', 'AU', 'HG', 'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', &
        'AC', 'TH', 'PA', 'U ', 'NP', 'PU', 'AM', 'CM', 'BK', 'MI', 'XX', 'FM'&
        , 'MD', 'CB', '++', '+', '--', '-', 'TV'/  
      do j = 1, 107
        do k = 1, 37
          if (k == 35) cycle
          if (k == 34) then
           il = 1
           iu = 98
          else
           il = 0
           iu = 0
          end if
          do ii = il, iu
            do jj = 0, il
              do i = 1, lpars
                iparam = ijpars(2, i)
                ielmnt = ijpars(1, i)
                jelmnt = mod(ielmnt,200)
                l = ielmnt/200
                if (iparam == k + jj .and. jelmnt == j .and. l == ii) then
                  if (lold) then
                    write (iw, "(//,8X,A)") " Parameters read in"
                    write (iw, "(/5X, ' Parameter Type  Element    Parameter')")
                    lold = .false.
                  end if              
                  if(l /= 0) then
                    elemnt2 = elemnt(l)
                    if(elemnt2(2:2) == " ")elemnt2 = elemnt2(1:1)
                  else
                    elemnt2 = " "
                  endif
                  write (line, "(12X,A7,7X,A2,F16.6)") partyp (iparam)//elemnt2, elemnt(jelmnt), parsij(i)
                  if (ielmnt /= 99 .or. (iparam < 22 .or. iparam > 33)) then
                    write (iw, "(a)") trim(line)
                  else
                    j1 = iparam - 21               
                    if (j1 < 10) then
                      write(iw,"(12x,'PAR',i1,13x,a)") j1, trim(line(30:))
                    else
                      write(iw,"(12x,'PAR',i2,12x,a)") j1, trim(line(30:))
                    end if
                  end if
                end if
              end do
            end do
          end do
        end do
      end do
    end subroutine write_params
 
