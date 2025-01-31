      subroutine gettxt 
!...Translated by Pacific-Sierra Research 77to90  4.4G  10:47:18  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use upcase_I 
      use mopend_I 
      use chanel_C, only: ir, iw, isetup, input_fn
      use molkst_C, only: keywrd, koment, title, refkey,  gui, numcal, line, &
        in_house_only, verson, moperr
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k, l, ipath
      character :: filen*100, oldkey*1000, path*240
      logical :: aux, exists, setup_present, zero_scf
      character (len = 100), external :: get_text
!-----------------------------------------------
      koment = '    NULL  ' 
      title  = '    NULL  ' 
      refkey = '    NULL  '
      aux = (index(keywrd, "AUX") /= 0) 
      read (ir, '(A1000)', end=100, err=100) refkey(1)
      keywrd = refkey(1)
      oldkey = keywrd 
      call upcase (keywrd, len_trim(keywrd))
      zero_scf = (index(keywrd, "0SCF") /= 0) 
      do i = len_trim(input_fn), 2, -1
        if (input_fn(i:i) == "\" .or. input_fn(i:i) == "/") exit 
      end do
      ipath = i 
      if  (ipath > 2) then
        path = input_fn(:ipath)
        filen = trim(path)//'SETUP' 
      else
        filen = 'SETUP' 
      endif
      if (in_house_only) then
        if( index(keywrd,'PM6') == 0 .and. verson(7:7) == "M") filen = '/Users/jstewart/SETUP.txt' 
      end if
      inquire (file=filen, exist = exists)
      i = len_trim(keywrd) 
      setup_present = (index(keywrd,'SETUP') /= 0)
      if (setup_present .or. (in_house_only .and. exists)) then 
        i = index(keywrd,'SETUP=') 
        if (i /= 0) then 
          filen = get_text(oldkey, i + 6, 1) 
        else 
          if (in_house_only) then
            if (verson(7:7) == "M") then
              filen = 'SETUP' 
            else
              if  (ipath > 2) then
                filen = trim(path)//'SETUP' 
              else
                filen = 'SETUP' 
              endif
            end if
          else
            if  (ipath > 2) then
              filen = trim(path)//'SETUP' 
            else
              filen = 'SETUP' 
            endif
          end if
        endif 
        !call add_path(filen)
        inquire (file=filen, exist = exists)
        if (.not. exists) then
          inquire (file=trim(filen)//".txt", exist = exists)
          if (exists) filen = trim(filen)//".txt"
        end if
        if (.not. exists) then
          if (setup_present .and. .not. zero_scf) then
            write (iw, '(A)') ' SETUP FILE MISSING' 
            write(iw,'(a)') " (Setup file name: '"//trim(filen)//"')"
            numcal = 2
            if (.not. gui )write(0,'(//30x,a)')' SETUP FILE MISSING, EMPTY OR CORRUPT' 
            if (.not. gui )write(0,'(2x,a,//)')" (An attempt was made to open the Setup file named: '"//trim(filen)//"')"
            call mopend ('SETUP FILE MISSING') 
            return
          end if
        else 
          open(unit=isetup, file=filen, status='UNKNOWN', form='FORMATTED',position='asis', iostat=i) 
          if (i /= 0) then
            if (.not. zero_scf) then
              call mopend ('COULD NOT OPEN SETUP FILE: '//trim(filen)) 
              if (zero_scf) moperr = .false.
              return 
            end if
          end if
          rewind isetup           
          read (isetup, '(A)', end=50, err=50) refkey(2)
          close (isetup)
          call upcase (refkey(2), len(refkey(2)) ) 
!
!  Check for " -" signs in setup file
!
          if (refkey(2)(1:1) /= " ") refkey(2) = " "//refkey(2)(:len_trim(refkey(2)))
          do
            i = index(refkey(2), " -")
            if (i == 0) exit
            j = index(refkey(2)(i + 2:), " ") + i + 1
            do 
              k = index(" "//keywrd, " "//refkey(2)(i + 2: j - 1))
              if (k == 0) exit
              l = index(keywrd(k + 1:)," ") + k + 1
              keywrd(k:) = keywrd(l:)
            end do
            refkey(2)(i:) = refkey(2)(j:)
          end do
!
!   Check for keywords in SETUP that are present in the keyword line, and delete them
!   ( Keywords on keyword line take precedence.)
!
          i = 1
          do
            i = i + 1
            if (i > len_trim(refkey(2))) exit
            if (refkey(2)(i - 1:i - 1) == " " .and. refkey(2)(i:i) /= " ") then
!
!  Found a keyword. Now look for the end of the keyword
!
              do j = i + 1, len_trim(refkey(2))
                if (refkey(2)(j:j) == " ") exit
              end do
              line = refkey(2)(i:j) 
              k = min(6, j - i) 
              if (index(keywrd, " "//line(:k)) > 0) then
                refkey(2) = refkey(2)(:i - 1)//refkey(2)(j + 1:)
                i = i - 1
              end if
            end if
          end do
          i = len_trim(keywrd)
          keywrd(i + 1:) = " "//refkey(2)(:999 - i)
!
!  Check for " -" signs in keyword line
!
          do
            i = index(keywrd, " -")
            if (i == 0) exit
            j = index(keywrd(i + 2:), " ") + i + 1
            do 
              k = index(" "//keywrd, " "//keywrd(i + 2: j - 1))
              if (k == 0) exit
              l = index(keywrd(k + 1:)," ") + k + 1
              keywrd(k:) = keywrd(l:)
            end do
            i = index(keywrd, " -")
            j = index(keywrd(i + 2:), " ") + i + 1
            keywrd(i:) = keywrd(j:)
          end do
!
!  Remove excess EXTERNALS
!
          if (in_house_only) then
            do
              i = index(keywrd, " EXTER")
              if (i /= 0) then
                j = index(keywrd(i+3:), " EXTER")
                if (j /= 0) then
                  j = index(keywrd(i + 1:), " ") + i 
                  keywrd(i:) = " "//keywrd(j:)
                else
                  exit
                end if
              else
                exit
              end if
            end do
          end if
          refkey(1) = trim(keywrd)
          refkey(2) = refkey(3)
          if (keywrd(i + 1:) == ' ') go to 50 
         
        end if
        read (ir, '(A)', end=100, err=100) koment, title 
      else if (index(keywrd(1:i),' +') /= 0) then 
!
!  READ SECOND KEYWORD LINE
!
        i = index(keywrd(1:i),' +')
        keywrd(i:i + 1) = " "
        i = len_trim(keywrd)
        read (ir, '(A)', end=100, err=100) refkey(2)
        keywrd(i + 2:) = refkey(2)(:999 - i)
        oldkey = keywrd 
        call upcase (keywrd(i + 1:), len(keywrd) - i) 
        if (index(keywrd,'SETUP') /= 0) then 
          i = index(keywrd,'SETUP=') 
          if (i /= 0) then 
            j = index(keywrd(i:),' ') 
            filen = oldkey(i+6:i+j-1) 
            keywrd(i: i + j) = " "
          else 
            filen = 'SETUP' 
            i = index(keywrd,'SETUP') 
            keywrd(i: i + 5) = " "
          endif 
          keywrd(i:i+6) = " "
          !call add_path(filen)
          open(unit=isetup, file=filen, status='UNKNOWN', form='FORMATTED', &
            position='asis') 
          rewind isetup 
          read (isetup, '(A)', end=30, err=30) refkey(2)
          close(isetup)
          i = len_trim(keywrd) + 1
          keywrd(i:) = refkey(2)(:1001 - i)
          call upcase (keywrd, len_trim(keywrd)) 
   30     continue 
        else if (index(keywrd(i + 1:),' +') /= 0) then 
!
!  READ THIRD KEYWORD LINE
!
          read (ir, '(A)', end=100, err=100) refkey(3)
          i = index(refkey(3), " + ")
          if (i /= 0) then
            write(iw,"(a)")" A maximum of three lines of keywords are allowed."
            write(iw,"(a)")" On the third line of keywords is a '+' sign, implying more lines of keywords."
            write(iw,"(a)")" Remove the '+' sign from the third line of keywords, and re-run."
            call web_message(iw,"plus.html")
            call mopend("A maximum of three lines of keywords are allowed.")
            return
          end if
          i = index(keywrd(1:len_trim(keywrd)),' +')
          keywrd(i:i + 1) = " "
          i = len_trim(keywrd)
          keywrd(i + 2:) = refkey(3)(:1001 - i)
          call upcase (keywrd, len_trim(keywrd)) 
        endif 
!
!  READ TITLE LINE
!
        read (ir, '(A)', end=100, err=100) koment, title 
      else if (index(keywrd,'&') /= 0) then
        i = index(keywrd,'&')
        keywrd(i:i) = ' '
        i = len_trim(keywrd)
        read (ir, '(A)', end=100, err=100) refkey(2)
        keywrd(i + 1:) = " "//refkey(2)(:1001 - i)
        oldkey = keywrd
        call upcase (keywrd, len_trim(keywrd)) 
        if (index(keywrd,'SETUP') /= 0) then 
          i = index(keywrd,'SETUP=') 
          if (i /= 0) then 
            j = index(keywrd(i:),' ') 
            filen = oldkey(i + 6:i + j) 
            keywrd(i: i + j) = " "
          else 
            filen = 'SETUP' 
            i = index(keywrd,'SETUP') 
            keywrd(i:i + 6) = " "
          endif 
          !call add_path(filen)
          open(unit=isetup, file=filen, status='UNKNOWN', form='FORMATTED', &
            position='asis') 
          rewind isetup 
          read (isetup, '(A)', end=40, err=40) keywrd(len_trim(keywrd) + 2:) 
          close (isetup)
          call upcase (keywrd, len_trim(keywrd)) 
          read (ir, '(A)', end=100, err=100) title 
   40     continue 
        else if (index(keywrd(i:),'&') /= 0) then 
          read (ir, '(A)', end=100, err=100) keywrd(len_trim(keywrd):) 
          call upcase (keywrd, len_trim(keywrd)) 
        else 
          read (ir, '(A)', end=100, err=100) title 
        endif 
      else 
        read (ir, '(A)', end=100, err=100) koment, title 
      endif 
      go to 60 
50    continue 
      if (zero_scf) go to 60
      numcal = 2
      write (iw, '(A)') ' SETUP FILE MISSING, EMPTY OR CORRUPT' 
      write(iw,'(a)') " (Setup file name: '"//trim(filen)//"')"
      call mopend ('SETUP FILE MISSING, EMPTY OR CORRUPT') 
      return  
   60 continue 
      call upcase (keywrd, len_trim(keywrd)) 
      if (gui) then
        i = index(keywrd,"PM3")  ! Convert PM3 to PM6 for CAChe only
        if (i /= 0) then
          keywrd(i:i+2) = "PM6"
           write (iw, '(A)') ' Keyword PM3 was supplied. PM3 is not supported, so keyword converted to PM6' 
        end if
        i = index(keywrd,"PM5")  ! Convert PM5 to PM7 for CAChe only
        if (i /= 0) then
          keywrd(i:i+2) = "PM7"
           write (iw, '(A)') ' Keyword PM5 was supplied. PM5 is not supported, so keyword converted to PM7' 
        end if
      end if
      return  
  100 continue 
      if (numcal > 1) then
        if (index(keywrd,"OLDGEO") /= 0) then ! User forgot to add extra lines for title and comment
          return
        end if
        if (aux) keywrd = " AUX"
        line = "JOB ENDED NORMALLY"
      else if (index(keywrd, "GEO_DAT") == 0) then
        line = ' ERROR IN READ OF FIRST THREE LINES' 
      end if     
      if (line /= " ") call mopend (trim(line)) 
      return  
  end subroutine gettxt
  character (len = 100) function get_text(line, i_start, zero)
    implicit none
    character :: line*(*)
    integer :: i_start, zero
!
!  Return text between character i_start and the next space.
!  If character i_start is '"' or ''', return text between character i_start + 1 and the closing character.
!
    integer, parameter :: num_lim = 2
    character :: limit(num_lim)*1, ch*1
    integer :: i, j
    data limit /"""","'"/
    do i = 1, num_lim
      if (line(i_start:i_start) == limit(i)) exit
    end do
    j = i_start
    if (i > num_lim) then
      ch = " "
      i = i_start
    else
      j = j + 1
      ch = limit(i)
      i = i_start + 1
    end if
    do
      if (line(i + 1:i + 1) == ch) exit
      i = i + 1
    end do
    get_text = line(j:i)
    if (zero == 0) line(i_start:i + 1) = " "
    return  
  end function get_text
  
