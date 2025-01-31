subroutine wrtkey 
  use molkst_C, only : moperr, allkey
  implicit none
!**********************************************************************
!
!  WRTKEY CHECKS ALL KEY-WORDS AND PRINTS THOSE IT RECOGNIZES.  IF IT
!  FINDS A WORD IT DOES NOT RECOGNIZE THE PROGRAM WILL BE STOPPED.
!
!**********************************************************************
!
!   Write out the control keywords
!
  call wrtcon (allkey)
  if (moperr) return
!
!   Write out the keywords controlling the working 
!
  call wrtwor (allkey)  
!
!   Write out the keywords controlling the output
!
  call wrtout (allkey)
!
!   Check the keywords to see if two or more conflict
!
  call wrtchk (allkey)
end subroutine wrtkey



subroutine wrtchk (allkey)
  use molkst_C, only: keywrd, id, is_PARAM, uhf, method_mndo, method_am1, & 
  method_pm3, method_mndod, method_pm6, method_rm1, rhf, mozyme, line, &
  method_pm7, koment, title, refkey
  use chanel_C, only: input_fn
  use myword_I
  implicit none
  logical :: birad, exci, ci, trip
  integer :: i, j, k, l
  character :: ch
  character (len=1000), intent (inout) :: allkey
  birad = (Index (keywrd, " BIRAD") /= 0)
  exci  = (Index (keywrd, " EXCI") /= 0)
  ci    = (Index (keywrd, " C.I.") /= 0)
  trip  = (Index (keywrd, " TRIP") /= 0)
  uhf   = (Index (keywrd, " UHF") /= 0)
  rhf   = (Index (keywrd, " RHF") + index(keywrd, " MECI") /= 0  &
    .or. ci .or. (.not. uhf .and. Index(keywrd, " OPEN") /= 0))
  if (.not. (mozyme .or. (index(keywrd," PDBOUT") + index(keywrd," RESID") + index(keywrd," ADD-H")/= 0)) &
    .and. index(keywrd, " 0SCF") == 0) then
!
! Check for keywords that require MOZYME to be present
!
    if (index(keywrd," PDBOUT") /= 0) then
      call mopend("Keyword PDBOUT only works when MOZYME or 0SCF is also present")
      return
    end if
    if (index(keywrd," RESID") /= 0) then
      call mopend("Keyword RESIDUES only works when MOZYME or 0SCF is also present")
      return
    end if
    if (index(keywrd," CVB") /= 0) then
      call mopend("Keyword CVB only works with MOZYME")
      return
    end if
    if (index(keywrd," SETPI") /= 0) then
      call mopend("Keyword SETPI only works with MOZYME ")
      return
    end if
  end if
  if (mozyme) then
    if (index(keywrd, " UHF") /= 0) then
      call mopend("Keyword UHF cannot be used with MOZYME")
    end if
    if (index(keywrd, " ENPART") /= 0)  then
      call mopend("Keyword ENPART is not available with MOZYME")
    end if  
    if (index(keywrd, " LOCAL") /= 0)  then
      call mopend("Keyword LOCAL is not available with MOZYME")
    end if  
    if (index(keywrd," RE-LOC") /= 0) then
      i = index(keywrd," RE-LOC") 
      j = index(keywrd(i + 7:), " ") + i + 7
      if (index(keywrd(i:j), "=")  + index(keywrd," DENOUT") + &
        index(keywrd," VEC") + index(keywrd," ALLVE") == 0) then
        call mopend("Keyword RE-LOC only has meaning if DENOUT or VECTORS or ALLVECTORS is present")
      end if                                               
    end if  
    if (id /= 0) then
      if (index(keywrd, " CUTOF") /= 0 ) then
        call mopend("CUTOFx=n.nn -type keywords do not work with MOZYME for infinite systems")
        return
      end if
    end if
    if (index(keywrd," 1SCF") /= 0 .and. index(keywrd," RAPID") /= 0) then
      call mopend("RAPID cannot be used with 1SCF")
      return     
    end if
  end if
!
!   MULLIK does not work with UHF
!
  if (Index (keywrd, " MULLIK") /= 0 ) then
    if (uhf) then
      call mopend ("MULLIKEN POPULATION NOT AVAILABLE WITH UHF")
      return
    end if
  end if
!
!   Check UHF with nonallowed keywords
!
  if (uhf) then
    if (rhf) then
      call mopend("UHF and RHF cannot both be used")
    end if
    if (birad .or. exci .or. ci) then
      write (*, "(//10X, ' UHF USED WITH EITHER BIRAD, EXCITED OR C.I. ')")
      write (*, '(/ / 10 x, " IMPOSSIBLE OPTION REQUESTED,")')
      go to 1020
!
!  POLAR and UHF do not work
!
    else if (Index (keywrd, " POLAR") /= 0) then
      call mopend("POLAR does not work with UHF")
      go to 1020
    end if
  else if (exci .and. trip) then
    write (*, "(//10X,' EXCITED USED WITH TRIPLET')")
    write (*, 10000)
    go to 1020
  end if
!
!   PMEP only works with AM1
!
  if (Index (keywrd, " PMEP") /= 0 ) then
    if (.not. method_am1) then
      write (*, "(A)") " PMEP only works with AM1"
      call mopend ("PMEP only works with AM1")
      return
    end if
  end if
!
!  XYZ and INT cannot simultaneously be present
!
  if (Index (keywrd, " INT ") /= 0 .and. Index (keywrd, " XYZ") /= 0) then
    call mopend ("INT cannot be used with XYZ")
  end if
!
!  COSMO does not work with polymers or other infinite systems
!
  if (id > 0 .and. Index(keywrd, "EPS=") /= 0) then
    if (id == 1) write(*,*) " COSMO cannot be used with polymers"
    if (id == 2) write(*,*) " COSMO cannot be used with layer systems"
    if (id == 3) write(*,*) " COSMO cannot be used with solids"
    call mopend ("COSMO cannot be used with systems with Tv")
  end if
!
!   Check T-PRIO MUST have DRC
!
  if (Index (keywrd, " T-PRIO") /= 0 .and. Index (keywrd, " DRC") == 0) then
    write (*, "(//10X,'T-PRIO AND NO DRC')")
    write (*, 10000)
    go to 1020
  end if
!
!   Check that only one method is used
!
    i = 0
    if (method_am1)     i = i + 1
    if (method_pm3)     i = i + 1
    if (method_pm6)     i = i + 1
    if (method_pm7)     i = i + 1
    if (method_mndo)    i = i + 1
    if (method_mndod)   i = i + 1
    if (method_rm1)     i = i + 1
    if (i > 1) then
      write (*, '(//10 x, " ONLY ONE OF MNDO, MNDOD, AM1, PM3, RM1, AND PM6 ALLOWED")')            
10000 format(/ / 10 x, "IMPOSSIBLE OPTION REQUESTED")
      call mopend("IMPOSSIBLE OPTION REQUESTED")
      return
    else
       !
       !   Check that only one geometry option is used
       !
      i = 0
      line = trim(keywrd)
      j = 0
      l = len_trim(line)
      do 
        j = j + 1
        if (line(j:j) == '"') then
          k = j
          do 
            j = j + 1
            if (line(j:j) == '"') exit
          end do
          line(k:j) = " "
        end if
        if (j >= l) exit
      end do       
      if (Index (line, " BFGS") /= 0)    i = i + 1
      if (Index (line, " LBFGS") /= 0)   i = i + 1
      if (Index (line, " EF") /= 0)      i = i + 1
      if (Index (line, " TS") /= 0)      i = i + 1
      if (Index (line, " SIGMA") /= 0)   i = i + 1
      if (Index (line, " NLLSQ") /= 0)   i = i + 1
      if (Index (line, " FORCE") + Index (line, " IRC") + Index (line, " DRC") /= 0)   i = i + 1
      if (i > 1) then
        call mopend ("MORE THAN ONE GEOMETRY OPTION HAS BEEN" // &
       & " SPECIFIED. CONFLICT MUST BE RESOLVED BEFORE JOB WILL RUN.")
        return
      end if
      if (Index (keywrd, " HESSIAN") /= 0) then
        if (Index (keywrd, " EF") == 0) then
          line = " Keyword EF must be present if HESSIAN is used"
          write(*,"(a)")trim(line)
          call mopend(trim(line))
          return
        end if 
      end if
    end if
    if (Index (keywrd, " PKA") /= 0 .and. .not. method_pm6) then
      write(*,"(a)")" *"," *  The pKa option only works with PM6."
      write(*,"(a)")" *  Either remove keyword pKa or change method to PM6."," *"
      call mopend("The pKa option only works with PM6. Calculation abandoned.")
        return
    end if
       !******************************************************************
       !
       !  Check to see that all keywords have been recognized
       !
       !*******************************************************************
       !
       !    DUMMY IF STATEMENT TO REMOVE AMPERSAND, PLUS SIGNS AND OBSOLETE KEYWORDS, IF PRESENT
       !
      if (myword(allkey(160:), " SETUP")) i = 1
      if (myword(allkey, "&"))            i = 2
      if (myword(allkey, " +"))           i = 3
      if (myword(allkey, " DIIS"))    write (*,'(" *  DIIS       - THIS IS AN OBSOLETE KEYWORD, IT WILL BE IGNORED")')
      if (myword(allkey, " NODIIS"))  write (*,'(" *  NODIIS     - THIS IS AN OBSOLETE KEYWORD, IT WILL BE IGNORED")')
      if (myword(allkey, " ROT"))     write (*,'(" *  ROT        - THIS IS AN OBSOLETE KEYWORD, IT WILL BE IGNORED")')
      if (myword(allkey, " HEADER") .or. myword(allkey, " USER"))  then
        write (*,'(" *  HEADER     - DATA SET IS IN PROTEIN DATA BANK FORMAT")')
        i = index(keywrd, " HEADER")
        keywrd = " ADD-H PDBOUT "//keywrd(:i)
        refkey(1) = trim(keywrd)
        refkey(2) = " NULL "
        koment = " From PDB file: '"//input_fn(1:len_trim(input_fn)   - 5)//"'"
        title = "   "
        return
      end if

      if (allkey == " ") return
!
!  Compress unrecognized key-words
!
   if (allkey /= ' ' .and. .not. is_PARAM) then 
        j = 0 
        do i = 1, 999 
          if (allkey(i:i+1) == '  ') cycle  
          j = j + 1 
          ch = allkey(i:i) 
          allkey(j:j) = ch 
        end do 
        j = max(1,j) 
        l = index(keywrd,' DEBUG') 
        if (l /= 0) then 
          write (*, '('' *  DEBUG KEYWORDS USED:  '',A)') allkey(:j) 
        else
          !write (line, "('UNRECOGNIZED KEY-WORDS: (',A,')')") allkey (:j)
          !call mopend (trim(line))
          !call mopend ("IF THESE ARE DEBUG KEYWORDS, ADD THE KEYWORD ""DEBUG"".")
        end if
      return
    end if
1020 if (is_PARAM) return
  write (*, "(//10X,' CALCULATION ABANDONED, SORRY!')")
  call mopend ("CALCULATION ABANDONED, SORRY!")
end subroutine wrtchk


subroutine wrtcon (allkey)
  use molkst_C, only: keywrd, numat, pressure, id, mozyme, mers, natoms, id, &
    line, old_chrge
  use cosmo_C, only : fepsi, nspa
  use meci_C, only : nmos
  use conref_C, only : fpcref
  use common_arrays_C, only :  lopt
  use myword_I
  use reada_I
  implicit none
  character (len=1000), intent (out) :: allkey
  double precision :: epsi
  integer :: i, ielec, ilevel, j, k, l
  logical :: l_add_H = .false., l_temp
  character (len = 100), external :: get_text
  save :: l_add_H
  i = index(keywrd, " LOCATE_TS") 
  if (i /= 0) keywrd(i:i + 9) = " LOCATE-TS"
  l_temp = (index(keywrd," MOZ") + index(keywrd," LEWIS") + index(keywrd," LOCATE-TS") + &
    index(keywrd, " RESEQ") + index(keywrd, " CHARGES") + index(keywrd," REFINE-TS") + &
    index(keywrd, " RAPID") + index(keywrd, " SITE=") + index(keywrd, " SITE(")  /= 0)   
  mozyme = (l_temp .or. + index(keywrd, " PDBOUT") /= 0)   
  if (index(keywrd, "CHARGE=") == 0 .and. old_chrge /= 0 .and. index(keywrd, " OLDGEO") /= 0) then
!
!  Use the charge from the previous calculation
!
    i = index(keywrd,"           ")
    if (old_chrge > 99) then
      write(keywrd(i:i+11),'(" CHARGE=",i3)')old_chrge
    else if (old_chrge > 9) then
      write(keywrd(i:i+11),'(" CHARGE=",i2)')old_chrge
    else if (old_chrge > 0) then
      write(keywrd(i:i+11),'(" CHARGE=",i1)')old_chrge
    else if (old_chrge > -10) then
      write(keywrd(i:i+11),'(" CHARGE=",i2)')old_chrge
    else 
      write(keywrd(i:i+11),'(" CHARGE=",i3)')old_chrge
    end if    
  else
    old_chrge=0 
  end if
  allkey = trim(keywrd)
  if (myword(allkey, " CCDC "))  i = 0 ! Dummy assignment   - to clear CCDC
  if (myword(allkey, " MNDO "))  write (*, '(" *  MNDO       - The MNDO Hamiltonian to be used")')
  if (myword(allkey, " AM1 "))   write (*, '(" *  AM1        - The AM1 Hamiltonian to be used")')
  if (myword(allkey, " PM3 "))   write (*, '(" *  PM3        - The PM3 Hamiltonian to be used")')
  if (myword(allkey, " MNDOD"))  write (*, '(" *  MNDOD      - The MNDOD Hamiltonian to be used")')
  if (myword(allkey, " PM6 "))   write (*, '(" *  PM6        - The PM6 Hamiltonian to be used")')
  if (myword(allkey, " PM7-TS")) write (*, '(" *  PM7-TS     - Calculate barrier height using PM7-TS")')
  if (myword(allkey, " PM7"))    write (*, '(" *  PM7        - The PM7 Hamiltonian to be used")')
  if (myword(allkey, " SPARKL")) write (*, '(" *  SPARKLE    - Use SPARKLES when they exist.")')
  if (myword(allkey, " RM1 "))   write (*, '(" *  RM1        - The RM1 Hamiltonian to be used")')
  if (myword(allkey, " PM5 "))   write (*, '(" *  PM5        - This Hamiltonian is not supported")')
!
!  The lack of space before QMMM on the next line is deliberate   - at allows MOL_QMMM as an option
!
  if (myword(allkey, "QMMM "))   write (*, '(" *  QMMM       - Generate energies and gradients for use in MM codes")')
  if (myword(allkey, " BZ"))     write (*, '(" *  BZ         - Write file <name>.brz for use by program BZ")')
  if (myword(allkey, " BIRAD"))  write (*, '(" *  BIRADICAL- SYSTEM HAS TWO UNPAIRED ELECTRONS")')
  if (myword(allkey, " EXCI"))   write (*, '(" *  EXCITED    - FIRST EXCITED STATE IS TO BE OPTIMIZED")')
  if (myword(allkey, " VELO"))   write (*, '(" *  VELOCITY   - INPUT STARTING VELOCITIES FOR DRC")')
  if (myword(allkey, " GEO-OK")) write (*, '(" *  GEO-OK     - OVERRIDE INTERATOMIC DISTANCE AND OTHER SAFETY CHECKS")')
  if (myword(allkey, " CHECK"))  write (*, '(" *  CHECK      - RUN EXTRA INTERATOMIC DISTANCE CHECKS")')
  if (myword(allkey, " PM6-D")) then 
    if (index(keywrd, "PM6-DH+") /= 0) then
      write (*, '(" *  PM6-DH+    - CORRECT DISPERSION AND HYDROGEN BOND TERMS USING PM6-DH+")')
    elseif (index(keywrd, "PM6-DH2X") /= 0) then
      write (*, '(" *  PM6-DH2X   - CORRECT DISPERSION, HALOGEN AND HYDROGEN BOND TERMS USING PM6-DH2X")')
    elseif (index(keywrd, "PM6-D2X") /= 0) then
      write (*, '(" *  PM6-D2X    - CORRECT DISPERSION AND HALOGEN BOND TERMS USING PM6-D2X")')
    elseif (index(keywrd, "PM6-D3 ") /= 0) then
      write (*, '(" *  PM6-D3     - CORRECT DISPERSION USING GRIMME''s D3 METHOD")')
    elseif (index(keywrd, "PM6-D3H4") /= 0) then
      write (*, '(" *  PM6-D3H4   - CORRECT DISPERSION AND HYDROGEN BOND TERMS USING THE D3H4 METHOD")')
    elseif (index(keywrd, "PM6-D3(H4)") /= 0) then
      write (*, '(" *  PM6-D3(H4) - CORRECT DISPERSION USING THE D3H4 METHOD")')
    elseif (index(keywrd, "PM6-DH2") /= 0) then
      write (*, '(" *  PM6-DH2    - CORRECT DISPERSION AND HYDROGEN BOND TERMS USING PM6-DH2")')
    else 
      write (*, '(" *  PM6-D      - CORRECT DISPERSION TERMS USING PM6-D")')
    end if
  else if (myword(allkey, " PM6-H"))  then  
    write (*, '(" *  PM6-H      - CORRECT HYDROGEN BOND TERMS USING PM6-H")')
  end if
  if (myword(allkey, " JCTC")) i = 6
  if (myword(allkey, " AIGIN"))  write (*, '(" *  AIGIN      - GEOMETRY MUST BE IN GAUSSIAN FORMAT")')
  if (myword(allkey, " ESR"))    write (*, '(" *  ESR        - RHF SPIN DENSITY CALCULATION REQUESTED")')
  if (myword(allkey, " NOMM"))   write (*, '(" *  NOMM       - DO NOT MAKE MM CORRECTION TO CONH BARRIER")')
  if (myword(allkey, " MMOK"))   write (*, '(" *  MMOK       - APPLY MM CORRECTION TO CONH BARRIER")')
  if (myword(allkey, " CIS "))   write (*, '(" *  CIS        - C.I. USES 1 ELECTRON EXCITATIONS ONLY")')
  if (myword(allkey, " CISD "))  write (*, '(" *  CISD       - C.I. USES 1 AND 2 ELECTRON EXCITATIONS")')
  if (myword(allkey, " CISDT ")) write (*, '(" *  CISDT      - C.I. USES 1, 2 AND 3 ELECTRON EXCITATIONS")')
  if (myword(allkey, " SING"))   write (*, '(" *  SINGLET    - SPIN STATE DEFINED AS A SINGLET")')
  if (myword(allkey, " DOUB"))   write (*, '(" *  DOUBLET    - SPIN STATE DEFINED AS A DOUBLET")')
  if (myword(allkey, " TRIP"))   write (*, '(" *  TRIPLET    - SPIN STATE DEFINED AS A TRIPLET")')
  if (myword(allkey, " QUAR"))   write (*, '(" *  QUARTET    - SPIN STATE DEFINED AS A QUARTET")')
  if (myword(allkey, " QUIN"))   write (*, '(" *  QUINTET    - SPIN STATE DEFINED AS A QUINTET")')
  if (myword(allkey, " SEXT"))   write (*, '(" *  SEXTET     - SPIN STATE DEFINED AS A SEXTET")')
  if (myword(allkey, " SEPT"))   write (*, '(" *  SEPTET     - SPIN STATE DEFINED AS A SEPTET")')
  if (myword(allkey, " OCTE"))   write (*, '(" *  OCTET      - SPIN STATE DEFINED AS A OCTET")')
  if (myword(allkey, " NONE"))   write (*, '(" *  NONET      - SPIN STATE DEFINED AS A NONET")')
  if (myword(allkey, " COSCCH")) write (*, '(" *  COSCCH     - ADD IN COSMO CHARGE CORRECTIONS")')
  if (myword(allkey, " FIELD"))  write (*, '(" *  FIELD      - AN EXTERNAL ELECTRIC FIELD IS TO BE USED")')
  if (myword(allkey, " NOREOR")) write (*, '(" *  NOREOR     - DO NOT REORIENTATE SYSTEM WHEN WORKING OUT THE SYMMETRY")')
  if (myword(allkey, " INVERT")) write (*, '(" *  INVERT     - REVERSE ALL OPTIMIZATION FLAGS")')
  if (myword(allkey, " ESP "))   write (*, '(" *  ESP        - ELECTROSTATIC POTENTIAL CALCULATION")')
  if (myword(allkey, " NSURF"))  write (*, '(" *  NSURF      - NUMBER OF LAYERS")')
  if (myword(allkey, " NOGPU"))  write (*, '(" *  NOGPU      - DO NOT USE GPU ACCELERATION")')
  if (myword(allkey, " SCALE"))  write (*, '(" *  SCALE      - SCALING FACTOR FOR VAN DER WAALS DISTANCE")')
  if (myword(allkey, " SCINCR")) write (*, '(" *  SCINCR     - INCREMENT BETWEEN LAYERS")')
  if (myword(allkey, " SLOPE"))  write (*, '(" *  SLOPE      - SLOPE   - USED TO SCALE MNDO ESP CHARGES")')
  if (myword(allkey, " PDBOUT")) write (*, '(" *  PDBOUT     - PRINT GEOMETY IN PDB FORMAT")')
  if (myword(allkey, " XENO"))   write (*, '(" *  XENO       - FRAGMENTS ARE ATTACHED TO RESIDUES")')
  if (myword(allkey, " CHARGES"))write (*, '(" *  CHARGES    - IDENTIFY AND PRINT CHARGED ATOMS")')
  if ( .not. l_add_H) l_add_H = (index(keywrd, " ADD-H") /= 0)
  if (myword(allkey, " RESEQ"))  write (*, '(" *  RESEQ      - RESEQUENCE Z-MATIX INTO NORMAL PDB FORMAT")')
  if (myword(allkey, " NORES"))  write (*, '(" *  NORESEQ    - DO NOT RESEQUENCE Z-MATIX INTO NORMAL PDB FORMAT")')
  i = index(keywrd," START_RES") + 1
  if (i > 1) then
    j = index(keywrd(i:),")") 
    if (j == 0) then
      write (*, '(" ***START_RES must be followed by open and close parentheses, e.g., START_RES=(1 10)")')
      call mopend("Keyword START_RES not properly defined")
      return
    end if
    j = j + i
    allkey(i:j) = " "
    write (*, '(" *  START_RES- STARTING RESIDUE NUMBERS DEFINED")')
    write (*, '(" *  Keyword:   ",a)')keywrd(i:j)
    if(myword(allkey, "START_RES")) then
      call mopend("Only one keyword START_RES allowed")
      return
    end if
  end if   
  if (myword(allkey, " GEO_DAT")) then
    i = index(keywrd," GEO_DAT")
    j = index(keywrd(i + 10:),'"') + i + 9 
    write (*, '(" *  GEO_DAT    - DATA SET GEOMETRY IS IN FILE """,a,"""")')keywrd(i+10:j   - 1)
    allkey(i:j) = " "
    if (index(keywrd," SETPI") /= 0 .and. index(keywrd," SETPI=") == 0) then
      write(line,"(a)")" When GEO_DAT is used, SETPI must be in the form SETPI=<file>.txt "
      call mopend(trim(line))
      return
    end if
  end if
  if (myword(allkey, " GEO_REF")) then
    if (index(keywrd, " RESID") /= 0) then
      call mopend("RESIDUES must not be used with GEO_REF")
      return
    end if
    i = index(keywrd," GEO_REF")
    j = index(keywrd(i + 10:),'"') + i + 9 
    if (index(keywrd(i:j), "SELF") /= 0) then
       write (*, '(" *  GEO_REF=""SELF""    - USE MOPAC DATA SET (<file>.mop) AS REFERENCE GEOMETRY")')
    else   
      write (*, '(" *  GEO_REF    - REFERENCE GEOMETRY IS IN FILE """,a,"""")')keywrd(i+10:j   - 1)
    end if
    j = j + 1
    k = j
    if (allkey(j:j) /= " ") then                               
      j = index(keywrd(j:), " ") + j   - 2
      write (*, '(3a)')" *  The constant 'c' in perturbation = sum( c*(deltax)**2 is set to "//allkey(k:j)//" kcal.mol-1.A-2"
    end if
    allkey(i:j) = " "
  end if
  i = index(allkey," CHAIN") + 1
  if (i > 1) then
    j = index(keywrd(i:),")") 
    if (j == 0) then
      write (*, '(" ***CHAIN must be followed by open and close parentheses, e.g., CHAIN=(ABCD)")')
      call mopend("Keyword CHAIN not properly defined")
      return
    end if
    j = j + i
    allkey(i:j) = " "
                                 write (*, '(" *  CHAIN      - PDB CHAIN LETTERS EXPLICITLY DEFINED")')
                                 write (*, '(" *  Keyword:   ",a)')keywrd(i:j)
  end if     
  if (myword(allkey, " GEOCHK")) write (*, '(" *  GEOCHK     - PRINT WORKING IN SUBROUTINE GEOCHK")') 
  if (myword(allkey, " LEWIS"))  write (*, '(" *  LEWIS      - PRINT OUT LEWIS STRUCTURE, THEN STOP")')
  if (myword(allkey, " SETPI"))  write (*, '(" *  SETPI      - SOME OR ALL PI BONDS EXPLICITLY SET BY USER")')
  if (myword(allkey, " PDB "))   write (*, '(" *  PDB        - INPUT GEOMETRY IS IN PDB FORMAT")')
  if (myword(allkey, " PDB("))   write (*, '(" *  PDB(txt)   - SYMBOLS IN PDB FILE ARE DEFINED BY USER")')
  if (myword(allkey, " MOZ"))    write (*, '(" *  MOZYME     - USE LOCALIZED M.O.s IN SOLVING THE SCF EQUATIONS")')
  if (myword(allkey, " RAPID"))  write (*, '(" *  RAPID      - IN MOZYME SCF, USE ATOMS BEING OPTIMIZED ONLY")')
  if (myword(allkey, " PINOUT")) write (*, '(" *  PINOUT     - WRITE OUT THE LMOs WHEN READING OR WRITING A ''.DEN'' FILE")')
  if (myword(allkey, " REORTH")) write (*, '(" *  REORTH     - RE-ORTHOGONALIZE LMOs EACH 10 SCF CALCULATIONS")')
  if (myword(allkey, " RE-LOCAL")) write(*,'(" *  RE-LOCAL   - RE-LOCALIZE LMOs AT END OF CALCULATION")')
  if (myword(allkey, " SWAP"))   write(*,  '(" *  SWAP       - THIS KEYWORD IS NOW REDUNDANT. SEE KEYWORD NOSWAP")')
  if (myword(allkey, " NOSWAP")) write(*,  '(" *  NOSWAP     - DO NOT SWAP ATOMS EVEN IF IT WILL IMPROVE OVERLAP IN GEO_REF")')
  if (myword(allkey, " A0 "))    write (*, '(" *  A0         - INPUT GEOMETRY IS IN ATOMIC UNITS (A0)")')
  if (myword(allkey, " ANG"))    write (*, '(" *  ANGSTROMS- INPUT GEOMETRY IS IN ANGSTROMS")')
  if (myword(allkey, " ADD-H"))  write (*, '(" *  ADD-H      - ADD HYDROGEN ATOMS TO SATISFY VALENCE")')
  if (myword(allkey, " CONNOL")) write (*, '(" *  CONNOLLY   - USE CONNOLLY SURFACE")')
  if (myword(allkey, " ESPRST")) write (*, '(" *  ESPRST     - RESTART OF ELECTRIC POTENTIAL CALCULATION")')
  if (myword(allkey, " UHF"))    write (*, '(" *  UHF        - UNRESTRICTED HARTREE-FOCK CALCULATION")')
  if (myword(allkey, " RHF"))    write (*, '(" *  RHF        - RESTRICTED HARTREE-FOCK CALCULATION")')
  if (myword(allkey, " STATIC")) write (*, '(" *  STATIC     - CALCULATE STATIC FIELD POLARIZABILITIES")')
  i = index(allkey," SETUP=")
  if (i /= 0) then
    line = get_text(allkey, i + 7, 0)  ! delete file name plus delimiter, if any.
    allkey(i:i + 7) = " "
                                 write (*, '(" *  SETUP      - EXTRA KEYWORDS TO BE READ FROM FILE '''//trim(line)//'''")')
  else
    if (myword(allkey, " SETUP")) write (*, '(" *  SETUP      - EXTRA KEYWORDS TO BE READ FROM FILE SETUP")')
  end if  
  if (myword(allkey, " MAX"))    write (*, '(" *  MAX        - GRID SIZE 23*23 ")')
  if (myword(allkey, " COSWRT")) write (*, '(" *  COSWRT")')
  if (myword(allkey, " OLDCAV")) write (*, '(" *  OLDCAV")')
  if (myword(allkey, " SYM ") .or. myword(allkey, " SYMM")) &
    write (*, '(" *  SYMMETRY   - SYMMETRY CONDITIONS TO BE IMPOSED")')
  if (myword(allkey, " POLAR"))  write (*, '(" *  POLAR      - CALCULATE FIRST, SECOND AND THIRD-ORDER POLARIZABILITIES")')
  if (myword(allkey, " RESI")) then
    write (*, '(" *  RESIDUES   - DETERMINE THE SEQUENCE OF RESIDUES")')
    if (numat /= natoms - id) then
      call mopend("DUMMY ATOMS CANNOT BE PRESENT WHEN KEYWORD ""RESIDUE"" IS PRESENT")
      write(*,'(10x,a)')"(A simple way to remove dummy atoms is to add keyword ""XYZ"")"
    end if
    if (index(keywrd," RESIDUES0") /= 0) then
      if (index(keywrd," ADD-H") + index(keywrd," RESEQ") + index(keywrd, " SITE=") + index(keywrd, " SITE(") /= 0) then
        call mopend("Keyword RESIDUES0 cannot be used when a keyword that changes the geometry is present")
        if (index(keywrd," ADD-H") /= 0) write(*,'(10x,a)')"Keyword ADD-H changes the geometry."
        if (index(keywrd," RESEQ") /= 0) write(*,'(10x,a)')"Keyword RESEQ changes the geometry."
        if (index(keywrd," SITE=") /= 0) write(*,'(10x,a)')"Keyword SITE changes the geometry."
        if (index(keywrd," SITE(") /= 0) write(*,'(10x,a)')"Keyword SITE changes the geometry."
      end if
    else
      if (index(keywrd," 0SCF") + index(keywrd," MOZ") + index(keywrd," LEWIS") + index(keywrd," ADD-H") + &
      index(keywrd," CHARGES") + index(keywrd," RESEQ") + index(keywrd, " SITE=") + &
      index(keywrd, " SITE(") == 0) then
      call mopend("Keyword RESIDUES can only be used when one of 0SCF, MOZYME, LEWIS, CHARGES, or RESEQ is present")
      end if
    end if
  end if
  if (myword(allkey, " NORES"))  &
    write (*, '(" *  NORES      - THIS IS THE DEFAULT.  USE ""RESIDUES"" IF RESIDUES ARE TO BE CALCULATED")')
  i = index(keywrd," SITE(")
  if (i > 0) keywrd(i:) = " SITE="//keywrd(i + 5:)
  if (myword(allkey, " SITE=") .or. myword(allkey, " SITE(") )   then
    write (*, '(" *  SITE       - SET IONIZATION LEVELS OF IONIZABLE RESIDUES ")')
    if (id /= 0) then
       write(*,'(1x,a)')"*","*","***Keyword SITE=(An(x)) can only be used for systems without Tv","*","*"
       call mopend("Keyword SITE=(An(x)) can only be used for systems without Tv")
    end if
    i = index(keywrd, " SITE=")
    do j = i, len_trim(keywrd)
      if (keywrd(j:j + 1) == ") ") exit
    end do
    allkey(i:j) = " "
  end if
  if (myword(allkey, " NOOPT"))   then
    line = trim(keywrd)
    do
      i = index(line, " NOOPT")
      if (i == 0) exit
      j = index(line(i + 5:)," ") + i + 3
      if (j   - i == 8) then
        line(j:j) = char( ichar(line(j:j)) + ichar("a")   - ichar("A"))
        write (*, '(" *  NOOPT-",a2,"   - DO NOT OPTIMIZE COORDINATES OF ATOMS OF TYPE ",a2)')line(j - 1:j), line(j   - 1:j)
      else if (j   - i == 7) then
        write (*, '(" *  NOOPT-",a1,"    - DO NOT OPTIMIZE COORDINATES OF ATOMS OF TYPE ",a1)')line(j:j), line(j:j)
      else if (j   - i == 5) then
        write (*, '(" *  NOOPT      - DO NOT OPTIMIZE ANY COORDINATES")')
      else  
        write (*, '(" ***NOOPT      - IMPROPER USE OF NOOPT KEYWORD, KEYWORD USED = ''", a,"''")')line(i + 1:j)
        call mopend("IMPROPER USE OF OPT KEYWORD")
        return
      end if
      line(i:j) = " "
    end do 
  end if
  line = trim(allkey)
  if (myword(allkey, " OPT-") .or. myword(allkey, " OPT "))   then    
    do
      i = index(line, " OPT")
      if (i == 0) exit
      j = index(line(i + 4:)," ") + i + 2
      if (j   - i == 6) then
        line(j:j) = char( ichar(line(j:j)) + ichar("a")   - ichar("A"))
        write (*, '(" *  OPT-",a2,"     - OPTIMIZE COORDINATES OF ATOMS OF TYPE ",a2)')line(j   - 1:j), line(j   - 1:j)
      else if (j   - i == 5) then
        write (*, '(" *  OPT-",a1,"      - OPTIMIZE COORDINATES OF ATOMS OF TYPE ",a1)')line(j:j), line(j:j)
      else if (j   - i == 3) then
        write (*, '(" *  OPT        - OPTIMIZE COORDINATES OF ALL ATOMS")')
      else  
        write (*, '(" ***OPT        - IMPROPER USE OF OPT KEYWORD, KEYWORD USED = ''", a,"''")')line(i + 1:j)
        call mopend("IMPROPER USE OF OPT KEYWORD")
        return
      end if
      line(i:j) = " "
    end do 
  end if
!
!  Keywords involving equals sign
!
!**********************************************************************
!
    !  Tension or pressure in a solid-state calculation
    !
    pressure = 0.d0
    if (myword(allkey, " P=")) then
      pressure = reada (keywrd, Index (keywrd, " P="))
      if (id == 1) then
        write (*, '(" *  P=         - TENSION IN POLYMER=", g13.6, &
       & " NEWTONS PER MOLE")') pressure
        pressure = -pressure * 10.d0 ** (-13) / 4.184
      else if (id == 2) then
      else if (id == 3) then
        i = Index (keywrd, " P=")
        j = Index (keywrd(i + 3:), " ") + i + 3
        if (index(keywrd(i + 3: j),"GP") /= 0) then
         write (*, '(" *  P=         - PRESSURE ON SOLID=", f7.3, &
             & " GIGAPASCALS")') pressure
          pressure = pressure*1.d9
        else
          write (*, '(" *  P=         - PRESSURE ON SOLID=", g13.6, &
             & " NEWTONS PER SQUARE METER")') pressure
        endif
!
!  Multiply by N, Avogadro's Number, to convert from J/M**3 per molecule 
!  to J/M**3/mol.
!  Divide by 4184 to convert from J/M**3/mol to Kcal/M**3/mol
!  Divide by 10**30 to convert from Kcal/M**3/mol to Kcal/Angstrom**3/mol
!
        pressure = -(fpcref(1,10)*pressure) / (4184.d0*10.d0**30)
      else
        write (*, *) " Keyword 'P=n.nn' is not allowed here"
        call mopend("Keyword 'P=n.nn' is not allowed here")
        return
      end if
    end if
!
!                       ESP grid
!
    i = -100
    if (myword(allkey, " CUBE")) then
      write (*,'(" *  CUBE       - USE GAUSSIAN CUBE FILE FOR ELECTROSTATICS")')    
    else if (myword(allkey, " ESPGRID")) then
      i = Nint (reada (keywrd, Index (keywrd, " ESPGRID")))
      write (*,'(" *  ESPGRID    - GENERATE ELECTROSTATICS CUBE FILE  =", i5," POINTS ON SIDE")') i
    end if
    if (i /= -100 .and. i < 2 .or. i > 100) then
      if (i < 2) write (*,'(//," *   Value of number of points is too small",//)')
      if (i > 100) write (*,'(//," *   Value of number of points is too large",//)')
      call mopend("Number of points in electrostatics is outside limits")
      return
    end if 
     
!
!
!               Electronic quantities (mainly C.I.)
!
!                       CHARGE
!
  if (myword(allkey, " CHARGE=")) then
    i = Nint (reada (keywrd, Index (keywrd, " CHARGE=")))
    if (i == 0) then
       write (*,'(3(" *", /), " *", 15 x, "  CHARGE ON SYSTEM = 0", 3 (/, " *"))')
    else
        write (*,'(3(" *", /), " *", 15 x, "  CHARGE ON SYSTEM =", SP,i6, 3 (/, " *"))') i
    end if
  end if
!
!                       C.I.=(n1,n2)
!
  if (Index(allkey, " C.I.") /= 0) then
    j = Index (keywrd, " C.I.=(")
    if (j /= 0) then
      j = Index (keywrd(j:j+10), ",") + j   - 1
      nmos = Nint (reada (keywrd, Index (keywrd, "C.I.=(")+5))
      write (*,'(" * C.I.=(N,M)-", i3, " DOUBLY FILLED LEVELS USED IN A",/, &
     & " *             C.I. INVOLVING ", i2, " M.O.''S")')  &
     Int (reada (keywrd, j)), nmos
     if (.not. myword(allkey, " C.I.=(")) return ! Impossible option used to delete keyword
    else if (myword(allkey, " C.I.=")) then
      nmos = Int (reada (keywrd, Index (keywrd, "C.I.")+5))
      write (*,' (" *  C.I.=N   -", i2, " M.O.S TO BE USED IN C.I.")') nmos  
    else    
      write (*, "(' *',/,a)") " *      C.I. keyword must be of form 'C.I.=n' or 'C.I.=(n1,n2)' (See manual)"
      call mopend("C.I. keyword must be of form 'C.I.=n' or 'C.I.=(n1,n2)'")
      return
    end if
    if (nmos > 29) then
      write (*, "(' *',/,a)") " *      Maximum size of open space = 29 M.O.s"
      write (*, "(a, i3,/,' *')") " *      Size requested:", nmos
      call mopend("Maximum size of open space = 29 M.O.s")
      return
    else if (nmos < 1) then
      write (*, "(' *',/,a)") " *      Minimum size of open space = 1 M.O.s"
      write (*, "(a, i3,/,' *')") " *      Size requested:", nmos
      call mopend("Minimum size of open space = 1 M.O.s")
      return
    end if
  end if
!
!                       MECI Microstates read in
!
  if (myword(allkey, " MICROS")) &
    write (*,'(" *  MICROS=N -", i4, " MICROSTATES TO BE SUPPLIED FOR C.I.")') &
    Int (reada (keywrd, Index (keywrd, " MICROS")))
!
!                       Fractionally occupied degenerate Open shell
!
  if (myword(allkey, " OPEN")) then
    i = Index (keywrd, " OPEN")
    j = Index (keywrd(i:i+10), ",") + i   - 1
    ilevel = Nint (reada (keywrd, j))
    ielec = Nint (reada (keywrd, Index (keywrd, " OPEN")+6)) 
    write (*,'(" *  OPEN(N,N)- THERE ARE", i2, " ELECTRONS IN", i2, " LEVELS")') &
     ielec, ilevel
  end if
!
!                       Magnetic component of spin
!
  if (myword(allkey, " MS=")) &
    write (*,'(" *  MS=        - IN MECI, MAGNETIC COMPONENT OF SPIN =", i3)') &
    Nint (reada (keywrd, Index (keywrd, " MS=")))
!
!   Select root of C.I. matrix
!
  if (myword(allkey, " ROOT")) &
    write (*,'(" *  ROOT       - IN A C.I. CALCULATION, ROOT", i2, &
   & " TO BE OPTIMIZED.")') Nint (reada (keywrd, Index (keywrd, " ROOT")))
!**********************************************************************
!
!   SOLVATION KEYWORDS
!

  if (myword(allkey, " EPS=") .or. (index(keywrd," PKA") .ne. 0)) then
    if (Index(keywrd," EPS=CRSDEF") /= 0) then
      fepsi = 1.d0
    else
      if (index(keywrd," PKA") .ne. 0) then
        epsi = 78.4d0
      else
        epsi = reada (keywrd, Index (keywrd, "EPS="))
      end if
      write (*, "(a, f6.2, a)")" *  EPS=", epsi," - USE ANDREAS KLAMT'S COSMO IMPLICIT SOLVATION MODEL"
      fepsi = (epsi-1.d0) / (epsi+0.5d0)
    end if
    if (epsi < 1.d0) then
      write(*,"(a)") " *  "
      write(*,"(a)") " *  The lower bound of a dielectric is 1.00 = vacuum"
      write(*,"(2a)")" *  Keywords: ",keywrd(:len_trim(keywrd))
      write(*,"(a)") " *  "
      call mopend("Nonsense value of dielectric constant supplied")
      return
    end if
  else if (myword(allkey, " COSMO")) then
    epsi = 999.d0
    write (*, "(a, f6.2, a)")" *  EPS=", epsi," - USE ANDREAS KLAMT'S COSMO IMPLICIT SOLVATION MODEL"
    fepsi = (epsi-1.d0) / (epsi+0.5d0)
  else if (numat < 100) then  !  Limit surface area calculation to 100 atoms   - saves time
    epsi = 78.4d0   !  For surface area and volume only
    fepsi = (epsi-1.d0) / (epsi+0.5d0)
  end if
  if (myword(allkey, " DIPL"))  write (*,'(" *  DIPL=", f7.3)') reada (keywrd, Index (keywrd, " DIPL"))
  if (myword(allkey, " RSOLV")) write (*,'(" *  RSOLV=", f7.3)') reada (keywrd, Index (keywrd, " RSOLV"))
  if (myword(allkey, " DELSC")) write (*,'(" *  DELSC=", f7.3)') reada (keywrd, Index (keywrd, " DELSC"))
  if (myword(allkey, " DISEX")) write (*,'(" *  DISEX=", f7.3)') reada (keywrd, Index (keywrd, " DISEX"))
  if (myword(allkey, "N**2"))   write (*,'(" *  n**2 =", f7.3, " for COSMO-CI")') reada (keywrd, 4+Index (keywrd, "N**2"))
  if (myword(allkey, " NSPA")) then
    i = Nint (reada (keywrd, Index (keywrd, " NSPA")))
    write (*,'(" *  NSPA=", i5)') i
    nspa = i
  else
    nspa = 42
  end if
  if (myword(allkey, " ROTX"))  write (*,'(" *  ROTX=", i5)') Nint (reada (keywrd, Index (keywrd, " ROTX")))
!**********************************************************************
!
!                       Geometric quantities
!
!                       Grids, Reaction paths, etc.
!
 i = 0
 if (myword(allkey, " STEP1"))  then
   write (*,'(" *  STEP1      - FIRST STEP-SIZE IN GRID =", f7.2)') &
   reada (keywrd, Index (keywrd, "STEP1")+6)
   if (index(keywrd, " POINT1") == 0) then
     write (*,'("*",/," *             - **** KEYWORD POINT1 MISSING ****",/,"*")') 
     call mopend("KEYWORD POINT1 MISSING")
   end if
   i = 2
 end if
 if (myword(allkey, " STEP2"))  then
   write (*,'(" *  STEP2      - SECOND STEP-SIZE IN GRID =", f7.2)') &
   reada (keywrd, Index (keywrd, "STEP2")+6)
   if (index(keywrd, " POINT2") == 0) then
     write (*,'("*",/," *             - **** KEYWORD POINT2 MISSING ****",/,"*")') 
     call mopend("KEYWORD POINT2 MISSING")
   end if
 end if
 if (myword(allkey, " STEP="))  then
   write (*,'(" *  STEP       - STEP-SIZE IN PATH =", f7.2)') &
   reada (keywrd, Index (keywrd, " STEP") + 4)
   if (index(keywrd, " POINT") == 0) call mopend ("KEYWORD POINT MISSING")
   i = 1
 end if
 if (i > 0) then
!
! Check that the appropriate flags are set
!
   k = 0
   do j = 1, natoms
     do l = 1, 3
       if (lopt(l,j) < 0) k = k + 1
     end do
   end do
   if(i /= k) then
     if (index(keywrd," 0SCF") == 0) &
       call mopend("The number of optimization flags set to '-1' does not match the keywords used")
   end if
 end if
 if (myword(allkey, " MERS"))then
   j = index(keywrd," MERS")
   mers = 0
   k = 0
   i = Index (keywrd(j + 1:), " ") + j
   do l = 1, 3
     j = j + k
     if (l > 1 .and. k == 0) exit
     mers(l) = Nint (reada (keywrd(j:), 1))
     k = Index (keywrd(j:i), ",")
   end do
   i = mers(1)
   do l = 2, 3
     if (mers(l) == 0) exit
     i = i*mers(l)
   end do
   write (*,'(" *  MERS=N     - NUMBER OF FUNDAMENTAL UNIT CELLS USED:", i3)') i
 else
   mers = 0
 end if
 if (myword(allkey, " POINT1")) then
   write (*,'(" *  POINT1     - NUMBER OF ROWS IN GRID =", i3)') &
   Nint (reada (keywrd, Index (keywrd, "POINT1")+7))
   if (index(keywrd, " STEP1") == 0) then
     write (*,'("*",/," *             - **** KEYWORD STEP1 MISSING ****",/,"*")') 
     call mopend("KEYWORD STEP1 MISSING")
   end if
 end if
 if (myword(allkey, " POINT2")) then
   write (*,'(" *  POINT2     - NUMBER OF COLUMNS IN GRID =", i3)') &
   Nint (reada (keywrd, Index (keywrd, "POINT2")+7))
   if (index(keywrd, " STEP2") == 0) then
     write (*,'("*",/," *             - **** KEYWORD STEP2 MISSING ****",/,"*")') 
     call mopend("KEYWORD STEP2 MISSING")
   end if
 end if
 if (myword(allkey, " POINT")) then
   write (*,'(" *  POINT      - NUMBER OF POINTS IN PATH =", i3)') &
   Nint (reada (keywrd, Index (keywrd, "POINT")+6))
   if (index(keywrd, " STEP") == 0) then
     write (*,'("*",/," *             - **** KEYWORD STEP MISSING ****",/,"*")') 
     call mopend("KEYWORD STEP MISSING")
   end if
 end if
 if (myword(allkey, " KINE"))   write (*,'(" *  KINETIC=   - ", f9.3, " KCAL KINETIC ENERGY ADDED TO DRC")') &
    reada (keywrd, Index (keywrd, " KINE"))
!
!                       Specific intrinsic reaction coordinate
!
 if (myword(allkey, " IRC="))  then
    write (*,'(" *  IRC=N      - INTRINSIC REACTION COORDINATE", i3, " DEFINED")') &
    Nint (reada (keywrd, Index (keywrd, " IRC=")))
  else if (myword(allkey, " IRC")) then
    write (*,'(" *  IRC        - INTRINSIC REACTION COORDINATE CALCULATION")')
  end if
!
!                       DRC (may have half-life)
!
  if (myword(allkey, " DRC=")) then
    write (*,'(" *  DRC=       - HALF-LIFE FOR KINETIC ENERGY LOSS =", f9.2, " *10**(-15) SECONDS")') &
     reada (keywrd, Index (keywrd, " DRC="))
  else if (myword(allkey, " DRC")) then
    write (*,'(" *  DRC        - DYNAMIC REACTION COORDINATE CALCULATION")')
  end if
   mozyme = l_temp
!
!                       External parameters read from file
!
  if ( .not. myword(allkey, " EXTERNAL")) return
  l = 1
  do
    i = Index (keywrd(l:), " EXTERNAL")
    if (i /= 0) then
      l = l + i
      i = l   - 1
      j = Index (keywrd(i:), "=") + i
      if (keywrd(j:j) == """") then
        i = Index (keywrd(j:), """ ") + j   - 1
        allkey(j:i) = " "
        j = j + 1
        i = i   - 1
      else
        i = Index (keywrd(j:), " ") + j   - 1
      end if
10870   format (" *  EXTERNAL   - USE ATOMIC PARAMETERS FROM THE ", &
     & "FOLLOWING FILE",/15 x, a)
      write (*, 10870) keywrd(j:i)
    else
      exit
    end if
  end do
  return
end subroutine wrtcon
subroutine wrtout (allkey)
  use molkst_C, only : keywrd, mozyme, maxtxt, line, prt_coords, prt_gradients, prt_cart, prt_charges, prt_pops, &
    prt_topo
  use chanel_C, only: log, input_fn
  use myword_I
  use reada_I
  implicit none
  character (len=1000), intent (inout) :: allkey
  integer :: i, j
  intrinsic Index, Nint
  if (myword(allkey, " PRTINT"))  write (*,'(" *  PRTINT     - INTERATOMIC DISTANCES TO BE PRINTED")')
  if (myword(allkey, " PRTCHAR")) write (*,'(" *  PRTCHARGE- PRINT CHARGES IN ARC FILE")')
  if (index(allkey, " AUX") > 0) then
    i = index(allkey, " AUX(") + 1
    if (i > 1) then
      j = index(allkey(i + 4:), ") ") + i + 3
      allkey(i:j) = " "
    else
      i = index(allkey, " AUX")
      allkey(i:i + 3) = " "
    end if      
                                  write (*,'(" *  AUX        - OUTPUT AUXILIARY INFORMATION")')
  end if
  if (index(allkey, " OUTPUT") > 0) then
    i = index(allkey, " OUTPUT(") + 1
    if (i > 1) then
      j = index(allkey(i + 4:), ") ") + i + 3
      line = allkey(i + 6:j - 1)
      prt_coords    = (index(line, "C") /= 0)
      prt_gradients = (index(line, "G") /= 0)
      prt_cart      = (index(line, "X") /= 0)
      prt_charges   = (index(line, "Q") /= 0)
      prt_pops      = (index(line, "P") /= 0) 
      prt_topo      = (index(line, "T") /= 0) 
                      write (*,'(" *  MOPAC      - USE OLD MOPAC CONVENTION FOR FIRST THREE ATOMS")')
      write (*,'(" *  OUTPUT     - REDUCE OUTPUT, BUT PRINT:")')
      if (prt_coords)    write (*,'(" *           C - COORDINATES")')
      if (prt_gradients) write (*,'(" *           G - GRADIENTS")')
      if (prt_cart)      write (*,'(" *           X - CARTESIAN COORDINATES")')
      if (prt_charges)   write (*,'(" *           Q - FRACTIONAL ATOMIC CHARGES")')
      if (prt_pops)      write (*,'(" *           P - ATOMIC ORBITAL POPULATIONS")')
      if (prt_topo)      write (*,'(" *           T - TOPOGRAPHY (ATOM CONNECTIVITY)")')
    else
      write (*,'(" *  OUTPUT     - MINIMIZE OUTPUT")')
      i = index(allkey, " OUTPUT")
      prt_coords = .FALSE.
      prt_gradients = .FALSE.
      prt_cart = .FALSE.
      prt_charges = .FALSE.
      prt_pops = .FALSE.
      prt_topo = .FALSE.
    end if
    if (.not. myword(allkey, " OUTPUT")) return ! An impossible option
  else
    prt_coords = .TRUE.
    prt_gradients = .TRUE.
    prt_cart = .TRUE.
    prt_charges = .TRUE.
    prt_pops = .TRUE. 
    prt_topo = .TRUE.
  end if
  if (myword(allkey, " MOPAC"))   write (*,'(" *  MOPAC      - USE OLD MOPAC CONVENTION FOR FIRST THREE ATOMS")')
  if (myword(allkey, " NOINT"))   write (*,'(" *  NOINTER    - THIS KEYWORD HAS BEEN REPLACED BY PRTINT")')
  if (myword(allkey, " ISOTOPE")) write (*,'(" *  ISOTOPE    - FORCE MATRIX WRITTEN TO DISK (CHAN. 9 )")')
  if (myword(allkey, " DENOUT"))  write (*,'(" *  DENOUT     - DENSITY MATRIX OUTPUT ON CHANNEL 10")')
  if (myword(allkey, " OLDENS"))  write (*,'(" *  OLDENS     - INITIAL DENSITY MATRIX READ OF DISK")')
  if (myword(allkey, " CIOSCI"))  write (*,'(" *  CIOSCI     - PRINT WORKING IN SUBROUTINE CIOSCI")')
  if (myword(allkey, " ENPART"))  write (*,'(" *  ENPART     - ENERGY TO BE PARTITIONED INTO COMPONENTS")')
  if (myword(allkey, " NOXYZ"))   write (*,'(" *  NOXYZ      - CARTESIAN COORDINATES NOT TO BE PRINTED")')
  if (myword(allkey, " PRTXYZ"))  write (*,'(" *  PRTXYZ     - PRINT CARTESIAN COORDINATES")')
  if (myword(allkey, " NOTXT"))   write (*,'(" *  NOTXT      - DO NOT PRINT TEXT ASSOCIATED WITH AN ATOM")')
  if (myword(allkey, " GRAD"))    write (*,'(" *  GRADIENTS- ALL GRADIENTS TO BE PRINTED")')
  if (myword(allkey, " MINI"))    write (*,'(" *  MINI       - REDUCE OUTPUT BY ONLY PRINTING FLAGGED ATOMS")')
  if (mozyme) then
    if (myword(allkey, " VEC"))     write (*,'(" *  VECTORS    - FINAL MOLECULAR ORBITALS TO BE PRINTED")')
  else
    if (myword(allkey, " VEC"))     write (*,'(" *  VECTORS    - FINAL EIGENVECTORS TO BE PRINTED")')
  end if
  if (myword(allkey, " EIGEN"))   write (*,'(" *  EIGEN      - CONVERT LMO''s INTO EIGENVECTORS")')
  if (myword(allkey, " CARTAB"))  write (*,'(" *  CARTAB     - PRINT ALL THE CHARACTER TABLES")')
  if (myword(allkey, " CHARST"))  write (*,'(" *  CHARST     - PRINT DETAILS OF SUBROUTINE CHARST")')
  if (myword(allkey, " DCART"))   write (*,'(" *  DCART      - PRINT DETAILS OF SUBROUTINE DCART")')
  if (myword(allkey, " DERI1"))   write (*,'(" *  DERI1      - PRINT DETAILS OF SUBROUTINE DERI1")')
  if (myword(allkey, " DERI2"))   write (*,'(" *  DERI2      - PRINT DETAILS OF SUBROUTINE DERI2 ")')
  if (myword(allkey, " DERITR"))  write (*,'(" *  DERITR     - PRINT DETAILS OF SUBROUTINE DERITR")')
  if (myword(allkey, " DERIV"))   write (*,'(" *  DERIV      - PRINT DETAILS OF SUBROUTINE DERIV")')
  if (myword(allkey, " DERNVO"))  write (*,'(" *  DERNVO     - PRINT DETAILS OF SUBROUTINE DERNVO")')
  if (myword(allkey, " MAKVEC"))  write (*,'(" *  MAKVEC     - PRINT DETAILS OF SUBROUTINE MAKVEC")')
  if (myword(allkey, " DIIS"))    write (*,'(" *  DIIS       - PRINT DETAILS OF SUBROUTINE DIIS  ")')
  if (myword(allkey, " FLEPO"))   write (*,'(" *  FLEPO      - PRINT DETAILS OF SUBROUTINE FLEPO ")')
  if (myword(allkey, " FMAT"))    write (*,'(" *  FMAT       - PRINT DETAILS OF SUBROUTINE FMAT  ")')
  if (myword(allkey, " DFORCE"))  write (*,'(" *  DFORCE     - PRINT FORCE MATRIX OVER CARTESIAN COORDINATES")')
  if (myword(allkey, " HCORE"))   write (*,'(" *  HCORE      - PRINT DETAILS OF SUBROUTINE HCORE ")')
  if (myword(allkey, " MOLDAT"))  write (*,'(" *  MOLDAT     - PRINT DETAILS OF SUBROUTINE MOLDAT")')
  if (myword(allkey, " FREQCY"))  write (*,'(" *  FREQCY     - PRINT DETAILS OF SUBROUTINE FREQCY")')
  if (myword(allkey, " ITER"))    write (*,'(" *  ITER       - PRINT DETAILS OF SUBROUTINE ITER  ")')
  if (myword(allkey, " LINMIN"))  write (*,'(" *  LINMIN     - PRINT DETAILS OF SUBROUTINE LINMIN")')
  if (myword(allkey, " MOLSYM"))  write (*,'(" *  MOLSYM     - PRINT DETAILS OF SUBROUTINE MOLSYM")')
  if (myword(allkey, " RAMA"))    write (*,'(" *  RAMA       - PRINT RAMACHANDRA ANGLES FOR PROTEIN RESIDUES")')
  if (myword(allkey, " SYMOIR"))  write (*,'(" *  SYMOIR     - PRINT DETAILS OF SUBROUTINE SYMOIR")')
  if (myword(allkey, " GROUP"))   write (*,'(" *  GROUP      - PRINT DETAILS OF SUBROUTINE GROUP ")')
  if (myword(allkey, " SYMTRZ"))  write (*,'(" *  SYMTRZ     - PRINT DETAILS OF SUBROUTINE SYMTRZ")')
  if (myword(allkey, " DENS"))    write (*,'(" *  DENSITY    - FINAL DENSITY MATRIX TO BE PRINTED")')
  if (myword(allkey, " SPIN"))    write (*,'(" *  SPIN       - FINAL UHF SPIN MATRIX TO BE PRINTED")')
  if (myword(allkey, " PRNT"))    write (*,'(" *  PRNT       - EXTRA PRINTING IN EF ROUTINE")')
  if (myword(allkey, " DISP"))    write (*,'(" *  DISP       - PRINT DISPERSION AND HYDROGEN BOND ENERGIES")')
  if (myword(allkey, " DEP ")) then
10410 format (" *  DEP        - OUTPUT FORTRAN CODE FOR BLOCK-DATA") 
10411 format (" *             THIS KEYWORD CANNOT BE USED IN THIS VERSION OF MOPAC.")
    write (*, 10410)
    write (*, 10411)
  end if
  if (myword(allkey, " TIMES"))  write (*,'(" *  TIMES      - TIMES OF VARIOUS STAGES TO BE PRINTED")')
  if (mozyme)then
    if (myword(allkey, " BONDS"))  write (*,'(" *  BONDS      - NON-HYDROGEN BOND-ORDERS TO BE PRINTED")')
    if (myword(allkey, " ALLBO"))  write (*,'(" *  ALLBONDS   - ALL SIGNIFICANT BOND-ORDERS TO BE PRINTED")')
  else
    if (myword(allkey, " BONDS"))  write (*,'(" *  BONDS      - FINAL BOND-ORDER MATRIX TO BE PRINTED")')
  end if
  if (myword(allkey, " SYBYL"))  write (*,'(" *  SYBYL      - OUTPUT SYBYL FILE")')
  if (myword(allkey, " FOCK"))   write (*,'(" *  FOCK       - LAST FOCK MATRIX TO BE PRINTED")')
  if (myword(allkey, " LARGE"))  write (*,'(" *  LARGE      - EXPANDED OUTPUT TO BE PRINTED")')
  log = .false.
  if (myword(allkey, " LOG") .or. &
    (index(keywrd," 0SCF") /= 0 .and. index(keywrd," OLDGEO") /= 0 .and. index(keywrd," PDBOUT") /= 0) .or. &
    index(keywrd, " ADD-H") /= 0 .or. index(keywrd, " HEADER") /= 0) then
    log = .true.
    if (index(keywrd, " LOG") /= 0) write (*,'(" *  LOG        - GENERATE A LOG FILE")')
  end if
  if (myword(allkey, " NOLOG")) then
    write (*,'(" *  NOLOG      - SUPPRESS LOG FILE")')
  end if
  if (myword(allkey, " HTML"))   then
    if ((index(keywrd, " STEP=") == 0 .or. index(keywrd, " POINT") == 0) .and. &
    index(keywrd, " IRC") + index(keywrd, " DRC") == 0) then      
      write (*,'(" *  HTML       - WRITE HTML SCRIPT TO READ PDB FILE USING JSMOL")')
    else
      write (*,'(" *  HTML       - WRITE HTML SCRIPT TO GENERATE ANIMATION USING JSMOL")')
      if (index(keywrd, " PDBOUT") /= 0) &
        call mopend("When STEP is present, PDBOUT must not be present.")        
    end if
!
!  Check that file-name is okay
!

  end if
  if (myword(allkey, " AIGOUT")) write (*,'(" *  AIGOUT     - IN ARC FILE, INCLUDE AB-INITIO GEOMETRY")')
  if (myword(allkey, " GRAP"))   write (*,'(" *  GRAPH      - GENERATE FILE FOR GRAPHICS")')
  if (myword(allkey, " 1ELEC"))  write (*,'(" *  1ELECTRON- FINAL ONE-ELECTRON MATRIX TO BE PRINTED")')
  if (myword(allkey, " INTERP")) write (*,'(" *  INTERP     - PRINT DETAILS OF CAMP-KING CONVERGER")')
  if (myword(allkey, " HTML"))   write (*,'(" *  HTML       - WRITE HTML SCRIPT TO READ PDB FILE USING JSMOL")')
  if (myword(allkey, " SUPER"))  write (*,'(" *  SUPER      - PRINT SUPERDELOCALIZABILITIES")')
  if (myword(allkey, " ALLVEC")) write (*,'(" *  ALLVEC     - PRINT ALL EIGENVECTORS")')
  if (myword(allkey, " H-PRIO")) write (*,'(" *  H-PRIOR    - HEAT OF FORMATION TAKES PRIORITY IN DRC")')
  if (myword(allkey, " X-PRIO")) write (*,'(" *  X-PRIOR    - GEOMETRY CHANGES TAKE PRIORITY IN DRC")')
  if (myword(allkey, " T-PRIO")) write (*,'(" *  T-PRIOR    - TIME TAKES PRIORITY IN DRC")')
  if (myword(allkey, " COMPFG")) write (*,'(" *  COMPFG     - PRINT HEAT OF FORMATION CALC''D IN COMPFG")')
  if (myword(allkey, " MECI"))   write (*,'(" *  MECI       - M.E.C.I. WORKING TO BE PRINTED")')
  if (myword(allkey, " HESSIAN"))write (*,'(" *  HESSIAN    - WRITE OUT HESSIAN FROM GEOMERY OPTIMIZATION")')
  if (myword(allkey, " LOCAL"))  write (*,'(" *  LOCALIZE   - LOCALIZED ORBITALS TO BE PRINTED")')
  if (myword(allkey, " MULLIK")) write (*,'(" *  MULLIK     - THE MULLIKEN ANALYSIS TO BE PERFORMED")')
  if (myword(allkey, " PI "))    write (*,'(" *  PI         - BONDS MATRIX, SPLIT INTO SIGMA-PI-DELL", &
   & " COMPONENTS, TO BE PRINTED")')
  if (myword(allkey, " ECHO"))   write (*,'(" *  ECHO       - ALL INPUT DATA TO BE ECHOED BEFORE RUN")')
  if (myword(allkey, " DEBUG ")) write (*,'(" *  DEBUG      - DEBUG OPTION TURNED ON")')
  if (myword(allkey, " POTWRT")) write (*,'(" *  POTWRT     - WRITE OUT ELECTRIC POT. DATA TO FILE 21")')
  if (myword(allkey, " PMEP"))   write (*,'(" *  PMEP       - GENERATE THE WANG-FORD ELECTROSTATIC POTENTIAL")')
  if (myword(allkey, " PRTMEP")) write (*,'(" *  PRTMEP     - MEP CONTOUR DATA OUTPUT TO <FILENAME>.mep")')
  if (myword(allkey, " MINMEP")) write (*,'(" *  MINMEP     - PRINT MEP MINIMA IN THE PLANE DEFINED")')
  if (myword(allkey, " PL"))     write (*,'(" *  PL         - MONITOR CONVERGENCE IN DENSITY MATRIX")')
  if (myword(allkey, " PKA"))    write (*,'(" *  PKA        - PRINT pKa FOR MOST ACIDIC HYDROGEN")')
  if (myword(allkey, " POPS"))   write (*,'(" *  POPS       - PRINT SCF ATOMIC ORBITAL POPULATIONS")')
  if (myword(allkey, " BCC"))    write (*,'(" *  BCC        - THE SYSTEM IS BODY-CENTERED CUBIC")')
  if (myword(allkey, " Z="))     write (*,'(" *  Z=",I2,"     - NUMBER OF FORMULA UNITS IN UNIT CELL")') &
   Nint(Reada(keywrd,Index(keywrd," Z=") + 2))
  if (myword(allkey, " EIGS"))   write (*,'(" *  EIGS       - PRINT ALL EIGENVALUES IN ITER")')
  if (myword(allkey, " HYPERF")) write (*,'(" *  HYPERFINE- HYPERFINE COUPLING CONSTANTS TO BE", " PRINTED")')
  if (myword(allkey, " THERMO")) write (*,'(" *  THERMO     - THERMODYNAMIC QUANTITIES TO BE CALCULATED")')
  if (myword(allkey, " K=")) then
    i = Index (keywrd, " K=")
    write (*,' (" *   K=        - BRILLOUIN ZONE STRUCTURE TO BE CALCULATED")')
    write (*,'(" *             STEP SIZE IN SAMPLING ZONE:", f8.4)') reada (keywrd, i)
    i = Index (keywrd(i:), ",") + i
    write (*,'(" *             NO. OF ATOMS IN FUNDAMENTAL UNIT CELL:", i6)') Nint (reada (keywrd, i))
  end if
  return
end subroutine wrtout
subroutine wrtwor (allkey)
  use molkst_C, only: tleft, tdump, keywrd, natoms, id, numat, maxtxt
  use common_arrays_C, only : labels
  use myword_I
  use reada_I
  implicit none
  character (len=1000), intent (inout) :: allkey
  character (len=20) :: spaces
  character :: ch
  character (len=7) :: chrono
  integer :: i, ii, j
  double precision :: time
  intrinsic Index, Min, Nint, Max
  spaces = " "
  if (myword(allkey, " EIGINV"))     write (*,'(" *  EIGINV     - USE HESSIAN EIGENVALUE REVERSION IN EF")')
  if (myword(allkey, " NONR"))       write (*,'(" *  NONR       - DO NOT USE NEWTON-RAPHSON STEP IN EF")')
  if (myword(allkey, " SNAP"))       write (*,'(" *  SNAP       - INCREASE PRECISION OF SYMMETRY ANGLES")')
  if (myword(allkey, " PULAY"))      write (*,'(" *  PULAY      - PULAY''S METHOD TO BE USED IN SCF")')
  if (myword(allkey, " CAMP"))       write (*,'(" *  CAMP,KING- THE CAMP-KING CONVERGER TO BE USED")')
  if (myword(allkey, " KING"))       write (*,'(" *  CAMP,KING- THE CAMP-KING CONVERGER TO BE USED")')
  if (myword(allkey, " LET"))        write (*,'(" *  LET        - OVERRIDE SOME SAFETY CHECKS")')
  if (myword(allkey, " OLDGEO"))     write (*,'(" *  OLDGEO     - PREVIOUS GEOMETRY TO BE USED")')
  if (myword(allkey, " OLDFPC"))     write (*,'(" *  OLDFPC     - OLD FUNDAMENTAL PHYSICAL CONSTANTS TO BE USED")')
  if (myword(allkey, " OLD_HESS"))   write (*,'(" *  OLD_HESS   - USE THE OLD HESSIAN MATRIX")')
  if (myword(allkey, " PREC"))       write (*,'(" *  PRECISE    - CRITERIA TO BE INCREASED BY 100 TIMES")')
  if (myword(allkey, " NOANCI"))     write (*,'(" *  NOANCI     - DO NOT USE ANALYTICAL C.I. DERIVATIVES")')
  if (myword(allkey, " DFP"))        write (*,'(" *  DFP        - USE DAVIDON FLETCHER POWELL OPTIMIZER")')
  if (myword(allkey, " XYZ"))        write (*,'(" *  XYZ        - CARTESIAN COORDINATE SYSTEM TO BE USED")')
  if (myword(allkey, " RESTART"))    write (*,'(" *  RESTART    - CALCULATION RESTARTED")')
  if (myword(allkey, " RSCAL"))      write (*,'(" *  RSCAL      - SCALE P-RFO STEP IN EF TO TRUST RADIUS")')
  if (myword(allkey, " SYMAVG"))     write (*,'(" *  SYMAVG     - AVERAGE SYMMETRY EQUIVALENT ESP CHARGES")')
  if (myword(allkey, " STO3G"))      write (*,'(" *  STO3G      - DEORTHOGONALIZE ORBITALS IN STO-3G BASIS")')
  if (myword(allkey, " FORCE "))     write (*,'(" *  FORCE      - FORCE CALCULATION SPECIFIED")')
  if (myword(allkey, " FORCETS"))    write (*,'(" *  FORCETS    - VERIFY THAT TRANSITION STATE IS GENUINE")')
  if (myword(allkey, " BFGS"))       write (*,'(" *  BFGS       - USE THE BFGS GEOMETRY OPTIMIZER")')
  if (myword(allkey, " LBFGS"))      write (*,'(" *  LBFGS      - USE THE LBFGS GEOMETRY OPTIMIZER")')
  if (myword(allkey, " EF"))         write (*,'(" *  EF         - USE EF ROUTINE FOR MINIMUM SEARCH")')
  if (myword(allkey, " TS"))         write (*,'(" *  TS         - USE EF ROUTINE FOR TS SEARCH")')
  if (myword(allkey, " LOCATE-TS"))  write (*,'(" *  LOCATE-TS  - LOCATE A TRANSITION STATE IN A PROTEIN MECHANISM")')
  if (myword(allkey, " REFINE-TS"))  write (*,'(" *  REFINE-TS  - REFINE A TRANSITION STATE IN A PROTEIN MECHANISM")')
  if (myword(allkey, " NOSYM"))      write (*,'(" *  NOSYM      - POINT-GROUP SYMMETRY SET TO C1")')
  if (myword(allkey, " SMOOTH"))     write (*,'(" *  SMOOTH     - IN A GRID CALCULATION, REMOVE COMPUTATIONAL ARTIFACTS")')
  if (myword(allkey, " AUTOSYM"))    write (*,'(" *  AUTOSYM    - SYMMETRY TO BE IMPOSED AUTOMATICALLY")')
  if (myword(allkey, " 1SCF "))      write (*,'(" *  1SCF       - DO 1 SCF AND THEN STOP ")')
  if (myword(allkey, " SIGMA"))      write (*,'(" *  SIGMA      - GEOMETRY TO BE OPTIMIZED USING SIGMA.")')
  if (myword(allkey, " NLLSQ"))      write (*,'(" *  NLLSQ      - GRADIENTS TO BE MINIMIZED USING NLLSQ.")')
  if (myword(allkey, " SADDLE"))     write (*,'(" *  SADDLE     - TRANSITION STATE TO BE OPTIMIZED")')
  if (myword(allkey, " 0SCF"))       write (*,'(" *  0SCF       - AFTER READING AND PRINTING DATA, STOP")')
  if (myword(allkey, " QPMEP "))     write (*,'(" *  QPMEP      - CHARGES DERIVED FROM WANG-FORD TYPE AM1", " MEP")')
  if (myword(allkey, " NEWGEO"))     write (*,'(" *  NEWGEO     - USE BOTH INT AND CARTESIAN COORDINATES")')
  if (myword(allkey, " BIGSCF"))     write (*,'(" *  BIGSCF     - DO INITIAL FULL SCF WHEN RESTARTING JOB")')
  if (myword(allkey, " WILLIA"))     write (*,'(" *  WILLIAMS   - USE WILLIAMS SURFACE")')
  if (myword(allkey, " INT "))       write (*,'(" *  INT        - INTERNAL COORDINATE SYSTEM TO BE USED")')
  if (myword(allkey, " PECI"))       write (*,'(" *  PECI       - SINGLE AND PAIRED ELECTRON", " EXCITATIONS ONLY")')
!
! Number of threads
!
  if (myword(allkey, " THREADS")) then
    i = Index (keywrd, " THREADS")
    i = nint(reada(keywrd, i))
    i = max(i,1)
    if (i == 1) then
      write (*,'(" *  THREADS=1  - MULTI-THREADING NOT USED")')
    else if (i < 10) then
      write (*,'(" *  THREADS    - USE A MAXIMUM OF", i2, " THREADS")') i
    else
      write (*,'(" *  THREADS    - USE A MAXIMUM OF", i3, " THREADS")') i
    end if
  end if
  !                   Times
  chrono = "SECONDS"
  time = 1.0d0
  tleft = 172800.d0 ! two days
  if (myword(allkey, " T=")) then
    i = Index (keywrd, " T=")
    tleft = reada (keywrd, i)
    do j = i + 3, 1000
      if (j == 1000 .or. keywrd(j+1:j+1) == " ") go to 1000
    end do
1010  if (tleft < 99999.9d0) then
      write (*,'(" *  T=         - A TIME OF", f9.1, " ", a7, " REQUESTED")') tleft, chrono
    else 
      write (*,'(" *  T=         - A TIME OF", g9.1, " ", a7, " REQUESTED")') tleft, chrono
    end if
!
!  Limit time to 9,999,999 seconds =115.74 days.
!
    tleft = Min (1.d7-1.d0, tleft*time)
    go to 1020
1000  ch = keywrd(j:j)
    if (ch == "M") then
      chrono = "MINUTES"
      time = 60.0d0
    end if
    if (ch == "H") then
      chrono = "HOURS"
      time = 3600.0d0
    end if
    if (ch == "D") then
      chrono = "DAYS"
      time = 86400.0d0
    end if
    if (ch == "W") then
      chrono = "WEEKS"
      time = 7*86400.0d0
    end if
    go to 1010
  else 
    write (*,'(" *  T=         - A TIME OF", f9.1, " ", a7, " REQUESTED")') tleft, "SECONDS"
  end if
1020 time = 1.0d0
  chrono = "SECONDS"
  tdump = 7200.d0
  if (myword(allkey, " DUMP")) then
    i = Index (keywrd, " DUMP")
    tdump = reada (keywrd, i)
    do j = i + 6, 1000
      if (j == 1000 .or. keywrd(j+1:j+1) == " ") go to 1030
    end do
1040  if (tdump < 99999.9d0) then
      write (*,'(" *  DUMP=N     - RESTART FILE WRITTEN EVERY", f10.3, " ", a7)') tdump, chrono
    else
      write (*,'(" *  DUMP=N     - RESTART FILE WRITTEN EVERY", g11.3, " ", a7)') tdump, chrono
    end if
    tdump = tdump * time
    go to 1050
1030  ch = keywrd(j:j)
    if (ch == "M") then
      chrono = "MINUTES"
      time = 60.d0
    end if
    if (ch == "H") then
      chrono = "HOURS"
      time = 3600.d0
    end if
    if (ch == "D") then
      chrono = "DAYS"
      time = 86400.d0
    end if
    if (ch == "W") then
      chrono = "WEEKS"
      time = 7*86400.d0
    end if
    go to 1040
  else 
    write (*, '(" *  DUMP=N     - RESTART FILE WRITTEN EVERY", f10.3, " ", a7)') tdump, "SECONDS"
  end if
1050 if (index(keywrd, " LOCATE-TS") /= 0) then
       tleft = 3.d6
       tdump = 3.d6
     end if
    if (myword(allkey, " CYCLES")) &
    write (*,'(" *  CYCLES=    - DO A MAXIMUM OF", i6, " STEPS")') &
    Nint (reada (keywrd, Index (keywrd, " CYCLES")))
  if (myword(allkey, " BIGCYCLES")) &
    write (*,'(" *  BIGCYCLES= DO A MAXIMUM OF", i6, " BIG STEPS")') &
    Nint (reada (keywrd, Index (keywrd, " BIGCYCLES")))
!
!**********************************************************************
!
!                   Geometric Quantities
!
  if (myword(allkey, " ALT_A")) then
    write (*, '(" *  ALT_A=", a1, "    - USE ALTERNATIVE ATOMS WHERE APPROPRIATE")') keywrd (Index (keywrd, " ALT_A")+7: &
     & Index (keywrd, " ALT_A")+7)
    end if
  if (myword(allkey, " ALT_R")) then
    write (*, '(" *  ALT_R=", a1, "    - USE ALTERNATIVE RESIDUES WHERE APPROPRIATE")') keywrd (Index (keywrd, " ALT_R")+7: &
     & Index (keywrd, " ALT_R")+7)
    end if
  if (myword(allkey, " VDW")) then
    i = Index (keywrd, " VDW")
    j = Index (keywrd(i:), ")") + i
    write (*, '(" *  ", a, "  USER-DEFINED ATOMIC OR VAN DER WAALS RADII")') keywrd(i + 1:j)
  end if
  if (myword(allkey, " METAL")) then
    i = Index (keywrd, " METAL")
    j = Index (keywrd(i:), ")") + i   - 2
    ii = max(15 + i   - j, 1)
    write (*, '(" *  ", a, "- ELEMENTS DEFINED AS FULLY IONIZED METALS FOR LEWIS STRUCTURE")') &
      keywrd(i + 7:j)//spaces(:ii)
  end if
  if (numat /= natoms .and. maxtxt == 26) then
    if (index(keywrd, " RESEQ") + index(keywrd, " ADD-H") + index(keywrd, " SITE") /= 0) then
      call mopend("DUMMY ATOMS CANNOT BE PRESENT WHEN GEOMETRY IS BEING CHANGED")
      if (index(keywrd, " RESEQ") /= 0) write(*,'(10x,a)')"Geometry is changed by keywrd: RESEQ"
      if (index(keywrd, " ADD-H") /= 0) write(*,'(10x,a)')"Geometry is changed by keywrd: ADD-H"
      if (index(keywrd, " SITE")  /= 0) write(*,'(10x,a)')"Geometry is changed by keywrd: SITE"
      write(*,'(10x,a)')"(A simple way to remove dummy atoms is to add keyword ""XYZ"")"
    end if
  end if
  i = index(allkey, " CVB")
  if (i > 0) then
    if (index(keywrd, " RESEQ") /= 0) then
      call mopend ("CVB cannot be used with RESEQ")
      write(*,*)
    else if (index(keywrd, " ADD-H") /= 0) then
      do j = 1, natoms
        if (labels(j) == 1) exit 
      end do
      do ii = j + 1, natoms
        if (labels(ii) /= 1) exit
      end do      
      if (index(keywrd, " NORESEQ") == 0 .or. ii <= natoms - id) then
        if (index(keywrd, " NORESEQ") == 0) &
          call mopend ("CVB cannot be used with ADD-H unless NORESEQ is present")
        if (ii <= natoms - id) then
          write(*,'(/2x,a, i6, a)')" Atom:",j, " is a hydrogen atom."
          write(*,'(2x,a, i6, a)')" Atom:",ii, " is a non-hydrogen atom with a higher atom number than the hydrogen atom."
          call mopend ("CVB cannot be used with ADD-H if there are non-hydrogen atoms after hydrogen atoms")
        end if
        write(*,*)
      end if
    end if
    j = index(allkey(i + 2:), ") ") + i + 1
    do i = i + 2, j
      if (allkey(i:i) == " ") allkey(i:i) ="_"
    end do
  end if
  if (myword(allkey, " CVB")) then
      i = Index (keywrd, " CVB")
      j = Index (keywrd(i:), ")") + i
      if (j - i < 24) then
        write (*, '(" *  ", a, "  USER-DEFINED EXPLICIT COVALENT BONDS")') keywrd(i + 1:j)
      else        
         write (*, '(" *  USER-DEFINED EXPLICIT COVALENT BONDS:")')
         write (*, '(" *    ",a)') keywrd(i + 1:j)
      end if        
    end if
  if (myword(allkey, " TRANS=")) then
    write (*, '(" *  TRANS=     - ", i4, " VIBRATIONS ARE TO BE DELETED FROM THE THERMO CALCULATION")') &
     Nint (reada (keywrd, Index (keywrd, " TRANS=")))
  else if (myword(allkey, " TRANS")) then
    write (*, '(" *  TRANS      - THE REACTION VIBRATION TO BE DELETED FROM THE THERMO CALCULATION")')
  end if
  if (myword(allkey, "MEP=")) then
    i = Nint (reada (keywrd, Index (keywrd, "MEP=")))
    if (i == 1) then
      write (*, '(" *  MEP=1      - MEP IN A CUBIC GRID")')
    else
      write (*, '(" *  MEP=2      - MEP IN CONNOLLY SURFACE")')
    end if
  end if
!
!   How hessian matrix is updated.
!
  if (myword(allkey, " IUPD")) then
    ii = Nint (reada (keywrd, Index (keywrd, " IUPD=")))
    if (ii == 0) then
10730   format (" *  IUPD=0     - HESSIAN WILL NOT BE UPDATED")
      write (*, 10730)
    end if
    if (ii == 1) then
10740   format (" *  IUPD=1     - HESSIAN WILL BE UPDATED USING POWELL")
      write (*, 10740)
    end if
    if (ii == 2) then
10750   format (" *  IUPD=2     - HESSIAN WILL BE UPDATED USING BFGS")
      write (*, 10750)
    end if
  end if
!
!   How starting hessian matrix is obtained.
!
  if (myword(allkey, " HESS=")) then
    ii = Nint (reada (keywrd, Index (keywrd, " HESS=")))
    if (ii == 0) then
10760   format (" *  HESS=0     - DIAGONAL HESSIAN USED AS INITIAL GUESS")
      write (*, 10760)
    end if
    if (ii == 1) then
10770   format (" *  HESS=1     - INITIAL HESSIAN WILL BE CALCULATED")
      write (*, 10770)
    end if
    if (ii == 2) then
10780   format (" *  HESS=2     - INITIAL HESSIAN READ FROM DISK")
      write (*, 10780)
    end if
    if (ii == 3) then
10790   format (" *  HESS=3     - INITIAL HESSIAN READ FROM INPUT")
      write (*, 10790)
    end if
  end if
!
!   Which "normal" mode to follow in geometry
!
  if (myword(allkey, " MODE=")) then
10800 format (" *  MODE=      - FOLLOW HESSIAN MODE", i3, " TOWARD TS")
    write (*, 10800) Nint (reada (keywrd, Index (keywrd, "MODE=")))
  end if
  if (myword(allkey, " RECALC")) then
10810 format (" *  RECALC=    - DO", i4, " CYCLES BETWEEN HESSIAN RECALC")
    write (*, 10810) Nint (reada (keywrd, Index (keywrd, "RECALC")))
  end if
  if (myword(allkey, " RMAX")) then
10820 format (" *  RMAX=      - MAX. CALC./PRED. ENERGY STEP IN TS", f7.3)
    write (*, 10820) reada (keywrd, Index (keywrd, " RMAX="))
  end if
  if (myword(allkey, " RMIN")) then
10830 format (" *  RMIN=      - MIN. CALC./PRED. ENERGY STEP IN TS", f7.3)
    write (*, 10830) reada (keywrd, Index (keywrd, " RMIN="))
  end if
  if (myword(allkey, " DDMAX")) then
10840 format (" *  DDMAX=     - MAXIMUM TRUST RADIUS", f7.3, " ANG/RAD")
    write (*, 10840) reada (keywrd, Index (keywrd, " DDMAX="))
  end if
  if (myword(allkey, " DDMIN")) then
10850 format (" *  DDMIN=     - MINIMUM TRUST RADIUS", f7.3, " ANG/RAD")
    write (*, 10850) reada (keywrd, Index (keywrd, " DDMIN="))
  end if
  if (myword(allkey, " DMAX")) then
10860 format (" *  DMAX=      - STARTING TRUST RADIUS", f7.3, " ANG/RAD")
    write (*, 10860) reada (keywrd, Index (keywrd, "DMAX="))
  end if
  if (myword(allkey, " OMIN")) then
10870 format (" *  OMIN=      - MINIMUM EIGENVECTOR OVERLAP IN TS", f7.3)
    write (*, 10870) reada (keywrd, Index (keywrd, "OMIN="))
  end if
  if (myword(allkey, " GNORM")) then
10880 format (" *  GNORM=     - EXIT WHEN GRADIENT NORM DROPS BELOW ", g10.3)
    write (*, 10880) reada (keywrd, Index (keywrd, " GNORM"))
  end if
  if (myword(allkey, " DELTA")) then
10881 format (" *  DELTA=     - EXIT WHEN ENERGY DIFFERENCE DROPS BELOW ", g10.3)
    write (*, 10881) reada (keywrd, Index (keywrd, " DELTA"))
  end if
  if (myword(allkey, " BAR")) then
10890 format (" *  BAR=       - REDUCE BAR LENGTH BY A MAX. OF", f7.2)
    write (*, 10890) reada (keywrd, Index (keywrd, " BAR"))
  end if
  if (myword(allkey, " SLOG=")) then
    i = Index (keywrd, " SLOG=")
10900 format (" *  SLOG=n     - DEFAULT STEP SIZE IN BFGS:", f5.2)
    write (*, 10900) reada (keywrd, i)
  else if (myword(allkey, " SLOG")) then
10910 format (" *  SLOG       - DEFAULT STEP SIZE IN BFGS: 0.25")
    write (*, 10910)
  end if
!**********************************************************************
!
!                      SCF quantities
!
!   SCF Criterion
!
  if (myword(allkey, " RELSCF")) then
      write (*, 10920) reada (keywrd, Index (keywrd, " RELSCF")), 1.d-4
    end if
10920   format (" *  RELSCF     - DEFAULT SCF CRITERION MULTIPLIED BY", g12.3,/, &
       & " *             (DEFAULT SCFCRT=", f8.5, ")")
  if (myword(allkey, " SCFCRT")) then
10930 format (" *  SCFCRT     - DEFAULT SCF CRITERION REPLACED BY", g12.3)
    write (*, 10930) reada (keywrd, Index (keywrd, " SCFCRT"))
  end if
  if (myword(allkey, " SHIFT")) then
10940 format (" *  SHIFT      - A DAMPING FACTOR OF", f8.2, " DEFINED")
    write (*, 10940) reada (keywrd, Index (keywrd, " SHIFT"))
  end if
  if (myword(allkey, " FILL")) then
10950 format (" *  FILL=      - IN RHF CLOSED SHELL, FORCE M.O.", i3, &
   & " TO BE FILLED")
    write (*, 10950) Nint (reada (keywrd, Index (keywrd, " FILL")))
  end if
  if (myword(allkey, " ITRY")) then
10960 format (" *  ITRY=      - DO A MAXIMUM OF", i6, " ITERATIONS FOR SCF")
    write (*, 10960) Nint (reada (keywrd, Index (keywrd, " ITRY")))
  end if
  if (myword(allkey, " DAMP=")) then
10970 format (" *  DAMP=      - DAMP SCF OSCILLATIONS USING DAMP=", f6.2)
    write (*, 10970) reada (keywrd, Index (keywrd, " DAMP=")+6)
  end if
  if (myword(allkey, " CUTOFF=")) then
10981 format (" *  CUTOFF=    - CUTOFF FOR NDDO AND DIPOLE INTERACTIONS =", f6.2, "A")
    write (*, 10981) reada (keywrd, Index (keywrd, " CUTOFF=")+8)
  end if
  if (myword(allkey, " CUTOF2=")) then
10980 format (" *  CUTOF2=    - CUTOFF FOR NDDO INTERACTIONS =", f6.2, "A")
    write (*, 10980) reada (keywrd, Index (keywrd, " CUTOF2=")+8)
  end if
  if (myword(allkey, " CUTOF1=")) then
10990 format (" *  CUTOF1=    - CUTOFF FOR DIPOLE INTERACTIONS =", f6.2, "A")
    write (*, 10990) reada (keywrd, Index (keywrd, " CUTOF1=")+8)
  end if
  if (myword(allkey, " CUTOFP=")) then
11010 format (" *  CUTOFP=    - CUTOFF FOR POLYMER ELECTROSTATICS =", f6.2, "A")
    write (*, 11010) reada (keywrd, Index (keywrd, " CUTOFP=")+8)
  end if
  if (myword(allkey, " NLMO=")) then
11020 format (" *  NLMO=N     - AVERAGE NUMBER OF ATOMS PER LMO =", i4)
    write (*, 11020) Nint (reada (keywrd, Index (keywrd, " NLMO=")))
  end if
  if (myword(allkey, " RELTHR")) then
11030 format (" *  RELTHR     - DEFAULT THRESHOLD FOR LMO's MULTIPLIED BY", &
   & f9.4,/, " *             (DEFAULT WAS 1.D-13)")
    write (*, 11030) Max (reada (keywrd, Index (keywrd, " RELTHR")), 1.d-12)
  end if
  if (myword(allkey, " THRESH")) then
11040 format (" *  THRESH     - DEFAULT THRESHOLD FOR LMO's:", g12.3)
    write (*, 11040) Max (reada (keywrd, Index (keywrd, " THRESH")), 1.d-25)
  end if
if (myword(allkey,' SETGPU=')) then
        i = index(keywrd,' SETGPU=')
     	j = nint(reada(keywrd,i))
    	write (*,'(" *  SETGPU=   - YOUR CALCULATION WILL RUN IN THE GPU NUM. = ",i2)') j           
  endif
  if (myword(allkey,' CPUTHREADS=')) then 
     i = index(keywrd,' CPUTHREADS=')
     j = nint(reada(keywrd,i))
     write (*,'(" *  CPUTHREADS=   - NUM. OF THREADS FOR CPU = ",i2)') j                               
  endif
  
  if (myword(allkey, ' FULLDIAG ')) write (*,'(" *  FULLDIAG   - USE ONLY FULL DIAGONALIZATIONS IN SCF ")')
  return
  end subroutine wrtwor
