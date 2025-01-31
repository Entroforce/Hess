subroutine geochk (computed_charge) 
!***********************************************************************
!                                                                      *
!  GEOCHK DOES SEVERAL VERY DIFFERENT THINGS:
!
!   (A) IT CHECKS THE GEOMETRY TO MAKE SURE IT CONFORMS TO THE LEWIS
!       STRUCTURE CONVENTIONS.
!
!   (B) IT IDENTIFIES ALL IONIZED ATOMS. LATER ON, MAKVEC WILL NEED TO
!       KNOW WHICH ATOMS ARE IONIZED.
!
!   (C) IT CALCULATES THE CHARGE ON THE SYSTEM.
!
!   (D) (OPTIONAL) IT RESEQUENCES THE ATOMS SO THAT ALL ATOMS OF EACH
!       RESIDUE ARE CONTIGUOUS.  IF THIS IS DONE, THE JOB CANNOT BE
!       CONTINUED, SO THE NEW GEOMETRY IS PRINTED OUT AND THE RUN
!       STOPPED.
!
!                                                                      *
!***********************************************************************
    use common_arrays_C, only : geo, coord, nfirst, nlast, breaks, &
       & labels, nat, na, nb, nc, nbonds, ibonds, txtatm, atmass, &
         txtatm1, chains, all_comments, coorda, loc, lopt, l_atom, break_coords
    use MOZYME_C, only : ions, icharges, angles, allres, iz, ib, tyr, allr, &
    iopt, nres, at_res, res_start, Lewis_tot, Lewis_elem, noccupied, nvirtual, &
    odd_h, tyres, start_res
!    
    use molkst_C, only: natoms, numat, nelecs, keywrd, moperr, maxtxt, &
      maxatoms, line, nalpha, nbeta, uhf, nclose, nopen, norbs, id, &
      refkey, ncomments, numat_old, nvar, prt_coords, prt_topo, allkey
    use chanel_C, only: iarc, archive_fn, log
    use mod_atomradii, only: atom_radius_covalent
    use parameters_C, only : ams, natorb, tore, main_group
    use elemts_C, only : elemnt
    use reada_I
    implicit none
    
    integer , intent(inout) :: computed_charge 
    
!
    integer, parameter :: max_sites = 400
    character :: padding*40, txtatm_1*26, txtatm_2*26, tmp*130
    character (len=20), allocatable :: Lewis_formatted(:,:)
    character (len=1), dimension (:), allocatable :: atom_charge
    character :: ion_names(-6:6)*12, charge(max_sites,3)*1, num
    double precision, dimension(:), allocatable :: radius
    logical, save :: debug, let, lres, lreseq, times, opend, charges, l_protein, &
      done, neutral(100), lsite, ter, residues, lbreaks, l_use_old_labels, &
    l_names(max_sites), first
    logical, dimension (:), allocatable :: ioptl
    integer :: i, ibad, ichrge, irefq, ires, ii, jj, m, nfrag, io, kk, &
   & j, jbad, k, l, large, n1, new, alloc_stat, uni_res, mres, near_ions(100), &
     maxtxt_store, nn1, n_new, new_res(max_sites), j2, mbreaks
    integer, dimension(:), allocatable ::  mb    
    integer, save :: numbon(3), num_ions(-6:6)
    integer, dimension(:), allocatable :: nnbonds
    integer, dimension(:,:), allocatable :: iibonds
    double precision :: sum, r_ions(100)
    character :: new_name(max_sites)*3, old_name(max_sites)*3, new_chain(max_sites)*1
    data numbon / 3 * 0 /
    data ion_names / &
      "Less than -5", &
      "Penta-anion ", &
      "Tetra-anion ", &
      "Tri-anion   ", &
      "Di-anion    ", &
      "Anion       ", &
      "(not used)  ", &
      "Cation      ", &
      "Di-cation   ", &
      "Tri-cation  ", &
      "Tetra-cation", &
      "Penta-cation", &
      "More than +5" /
!
    intrinsic Abs, Index, Min, Nint, Allocated
!
    if (allocated(ions))    deallocate (ions)
    if (allocated(iopt))    deallocate(iopt)
    if (allocated(at_res))  deallocate(at_res)
    if (Allocated (iz))     deallocate (iz)
    if (Allocated (ib))     deallocate (ib)
    allocate (ions(maxatoms), iz(maxatoms), ib(maxatoms), mb(maxatoms), at_res(maxatoms), &
            & atom_charge(maxatoms), ioptl(maxatoms), iopt(maxatoms), radius(maxatoms), &
            stat=alloc_stat)
    if (alloc_stat /= 0) then
      call mopend("Failed to allocate arrays in GEOCHK")
      go to 1100
    end if
    l_protein = .false.
    ib(:) = 0
    at_res = 0
    mres = 0
    ions = 0
    j = 0
    lbreaks = (breaks(1) /= -300)
    if (.not. lbreaks) mbreaks = 0
    do i = 1, ncomments
      if (index(all_comments(i), "REMARK   2") == 0) then
        j = j + 1
        all_comments(j) = all_comments(i)
      end if
    end do
    ncomments = j
!
    done = .false.
    n_new = 0
    line = trim(keywrd)
    i = Index (line, " XENO")
    if (i /= 0) then
      do
        if (line(i:i) == "(") exit
        line = line(2:)
      end do
      j = index(line(i:), ") ") 
      if (j == 0) then
        call mopend("Closing parenthesis for XENO not found")
        return
      end if      
      line = ","//line(i + 1: j + i)
      do
        if (line == " ") exit
        do j = 1, 100
          if (line(1:1) == "," .or. line(1:1) == ";" .or. line(1:1) == ")") exit
          line = line(2:)
        end do
        if (line(1:1) == ")") exit
        line = line(2:)
        k = index(line, "=")
        if (k == 0) then
          call mopend("Equals sign (""="") expected in XENO but not found")
          return
        end if  
        l = k + 2
        j = k + 4
        n_new = n_new + 1
        new_chain(n_new) = line(1:1)
        if (new_chain(n_new) >= "0" .and. new_chain(n_new) <= "9") new_chain(n_new) = "A"
        new_res(n_new) = nint(reada(line, 1))
        if (line(l:l) == "," .or. line(l:l) == ";" .or. line(l:l) == ")") then
!
!   A one-letter residue name has been detected.
!
          new_name(n_new) = line(k + 1:k + 1)
        else if (line(j:j) == "," .or. line(j:j) == ";" .or. line(j:j) == ")") then
!
!   A three-letter residue name has been detected.
!
          new_name(n_new) = line(k + 1:k + 3)
        else          
          exit
        end if
      end do
    end if      
    do
      i = Index (keywrd, " xeno")
      if (i == 0) exit
      keywrd(i:i+6) = " XENO=("
    end do
    do i = 1, n_new
      do j = i + 1, n_new
        if (new_chain(i) == new_chain(j) .and. new_res(i) == new_res(j) .and. &
          new_name(i) /= new_name(j)) then
          write(line,'(a,a,i4,a)')" XENO residue ",new_chain(i), new_res(i),&
          " occurs more than once.  Remove extra definition."
          call mopend(trim(line))
          return
        end if
      end do
      if (new_name(i)(2:3) == "  ") then
        do j = 1, 20
          if (tyr(j) == new_name(i)(1:1)) then
            new_name(i) = tyres(j)
            exit
          end if
        end do
      end if
    end do
!
!   Add or remove hydrogen atoms, as necessary.
!
    i =  index(keywrd," SITE=(IONIZE)")
    if (i > 0) then
      line = "SITE=(COO,NH3,ARG(+),SO4,PO4)"
      keywrd = keywrd(:i)//trim(line)//keywrd(i + 14:)
    end if      
    i = index(keywrd," SITE=")
    lsite = .false.
!
!   Identify sites that either should be ionized or should not be ionized.
!
!   These are:
!
!   1:   -COOH
!   2:   -COO(-)
!   3:   -NH3(+)
!   4:   -NH2
!   5:   -Arg(+)-
!   6:   -Arg-
!   7:   -His(+)-
!   8:   -His-
!   9:    SO4(=)
!  10:    PO4(=)
!
    
    i = index(keywrd," SITE=")
    if (i /= 0 .and. index (keywrd, " ADD-H") == 0) then
      call lewis(.true.)
    do
      i = index(keywrd," SITE=")
      if (keywrd(i + 7:i + 7) == ")") exit
      i = index(keywrd," SITE=")
      if (i == 0) exit
      do j = i + 1, len_trim(keywrd)
        if (keywrd(j:j + 1) == ") ") exit
      end do
      keywrd = keywrd(:i)//trim(keywrd(j + 1:))
    end do
    end if
!
!  Store charge, if present
!
    do i = 1, natoms
      atom_charge(i) = txtatm(i)(2:2)
!
!  Prevent atom number being mis-read as a charge.
!
      if (txtatm(i)(2:2) /= "+" .and. txtatm(i)(2:2) /= "-" .and. txtatm(i)(2:2) /= "0") &
      & atom_charge(i) = " "
    end do
!
!    ASSIGN LOGICALS USING KEYWRD
!
    lres = (Index (keywrd, " RESI") + Index (keywrd, " NEWGEO") + Index (keywrd, " RESEQ") /= 0)
    if (.not. lres) lres = (Index (keywrd, " PDBOUT") /= 0 .and. maxtxt /= 26)
    if (.not. lres) lres = (Index (keywrd, " ADD-H") /= 0 .and. Index (keywrd, " NORESEQ") == 0)
    lreseq = (Index (keywrd, " NORESEQ") == 0 .and. Index (keywrd, " RESEQ") /= 0)
    if (lreseq .and. maxtxt /= 26 .and. index(keywrd, "RESID") == 0) then
      write(line,'(a)')"RESEQ only works when the atom labels are in PDB format"
      call mopend(trim(line))
      write(*,'(a,/)')"(Before using RESEQ, run a job using keyword RESIDUES to add PDB atom labels.)"      
      return
    end if
    let = (Index (keywrd, " 0SCF")+Index (keywrd, " LET")+ &
   & Index (keywrd, " RESEQ")+Index (keywrd, " GEO-OK") /= 0)
    times = (Index (keywrd, " TIMES") /= 0)
    if (times) then
      call timer (" START OF GEOCHK")
    end if
    debug = (Index (keywrd, " GEOCHK") /= 0)
    if (Index (keywrd, " CHARGE=") /= 0) then
      irefq = Nint (reada (keywrd, Index (keywrd, " CHARGE=")))
    else
      irefq = 0
    end if
    call extvdw_for_MOZYME (radius, atom_radius_covalent)
    if (moperr) return
!
    if (Index (keywrd, " LARGE") /= 0) then
      large = 1000000
    else
      large = 20
    end if
!
!   WORK OUT WHAT ATOMS ARE BONDED TO EACH OTHER.
!
    call lewis (.true.)
    if (moperr) return   
    allocate (nnbonds(numat), iibonds(15,numat), stat=alloc_stat)
    if (alloc_stat /= 0) then
      call memory_error ("geochn")
      go to 1100
    end if
!
!  Store nbonds and ibonds in case they are modified within this subroutine
!
    nnbonds = nbonds
    iibonds = ibonds
    if (moperr) then
      if (Index (keywrd, " GEO-OK") == 0) then
        write (*,*) " GEOMETRY CONTAINS FAULTS. TO CONTINUE CALCULATION SPECIFY ""GEO-OK"""
        go to 1100
      else
        moperr = .false.
      end if
    end if
!
!   ZERO OUT IONS
!
!
!  FIND THE NITROGEN ATOM OF THE N END OF THE PROTEIN.
!
    do i = 1, numat
      ioptl(i) = .false.
    end do
    call findn1 (n1, ioptl, io)
    if (lreseq) then
!
!   RESEQUENCE THE ATOMS.  WHEN RESEQ IS CALLED, THE SEQUENCE OF
!   ATOMS IS CHANGED.  A CONSEQUENCE OF THIS IS THAT THE SCF CANNOT
!   BE RUN.
!
! First, delete all bonds between ATOMs and HETATMs
!
      do i = 1, numat
        if (txtatm(i)(:4) /= "ATOM") cycle
        do j = 1, nbonds(i)
          k = ibonds(j,i)
          if (txtatm(k)(:4) /= "ATOM") ibonds(j,i) = 0
        end do
        l = 0
        do j = 1, nbonds(i)
          k = ibonds(j,i)
          if (k > 0) then
            l = l + 1
            ibonds(l,i) = k
          end if
        end do
        nbonds(i) = l
      end do
      if (index(keywrd, "RESID") /= 0) txtatm1(:numat) = " "
      new = 0
      iz = -1000
      do
        if (n1 /= 0) call reseq (ioptl, iz, n1, new, io)
        if (moperr) go to 1100
        call findn1 (n1, ioptl, io)
        if (n1 == 0) exit
      end do
      if (new /= numat) then
        do i = 1, numat
          if (nat(i) /= 1 .and. .not. ioptl(i)) then
!
!   IDENTIFY ALL NON-PROTEIN MOLECULES IN THE SYSTEM
!
            call moiety (ioptl, iz, i, new)
          end if
        end do
        do i = 1, numat
          if (nat(i) == 1 .and. .not. ioptl(i)) then
!
!   IDENTIFY ALL HYDROGENS ATTACHED TO RESIDUE-LIKE SPECIES
!
            new = new + 1
            iz(new) = i
          end if
        end do
        if (new /= numat) then
          write (*,*) " THERE IS A FAULT IN RESEQ"
          write (*,'(a,i5)') "  Number of atoms found in data-set:  ", numat
          write (*,'(a,i5)') "  Number of atoms after re-sequencing:", new
          if (new < numat) then
            write(*,'(a)') "Atoms missing (Use original numbering system)"
            ioptl(:numat) = .true.
            do i = 1, new
              ioptl(iz(i)) = .false.              
            end do
            do i = 1, numat
              if (ioptl(i)) write(*,'(i5)')i
            end do         
          else
          end if
          call mopend("THERE IS A FAULT IN RESEQ")
          go to 1100
        end if
!
!  Unconditionally, convert geometry into Cartesian coordinates
!
        geo(:,:numat) = coord(:,:numat)
        do i = 1, numat
          if (na(i) /= 0) exit
        end do
        if (i <= numat) then
          call mopend("SOME COORDINATES WERE IN INTERNAL. THESE HAVE BEEN CHANGED TO CARTESIAN")
          moperr = .false.
        end if
        na = 0
      end if
      l = 1
      do i = 1, numat
        nfirst(i) = l
        j = iz(i)
        mb(j) = i
        do k = 1, 3
          geo(k, i) = coord(k, j)
        end do
        labels(i) = nat(j)
        nlast(i) = nfirst(i) + natorb(labels(i)) - 1
        l = nlast(i) + 1 
      end do
      nat(:numat) = labels(:numat)
      coord(:,:numat) = geo(:,:numat)
      done = .true.
!
!   DUMMY ATOMS ARE EXCLUDED, THEREFORE
!
      natoms = numat
!
!
!   REARRANGE ATOMS TO SUIT THE NEW NUMBERING SYSTEM
!
      do i = 1, numat
        j = iz(i)
        ib(i) = Min (nbonds(j), 4)
        l = Min (nbonds(i), 4)
        do k = 1, l
          ibonds(k+4, i) = mb(ibonds(k, i))
        end do
      end do
      do i = 1, numat
        nbonds(i) = ib(i)
        do k = 1, nbonds(i)
          ibonds(k, i) = ibonds(k+4, iz(i))
        end do
      end do
    end if
    noccupied = 0  
    if (Index (keywrd, " RESEQ") + Index (keywrd, " SITE=") + index(keywrd, " ADD-H") &
      + index(keywrd, " 0SCF") == 0) then
!
!  EXAMINE THE GEOMETRY - IDENTIFY THE LEWIS ELEMENTS (SIGMA BONDS,
!  LONE PAIRS, CATIONS, PI BONDS, ANIONS, OTHER CHARGES)
!
      numbon = 0
      call chklew (mb, numbon, l, large, debug)
      if (moperr) return
      l = 0
      do i = 1, numat
        l = l + Abs (iz(i))
      end do
      if (l /= 0) then         
!
!  THERE ARE IONS.  IDENTIFY THEM.
!
        call chkion (mb, numbon(2), atom_charge)
        if (lreseq) moperr = .false.
      end if
      do i = 1, numat
        ions(i) = nint(tore(nat(i)))
      end do
      do i = 1, Lewis_tot
        if (Lewis_elem(1,i) > 0) then
          noccupied = noccupied + 1
          j = Lewis_elem(1,i)
          if (Lewis_elem(2,i) > 0) then
            k = Lewis_elem(2,i)
            ions(k) = ions(k) - 1! one electron from a bond
            ions(j) = ions(j) - 1! one electron from a bond
          else
            ions(j) = ions(j) - 2! two electrons from a lone pair
          end if
        end if
      end do
      nvirtual = 0
      do i = 1, Lewis_tot
        if (Lewis_elem(2,i) > 0) nvirtual = nvirtual + 1
      end do
    else
      nvar = 0
      do i = 1, numat
        do j = 1, numat_old          
          if (abs(coord(1,i) - coorda(1,j)) > 0.1d0) cycle
          if (abs(coord(2,i) - coorda(2,j)) > 0.1d0) cycle
          if (abs(coord(3,i) - coorda(3,j)) > 0.1d0) cycle
          exit
        end do 
        if (j <= numat_old) then
          do k = 1, 3
            if (lopt(k,j) == 1) then
              nvar = nvar + 1
              loc(1,nvar) = i
              loc(2,nvar) = k
            end if
          end do
        else
          do k = 1, 3
            nvar = nvar + 1
            loc(1,nvar) = i
            loc(2,nvar) = k
          end do          
        end if
      end do
    end if 
!
!   Use the following block to debug the construction of the Lewis structure
!
    if ( .false.) then
!
!  Sanity check - are all atomic orbitals accounted for?
!
      k = abs(norbs - noccupied - nvirtual)
      do i = 1, numat
        if (iz(i) /= 0) then
          k = k + 1
        end if
        if (ib(i) /= 0) then
          k = k + 1
        end if
      end do
      iz = 0
      do i = 1, Lewis_tot
        j = Lewis_elem(1,i)
        if (j > 0) iz(j) = iz(j) + 1
        j = Lewis_elem(2,i)
        if (j > 0) iz(j) = iz(j) + 1
      end do
      do i = 1, numat
        if (iz(i) - natorb(nat(i)) /= 0) then
          k = k + 1
        end if
      end do
      if (k /= 0) then
        write(*,*)" An error has been detected."
      end if
    end if
!
!
    if (lres) then      
      txtatm(:) = " "
      angles = 0.d0
      allres = " "
      if (done) then        
!
!   THE EARLIER CALL TO RESEQ MEANS THAT N1 MIGHT HAVE MOVED.
!   SO FIND N1 AGAIN.
!
        ioptl(:numat) = .false.
        call findn1 (n1, ioptl, io)
      end if
      l_protein = (n1 /= 0)
      ib(:numat) = -100000
      nfrag = 0
      ires = 0
      uni_res = 0
      odd_h = .true.
!
!  Break all intra-chain bonds, so that the residues can easily be
!  identified.
!
      call lyse !
      allr = " "
      do
        nfrag = nfrag + 1
        if (start_res(nfrag) /= -200) ires = start_res(nfrag)

! General bug-trap: if the number of fragments is unreasonably large,
! assume the system is unrecognizable and exit
        if (nfrag > 99) then
          call mopend("STRUCTURE UNRECOGNIZABLE")
          inquire(unit=iarc, opened=opend) 
          if (opend) close(iarc, status = "DELETE")          
          go to 1100
        endif
!
        call names (ioptl, ib, n1, ires, nfrag, io, uni_res, mres)
        if (moperr) return
!    
        nn1 = n1
        call findn1 (n1, ioptl, io)
        if (.not. l_protein) nfrag = 0
        if (n1 == 0) exit  
        if (n1 == nn1) ioptl(n1) = .true.   
        if (.not. lbreaks) then
          mbreaks = mbreaks + 1
          breaks(mbreaks) = n1
!
! Find the last atom that has been defined.  Use this for defining the break
!
          do i = 1, numat
            if (.not. ioptl(i)) exit
          end do
          if (i - 1 > 0) break_coords(:,mbreaks) = coord(:,i - 1)          
        end if             
      end do
!
!  Re-evaluate all residues
!
      j = 1
      allres(j) = txtatm(1)(18:20)
      do i = 2, natoms
        if (txtatm(i) == " ")exit
        if (nat(i) /= 1 .and. txtatm(i)(23:26) /= txtatm(i - 1)(23:26)) then
        j = j + 1
        allres(j) = txtatm(i)(18:20)
        end if
      end do
      ires = j      
      iopt(:natoms) = ib(:natoms)
!
!   LABEL THE ATOMS IN ANY NON-PROTEIN MOLECULES IN THE SYSTEM
!
      nfrag = nfrag + 1
      if (start_res(nfrag) == -200) then
        if (.not. l_protein) ires = 0       
      else
        ires = start_res(nfrag)         
      end if
      call ligand (ires, start_res, nfrag)
!
!  If ligands are present, set a break at the end of the protein, to separate protein from ligands
!
      if (.not. lbreaks) then
        mbreaks = mbreaks + 1
        breaks(mbreaks) = n1
!
! Find the last atom that has been defined.  Use this for defining the break
!
        do i = 1, numat
          if (.not. ioptl(i)) exit
        end do
        if (i - 1 > 0) break_coords(:,mbreaks) = coord(:,i - 1)          
      end if         
      nres = uni_res
      maxtxt = 26
!
!  Add chain letters
!
      i = index(keywrd, " RESI")
      if (i > 0) then
        j = index(keywrd(i + 1:), " ") + i
        j = index(keywrd(i + 1: j), "0")
      end if
      if (i > 0 .and. j > 0) then
        do i = 1, numat
          txtatm(i)(13:16) = txtatm1(i)(13:16)
        end do
      end if
      if (moperr) return
      mbreaks = 1
      line = txtatm(1)(23:26)
      do i = 1, numat
        txtatm(i) = txtatm(i)(:21)//chains(mbreaks)//txtatm(i)(23:)
        if (i == breaks(mbreaks)) then
          mbreaks = mbreaks + 1
        else  
          if (txtatm(i)(23:26) == "    ") cycle          
          if (txtatm(i)(:6) == "HETATM") then
            if (txtatm(i)(23:26) /= line(1:4)) then
              if(txtatm(i - 1)(:4) /= "ATOM") mbreaks = min(26, mbreaks + 1)
            end if
            txtatm(i) = txtatm(i)(:21)//chains(mbreaks)//txtatm(i)(23:)            
          end if
        end if
        line = txtatm(i)(23:26)
      end do 
!
!   Check for unknowns
!
      do i = 1, numat
        do j = 1, 20
          if (txtatm(i)(18:20) == tyres(j)) exit
        end do
        if (j == 21) then
          if (txtatm(i)(15:16) /= "  ") cycle
          j = 1
          do k = i + 1, numat
            if (txtatm(k)(14:26) == txtatm(i)(14:26)) then
              j = j + 1
              if (j < 10) then
                write(txtatm(k)(15:15),'(i1)')j
              else
                 write(txtatm(k)(15:16),'(i2)')j
              end if                
            end if
          end do
          if (j > 1) write(txtatm(i)(15:15),'(i1)')1
        end if
      end do
      if (n_new /= 0) then
!
!  Re-name residues to use the XENO name
!
        mbreaks = 1
        old_name(:n_new) = "---"
        l_names = .false.
        do i = 1, natoms
          ter = (i == breaks(mbreaks))
          if (ter) mbreaks = mbreaks + 1
          j = nint(reada(txtatm(i)(23:),1))
          if (txtatm(i)(:6) == "HETATM") cycle
          do k = 1, n_new
            if (new_chain(k) == chains(mbreaks) .and. new_res(k) == j) then
              if (old_name(k) == "---") old_name(k) = txtatm(i)(18:20)
              txtatm(i) = txtatm(i)(:17)//new_name(k)//txtatm(i)(21:)
              l_names(k) = .true.
              exit
            end if
          end do
        end do
        first = .true.
        do i = 1, n_new
          if (l_names(i)) then
            if (old_name(i) ==  new_name(i)) then
              if (first) then
                write(*,'(/,a)') "      Residue names that have been changed"
                write(*,'(a)') "      Residue No.  Calc'd name   XENO name"
                first = .false.
              end if
              write(*,'(i3,i10,2x,a1,7x,a3,9x,a3, a)')i, new_res(i), new_chain(i), old_name(i), new_name(i), &
                "  Name not changed!"
            else
              write(*,'(i3,i10,2x,a1,7x,a3,9x,a3)')i, new_res(i), new_chain(i), old_name(i), new_name(i)
            end if 
          end if
        end do
        !if (log) then
        !  write(ilog,'(/,a)') "      Residue names that have been changed"
        !  write(ilog,'(a)') "      Residue No.  Calc'd name   XENO name"
        !  do i = 1, n_new
        !    if (old_name(i) ==  new_name(i)) then
        !      write(ilog,'(i3,i10,2x,a1,7x,a3,9x,a3, a)')i, new_res(i), new_chain(i), old_name(i), new_name(i), &
        !      "  Name not changed!"
        !    else
        !      write(ilog,'(i3,i10,2x,a1,7x,a3,9x,a3)')i, new_res(i), new_chain(i), old_name(i), new_name(i)
        !    end if        
        !  end do
        !end if
      end if  
      l_use_old_labels = (index(keywrd," SITE=") /= 0 .and. index(keywrd, " ADD-H") == 0)
      l_use_old_labels = .true.
      call update_txtatm(l_use_old_labels, .true.)
      call write_sequence
      if (index(keywrd, " RAMA") /= 0) then
        if (index(keywrd, " ADD-H") == 0 .and. uni_res > 1) then
          write(*,"(/10x,a)")"        Ramachandran Angles"
          write(*,"(10x,a,/)")"    Residue    phi    psi  omega"
        end if
        do i = 1, uni_res
          if (Abs(angles(1,i)) + Abs(angles(3,i))> 1.d-20 .and. res_start(i) > 0) &
          write(*,"(14x,a,3f7.1)")txtatm(res_start(i))(18:20)//txtatm(res_start(i))(23:26), (angles(j,i),j = 1,3)
        end do
        write(*,*)" "
      end if
    end if
    if (index(keywrd, " PDBOUT") /= 0) then
!
!  Identify atoms where chain breaks occur
!
      if (index(keywrd, " RESEQ") + index(keywrd, " ADD-H") == 0) then
        if (allocated(txtatm1) .and. index(keywrd, " RESID") /= 0) then
          if (txtatm1(1) /= " ") call compare_sequence(n_new)
        end if
      end if
    end if
!
!  Edit keywords to remove text that would not be used in the next calculation. 
!
    do i = 1, 6
      line = refkey(i)
      call upcase(line, len_trim(line))
      j = index(line, " SITE=")
      if (j > 0) then
        k = index(line(j:), ") ") + j
        refkey(i) = refkey(i)(:j)//refkey(i)(k:)
      end if
      j = index(line, " RESEQ")
      if (j /= 0) refkey(i) = refkey(i)(:j - 1)//trim(refkey(i)(j + 6:))
    end do
    line = " "//trim(refkey(1))
    call upcase(line, len_trim(line))
    i = index(line, " RESEQ")
    if (i /= 0) refkey(1) = refkey(1)(:i - 1)//trim(refkey(1)(i + 5:))
    if (index(keywrd, " ADD-H") /= 0) return
!
!  MODIFY IONS SO THAT IT REFERS TO ALL ATOMS (REAL AND DUMMY)
!
    j = 0
    iz = ions
    ions = 0
    do i = 1, natoms
      if (labels(i) == 99) then
        ions(i) = 0       
      else
         j = j + 1
        ions(i) = iz(j)
      end if
    end do
    ibad = 0
    if (lres .and. .not. lreseq) then
!
!   CHECK ALL IONS TO SEE IF ANY RESIDUE IS A DI-ION
!
      outer_loop: do i = 1, numat
        if(ions(i) /= 1 .and. ions(i) /= -1) cycle! WARNING
      end do outer_loop
      if (ibad /= 0) then
        write (*,*)
      end if
      jbad = ibad
      ibad = 0
      ibad = ibad + jbad
    end if
!
!
    charges = (lsite .or. index(keywrd, "CHARGES") /= 0)
    ichrge = -noccupied*2
    do i = 1, numat
      ichrge = ichrge + nint(tore(nat(i)))
    end do
    line = " "
    if (noccupied /= 0 .and. (index(keywrd," LEWIS") > 0 .or. noccupied*2 /= nelecs)) then
      maxtxt = 0
      do i = 1, numat
        maxtxt = max(maxtxt, len_trim(txtatm(i)))
      end do
      if (maxtxt == 0) then
        j = 1
      else
        j = maxtxt/2 + 2
      end if
      !call update_txtatm(.true., .false.)
      if (prt_topo) then
        write (*, "(/,A,/)") "   TOPOGRAPHY OF SYSTEM"
        write (*,*) "  ATOM No. "//line(:j)//"  LABEL  "//line(:j)//"Atoms connected to this atom"
        if (j == 0) then
          do i = 1, numat
            write (*, "(I7,9X,A,9I7)") i, elemnt (nat(i)) // "  ", (ibonds(j, i), j=1, nbonds(i))
          end do
        else
          if (maxtxt > 2) then
            do i = 1, numat
              write (*, "(I7,9X,A,9I7)") i, elemnt (nat(i)) // " (" // txtatm(i) (:maxtxt) // ") ", &
                (ibonds(j, i), j=1, nbonds(i))
            end do
          else
            do i = 1, numat
              write (*, "(I7,9X,A,9I7)") i, elemnt (nat(i)), (ibonds(j, i), j=1, nbonds(i))
            end do
          end if        
        end if 
      end if
    end if
    if (noccupied /= 0 .and. index(keywrd," LEWIS") > 0) then
      write (*, "(/37x,A,/)") "   Lewis Structure"
      if (index(keywrd, " LARGE") /= 0)  &
        write (*,"(23x,a,/)") "  ATOMS IN OCCUPIED LOCALIZED MOLECULAR ORBITALS"
      l = 4
      allocate(Lewis_formatted(Lewis_tot/l + 5, l))
      Lewis_formatted = " "
      k = 0
      j = 0
      do i = 1, Lewis_tot
        if (Lewis_elem(1,i) > 0) j = j + 1
      end do
      m = j/l + 1    
      write(*,"(4('     LMO  Atom  Atom    '))")  
      ii = 0
      jj = 1  
      do i = 1, Lewis_tot
        if (Lewis_elem(1,i) > 0) then
          k = k + 1
          ii = ii + 1
          if (ii > m) then
            ii = 1
            jj = jj + 1
          end if
          if (Lewis_elem(2,i) > 0) then
            write (Lewis_formatted(ii,jj), "(I8,I6,I6)") k, (Lewis_elem(j,i),j=1,2)
          else
            write (Lewis_formatted(ii,jj), "(I8,I6)") k, Lewis_elem(1,i)
          end if
        end if
      end do
      do i = 1, Lewis_tot
        if (Lewis_formatted(i,1) == " ") exit
        write(*,"(10a)")(Lewis_formatted(i,j)//"    ", j = 1, l)," "
      end do
      if (index(keywrd, " LARGE") /= 0) then
        write (*,"(/23x,a,/)") "  ATOMS IN UNOCCUPIED LOCALIZED MOLECULAR ORBITALS"
        Lewis_formatted = " "
        j = 0
        do i = 1, Lewis_tot
          if (Lewis_elem(2,i) > 0) j = j + 1
        end do
        m = j/l + 1
        k = 0
        ii = 0
        jj = 1  
        do i = 1, Lewis_tot
          if (Lewis_elem(2,i) > 0) then
            k = k + 1
            ii = ii + 1
            if (ii > m) then
              ii = 1
              jj = jj + 1
            end if
            if (Lewis_elem(1,i) > 0) then
              write (Lewis_formatted(ii,jj), "(I8,I6,I6)") k, (Lewis_elem(j,i),j=1,2)
            else
              write (Lewis_formatted(ii,jj), "(I8,I6)") k, Lewis_elem(2,i)
            end if
          end if
        end do
        do i = 1, Lewis_tot
          if (Lewis_formatted(i,1) == " ") exit
          write(*,"(10a)")(Lewis_formatted(i,j)//"    ", j = 1, l)," "
        end do
      end if
      deallocate(Lewis_formatted)
      write (*,*)
      if (noccupied > 0) then
        write (*,"(a,/)") "          Type          Number of Lewis structural elements identified"
      !  if (log) write (ilog,"(/,a,/)") "          Type          Number of Lewis structural elements identified"
      end if
      if (numbon(1) /= 0) then
        write (*, "(A,I6)") "         SIGMA BONDS   ", numbon (1)
      !  if (log) write (ilog, "(A,I6)") "         SIGMA BONDS   ", numbon (1)
      end if
      if (numbon(2) /= 0) then
        write (*, "(A,I6)") "         LONE PAIRS    ", numbon (2)
      !  if (log) write (ilog, "(A,I6)") "         LONE PAIRS    ", numbon (2)
      end if
      if (numbon(3) /= 0) then
        write (*, "(A,I6)") "         PI BONDS      ", numbon (3)
      !  if (log) write (ilog, "(A,I6)") "         PI BONDS      ", numbon (3)
      end if
      if (index(keywrd," LEWIS") > 0 .or. (noccupied*2 /= nelecs .and. noccupied /= 0)) then
        write(*,"(/,a,i6)")" Number of filled levels from atoms and charge:", nelecs/2
        write(*,"(a,i6)")" Number of filled levels from Lewis structure: ", noccupied
      !  if (log) write (ilog,"(/,a,i6)")" Number of filled levels from atoms and charge:", nelecs/2
      !  if (log) write (ilog,"(a,i6)")" Number of filled levels from Lewis structure: ", noccupied
      end if
      l = 0
      m = 0
      num = char(Int(log10(numat + 1.0)) + ichar("1") + 1) 
      do i = 1, numat
        if (.not. main_group(nat(i))) then
!
!  Element is a transition metal.  Work out its formal oxidation state
!
          k = ions(i)
          do j = 1, Lewis_tot
            if (Lewis_elem(1,j) /= 0 .and. Lewis_elem(2,j) /= 0) then
              if (Lewis_elem(1,j) == i) k = k + 1
              if (Lewis_elem(2,j) == i) k = k + 1
            end if
          end do
          if (m == 0) write(*,*)
          m = 1
          write(*,"(10x,a,i"//num//",a,sp,i3)")" Formal oxidation state of atom",i,", a "//trim(elemnt(nat(i)))//", is", k
          if (k < 0) l = 1
          if (k > 3) l = 1
        end if
      end do
    end if    
!
! Check for sulfate and phosphate
!
      do i = 1, numat
        if (nat(i) == 16 .and. txtatm(i)(18:20) == "SO4") then
          k = 2
          ions(i) = 0
          do j = 1,4
            l = ibonds(j,i)
            if (ions(l) == -1) then
              ions(l) = 0
              k = k - 1
              if (k == 0) exit
            end if
          end do
        end if
        if (nat(i) == 15 .and. txtatm(i)(18:20) == "PO4") then
          k = 1
          ions(i) = 0
          do j = 1,4
            l = ibonds(j,i)
            if (ions(l) == -1) then
              ions(l) = 0
              k = k - 1
              if (k == 0) exit
            end if
          end do
        end if
      end do
    num_ions = 0
    do i = 1, numat
      j = (Min(6, Max(-6, ions(i))))
      num_ions(j) = num_ions(j) + 1
    end do
    i = 0
    do j = 1,6
      i = i + num_ions(j) + num_ions(-j)
    end do
    if (i > 0) then
      if (index(keywrd," LEWIS") > 0) then
        write (*,*)
        write (*,"(a,/)") "          Type           Number of charged sites identified"
    !    if (log) write(ilog,"(/,a,/)") "          Type           Number of charged sites identified"
        do i = 1,6
        if (num_ions(i)  > 0) then
          write (*, "(9x,a,2x,i6)") ion_names(i) , num_ions(i)
    !      if (log) write(ilog, "(9x,a,2x,i6)") ion_names(i) , num_ions(i)
        end if
        if (num_ions(-i) > 0) then
          write (*, "(9x,a,2x,i6)") ion_names(-i), num_ions(-i)
    !      if (log) write(ilog, "(9x,a,2x,i6)") ion_names(-i), num_ions(-i)
        end if
        end do
        i = 0
        do j = 1,6
          i = i + num_ions(j) * j
        end do
        write(*,"(SP/,a,i5)")" SUM OF POSITIVE CHARGES", i
    !    if (log) write(ilog,"(SP/,a,i5)")" SUM OF POSITIVE CHARGES", i
        i = 0
        do j = 1,6
          i = i + num_ions(-j) * j
        end do
        write(*,"(a,i5)")" SUM OF NEGATIVE CHARGES", -i
    !    if (log) write(ilog,"(a,i5)")" SUM OF NEGATIVE CHARGES", -i
      end if
      padding = " "
      do i = 1, numat
        if (ions(i) /= 0) then
          if (nat(i) == 15 .or. nat(i) == 16) then
            kk = 0
            do jj = 1, nbonds(i)
              ii = ibonds(jj,i)
              if (nat(ii) == 8) then
                if (ions(ii) == -1) then
!
!  Found a PO4 or SO4.  Neutralize the P-O or S-O Zwitterion.
!
                  ions(i) = ions(i) - 1
                  ions(ii) = ions(ii) + 1
                  kk = 1
                  exit
                end if
              end if
            end do
            if (kk == 1) cycle
          end if
        end if
      end do
      do i = 1, numat
        if (ions(i) /= 0) exit
      end do
      maxtxt_store = maxtxt
      if (maxtxt < 0) maxtxt = 14
      if (i <= numat) then
        if (maxtxt > 1) then        
           write(line,"(a)")"   Ion Atom No.           Label               Charge"
           write(*,'(/,a)') trim(line)
        !   if (log) write(ilog,"(/,a)") trim(line)
           l = max(1,(17 - maxtxt/2))
           residues = (index(keywrd, " RESID") /= 0)
        else
           write(*,"(/,a)")"     Ion Atom No.  Type    Charge"
        !   if (log) write(ilog,"(/,a)")"    Ion Atom No.  Type    Charge"
        end if
      else
        write(*,'(/18x,a)') "NO CHARGED ATOMS FOUND."        
      end if
      do j = 6, -6, -1
        if (j == 0) cycle
        k = 0
        do i = 1, numat
          if (ions(i) == j) then
            if (k == 0) write(*,*)
            k = k + 1
            line = " "
            jj = 0
            if (j == 1) then
              do ii = 1, numat
                if (ions(ii) < 0) then
                  sum = (coord(1,i) - coord(1,ii))**2 + (coord(2,i) - coord(2,ii))**2 + (coord(3,i) - coord(3,ii))**2 
                  if (sum < 99.d0) then
                    jj = jj + 1
                    if (jj > 100) exit
                    near_ions(jj) = ii
                    r_ions(jj) = sqrt(sum)
                    line = " Angstroms from anion"
                  end if                  
                end if
              end do
            else if (j == -1) then
              do ii = 1, numat
                if (ions(ii) > 0) then
                  sum = (coord(1,i) - coord(1,ii))**2 + (coord(2,i) - coord(2,ii))**2 + (coord(3,i) - coord(3,ii))**2 
                  if (sum < 99.d0) then
                    jj = jj + 1
                    if (jj > 100) exit
                    near_ions(jj) = ii
                    r_ions(jj) = sqrt(sum)
                    line = " Angstroms from cation"
                  end if                  
                end if
              end do
            end if
            if (jj > 1) then
!
!  Sort near ions into increasing distance
!
              jj = min(jj,100)
              do ii = 1, jj
                do kk = ii + 1, jj
                  if (r_ions(kk) < r_ions(ii)) then
                    sum = r_ions(kk)
                    r_ions(kk) = r_ions(ii)
                    r_ions(ii) = sum
                    m = near_ions(kk)
                    near_ions(kk) = near_ions(ii)
                    near_ions(ii) = m
                  end if
                end do 
              end do
              do ii = 2, jj
                if (r_ions(ii) > 5.d0) then
                  jj = ii - 1
                  exit
                end if
              end do
            end if
            if (maxtxt > 1) then
              if (jj > 0) then
                if (residues) then
                  txtatm_1 = txtatm(i)
                  txtatm_2 = txtatm(near_ions(1))
                else
                  txtatm_1 = txtatm1(i)
                  txtatm_2 = txtatm1(near_ions(1))
                end if
                num = "5"
                if (elemnt(nat(near_ions(1)))(1:1) /= " ") num = "5"
                  m = len_trim(line) + 1
                  if (txtatm_2(maxtxt:maxtxt) == ")") then
                    ii = maxtxt + 1
                    txtatm_2 = "("//trim(txtatm_2)
                    kk = maxtxt - 1
                  else
                    ii = maxtxt
                    kk = maxtxt
                  end if
                  write(tmp,"(i5, 2x, a2, i5,3x, a, SP,i5,S,a,f3.1,a,i"//num//",a)") &
                  & k, elemnt(nat(i)), i, padding(:l-5)//"("//txtatm_1(1:kk)//")"//padding(:l-15), ions(i), &
                  & "  (",r_ions(1),line(:m)//elemnt(nat(near_ions(1))),near_ions(1),  &
                  &", Label: "//txtatm_2(:ii)//")"
                  write(*,'(a)')trim(tmp)
!
            !      if (log) write(ilog,'(a)')trim(tmp)
                do ii = 2, jj
                  if (residues) then
                    txtatm_2 = txtatm(near_ions(ii))
                  else
                    txtatm_2 = txtatm1(near_ions(ii))
                  end if
                  num = "5"
                  if (elemnt(nat(near_ions(ii)))(1:1) /= " ") num = "5"
                    write(tmp,"(49x, a, f3.1,a,i"//num//",a)")  "   (",r_ions(ii), &
                      line(:m)//elemnt(nat(near_ions(ii))), near_ions(ii), &
                      &", Label: "//txtatm_2(:maxtxt)//")"
                    write(*,'(a)')trim(tmp)
            !        if (log) write(ilog,"(a)") trim(tmp)
                end do
              else
                 if (residues) then
                  txtatm_1 = txtatm(i)
                else
                  txtatm_1 = txtatm1(i)
                end if
                kk = len_trim(txtatm_1)
                if (txtatm_1(kk:kk) == ")") then
                  kk = maxtxt - 1
                  ii = l - 6
                  j2 = l - 1
                else
                  kk = maxtxt
                  ii = l - 5
                  j2 = l - 5
                end if
                write(line,"(i5, 2x, a2, i5,3x, a,SP,i5,S,a)") &
                  k, elemnt(nat(i)), i, padding(:ii)//"("//txtatm_1(1:kk)//")"//padding(:j2), ions(i)
                write(*,'(a)') trim(line)
            !    if (log) write(ilog,"(a)") trim(line)
              end if
            else
              write(line,"(i5,3x,i5,5x,a2,SP,i9,S)")k, i, elemnt(nat(i)),ions(i)
              write(*,'(a)') trim(line)
            !  if (log) write(ilog,"(a)")trim(line)
            end if            
          end if
        end do
      end do
      maxtxt = maxtxt_store
    end if
    if (index(keywrd, " LEWIS") /= 0 .and. index(keywrd," 0SCF") == 0) &
      call mopend ("RUN STOPPED BECAUSE KEYWORD LEWIS WAS USED.")
    if (noccupied /= 0 .and. (charges .or. .not. lreseq .and. (ichrge /= 0 .or. irefq /= 0))) then
      num = char(ichar("3") + max(int(log10(abs(ichrge) + 0.05)),0))
      write (*, "(SP/17x,A,I"//num//")") " COMPUTED CHARGE ON SYSTEM:", ichrge
      computed_charge = ichrge
    !  if (log) write (ilog, "(SP/17x,A,I"//num//")") " COMPUTED CHARGE ON SYSTEM:", ichrge
    end if
!
    atmass(1:numat) = ams(nat(1:numat))
    if (done .and. .not. lreseq) then
      call xyzint (coord, numat, na, nb, nc, 1.d0, geo)
    end if
    if (lreseq) then
      geo(:,:numat) = coord(:,:numat)
      na = 0
    end if   
    if (lreseq .or. lsite) then    
      archive_fn = archive_fn(:len_trim(archive_fn) - 3)//"arc"
      inquire(unit=iarc, opened=opend) 
      if (opend) close(iarc, status = 'keep', iostat=i)  
      open(unit=iarc, file=archive_fn, status='UNKNOWN', position='asis') 
      rewind iarc 
      do i = 1, 6
        line = refkey(i)
        call upcase(line, len_trim(line))
        j = index(line, " PDBOUT")
        if (j > 0) refkey(i) = refkey(i)(:j)//refkey(i)(j + 7:)
      end do
      !call geout (iarc)
    end if
    if (irefq /= ichrge .and. .not. lreseq .or. charges) then
!
!  THE CALCULATED CHARGE DOES NOT MATCH THAT DEFINED BY CHARGE=N.
!  THEREFORE, THE USER HAS MADE A MISTAKE.  WRITE OUT CHARGES
!  FOUND HERE.
!
!
      if (icharges /= 0) write (*,*)   
      line = " "
      if (index(keywrd, " 0SCF") == 0) line = "JOB STOPPED BECAUSE"
      i = len_trim(line)
      if (i > 0) i = i + 1
      if (charges) then
        if (index(keywrd," CHARGES") == 0) then
          call mopend (line(:i)//"CHARGES MODIFIED BY SITE COMMAND")
        else
          call mopend (line(:i)//"CHARGES KEYWORD USED")
        end if
        nelecs = nelecs - ichrge
      end if    
      if (lreseq) call mopend(line(:i)//"GEOMETRY RESEQUENCED")
      if (charges .or. lreseq) return
      if (index(keywrd, " LEWIS") /= 0) then
        if (index(keywrd, " 0SCF") == 0) then
          if (.not. moperr) call mopend (line(:i)//"KEYWORD LEWIS USED")
          return
        end if
      end if 
!
!
      if (Index (keywrd, " 0SCF") + Index (keywrd, " RESEQ") + &
      & Index (keywrd, " LEWIS") == 0) then    
        write(*,*)" "
        num = char(ichar("2") + max(int(log10(abs(irefq) + 0.05)),0))
        write (*, "(10x,A,SP,I"//num//",A)") "In the data-set supplied, the charge specified (", irefq,") is incorrect."  
        if (index(keywrd," GEO-OK") /= 0 .and. index(keywrd," CHARGE=") /= 0) then
          write (*, "(10x,A)") " KEYWORD 'GEO-OK' WAS PRESENT, SO THE CHARGE HAS BEEN RESET."
          write (*, "(10x,A)") " IF THE NEW CHARGE IS INCORRECT, EITHER MODIFY THE STRUCTURE OR "//&
          "USE KEYWORD 'SETPI' TO CORRECT THE LEWIS STRUCTURE." 
        end if     
        nelecs = nelecs + irefq - ichrge     
        nclose = nelecs/2
        nopen = nclose
        nalpha = 0
        nbeta = 0
        uhf = .false.     
        if (irefq /= ichrge) then
          if (mod(irefq - ichrge, 2) == 0) then
            write(*,"(/10x,a)")"If the Lewis structure used by MOZYME is incorrect, "// &
            &"use keywords such as CVB or SETPI to correct it"
          else
            write(*,"(/10x,a)")"The charge keyword, Lewis structure, or the chemical formula is faulty"
          end if            
          if (Index(keywrd," GEO-OK") == 0 .and. index(keywrd," CHARGE=") /= 0 ) then  
          call mopend("CHARGE SPECIFIED IS INCORRECT. CORRECT THE ERROR BEFORE CONTINUING")
          write(*,"(/10x,a)")"If that is done, then the correct charge will be used."
          return  
        else if (id > 0) then
          write(*,"(10x,a)")" Infinite systems must have a zero charge on the unit cell."
          call mopend("Unit cell has a charge. Correct fault and re-submit ")
          return
        else 
          call fix_charges(ichrge)  
        end if
      end if
      end if
    end if
    if (done) then
      if (prt_coords) write (*, "(//10X,A,//)") " GEOMETRY AFTER RE-SEQUENCING"
      call update_txtatm(.true., .true.) 
      !if (prt_coords) call geout (*)
      if (index(keywrd, "0SCF") + index(keywrd, " RESEQ") == 0 .or. &
         index(keywrd, " PDBOUT") == 0) then
        inquire(unit=iarc, opened=opend) 
        if ( .not. opend) open(unit=iarc, file=archive_fn, status='UNKNOWN', position='asis')
        rewind(iarc)
        !call geout (iarc)
      end if
      line = "GEOMETRY RESEQUENCED"
      call mopend(trim(line))
      go to 1100
    end if
    if (ibad /= 0 .and. .not. let) then
      call mopend("ERROR")
      go to 1100
    end if
    if (Index (keywrd, " NEWGEO") /= 0) then
      call newflg ()
    end if
    if (times) then
      call timer (" END OF GEOCHK")
    end if
    if (lreseq .or. lsite) then     
      if (lsite) then
        call mopend ("Keyword SITE used")
        write(*,'(//,a)')" Run stopped because SITE used"
      else
        call mopend ("Keyword RESEQ used")
        write(*,'(//,a)')" Run stopped because RESEQ used"
      end if      
      return
    end if
!
!  MODIFY IONS SO THAT IT REFERS TO REAL ATOMS ONLY
!
    j = 0
    do i = 1, natoms
      if (labels(i) /= 99) then
        j = j + 1
        iz(i) = j
      end if
    end do
!
!  Restore charges, if present
!
    if (maxtxt /= 26) then
      do i = 1, natoms
        if(atom_charge(i) /= " ") txtatm(i)(2:2) = atom_charge(i)
      end do
    end if
!
!  Restore nbonds and ibonds in case they are modified within this subroutine
!
    nbonds(:numat) = nnbonds
    ibonds(:,:numat) = iibonds(:,:numat)
1100 continue
    if (Allocated (nnbonds))      deallocate (nnbonds, iibonds)
    if (Allocated (iz))           deallocate (iz)
    if (Allocated (ib))           deallocate (ib)
    if (Allocated (mb))           deallocate (mb)
    if (Allocated (atom_charge))  deallocate (atom_charge)
    if (Allocated (ioptl))        deallocate (ioptl)
    return
end subroutine geochk
subroutine extvdw_for_MOZYME (radius, refvdw)
    use molkst_C, only: numat, keywrd, line
    use common_arrays_C, only: nat
    use mod_atomradii, only: is_metal
    use elemts_C, only : cap_elemnt
    use reada_I
    implicit none
    double precision, dimension (numat), intent (out) :: radius
    double precision, dimension (106), intent (in) :: refvdw
    integer :: i, j, k
    double precision, dimension (106) :: vdw
    character (len=80) :: txt_rad
!
!  Modify Van der Waal's radii of various atoms.  The format of
! the keyword is: VDW(:chemical symbol=n.nn:chemical symbol=n.nn ...)
! e.g. VDWM(:H=1.0:Cl=1.7)
!
    txt_rad = " "
    i = index(keywrd," METAL")
    if (i /= 0) then
      j = index(keywrd(i:),")") + i
      txt_rad = keywrd(i:j)
      do k = 1,83
        if (.not. is_metal(k)) is_metal(k) = (index(txt_rad,cap_elemnt(k)) /= 0)
      end do
    end if
    i = index(keywrd," VDWM(")
    if (i == 0) then
      line = " "
    else
      if (keywrd(i + 6:i + 6) /= ";" .and. keywrd(i + 6:i + 6) /= ":") keywrd(i + 6:) = ";"//keywrd(i + 6:)
      j = index(keywrd(i + 6:),")") + i + 6
      do k = k + 5, j
         if (keywrd(k:k) == ":") keywrd(k:k) = ";"
         if (keywrd(k:k) == ",") keywrd(k:k) = ";"
      end do
      line = keywrd(i + 5:j - 1)
    end if
    vdw(:106) = refvdw(:106)
    if (line /= " ") then
      do i = 1, 106
        j = 2
        if (cap_elemnt(i)(2:2) == " ") j = 1
        k = index(line, ";"//cap_elemnt(i)(:j)//"=")
        if (k > 0) vdw(i) = reada(line, k)
      end do
    end if
!
!  Verify that all radii that will be used are, in fact, set correctly
!
    do i = 1, numat
      if (nat(i) > 102) cycle
      if (vdw(nat(i)) > 900.d0) then
        write (line,*) "MISSING VAN DER WAALS RADIUS FOR " // cap_elemnt (nat(i))
        call mopend(trim(line))
        j = 2
        if (cap_elemnt(i) (2:2) == " ") j = 1
        write (*, "(2x,3a)") "To correct this, add keyword 'VDWM(:", &
       & cap_elemnt (nat(i)) (1:j), "=n.nn)'"
        return
      end if
    end do
!
! Flag which elements are to be treated as metals by giving them negative
! atomic radii
    do i = 1, 102
      if (is_metal(i)) vdw(i) = -1.0d0
    end do
    do i = 1, numat
      radius(i) = vdw(nat(i))
    end do
!
!  Set VDW radius for individual atoms to -1
!
    if (txt_rad /= " ") then
      do i = 1, len_trim(txt_rad)
        if (ichar(txt_rad(i:i)) - ichar("0") <= 9 .and. ichar(txt_rad(i:i)) - ichar("0") > 0) then
          j = nint(reada(txt_rad(i:), 1))
          radius(j) = -1.d0
          do j = i, len_trim(txt_rad)
            if (ichar(txt_rad(i:i)) - ichar("0") > 9 .or. &
            ichar(txt_rad(i:i)) - ichar("0") < 0) exit
            txt_rad(i:i) = " "
          end do
        end if
      end do
    end if
    return
end subroutine extvdw_for_MOZYME
subroutine fix_charges(ichrge)  
  use molkst_C, only: refkey, keywrd, line
  implicit none
  integer, intent (in) :: ichrge
  integer :: i, j
!
!  Modify charges in keyword line so that the system will run using MOZYME.
!  Steps:
!     If "old" CHARGE=n keyword exists in the reference keyword line, delete it.
!     Write new CHARGE keyword, if charge is non-zero
!     Rewind the data set file.
!     Write out the new data set.
!     Rewind it again, so it is ready for reading.
!     Rewind the output file, so the user does not see the faulty output.
!
      line = trim(refkey(1))
      call upcase(line,len_trim(line))     
      i = Index(line," CHARGE=")
      if (i /= 0) then
        j = Index(refkey(1)(i+2:)," ")
        refkey(1)(i + 1:) = refkey(1)(i + j + 2:)
      end if
!                 12345678901 
      i = Index(refkey(1),"            ")
      if (ichrge /= 0) then
        if (ichrge > 99) then
          write(refkey(1)(i:i+11),'(" CHARGE=",i3)')ichrge
        else if (ichrge > 9) then
          write(refkey(1)(i:i+11),'(" CHARGE=",i2)')ichrge
        else if (ichrge > 0) then
          write(refkey(1)(i:i+11),'(" CHARGE=",i1)')ichrge
        else if (ichrge > -10) then
          write(refkey(1)(i:i+11),'(" CHARGE=",i2)')ichrge
        else 
          write(refkey(1)(i:i+11),'(" CHARGE=",i3)')ichrge
        end if
      end if
       i = Index(keywrd," CHARGE=")
      if (i /= 0) then
        j = Index(keywrd(i+2:)," ")
        keywrd(i + 1:) = keywrd(i + j + 2:)
      end if
      i = Index(keywrd,"            ")
      if (ichrge /= 0) then
        if (ichrge > 99) then
          write(keywrd(i:i+11),'(" CHARGE=",i3)')ichrge
        else if (ichrge > 9) then
          write(keywrd(i:i+11),'(" CHARGE=",i2)')ichrge
        else if (ichrge > 0) then
          write(keywrd(i:i+11),'(" CHARGE=",i1)')ichrge
        else if (ichrge > -10) then
          write(keywrd(i:i+11),'(" CHARGE=",i2)')ichrge
        else 
          write(keywrd(i:i+11),'(" CHARGE=",i3)')ichrge
        end if
      end if
end subroutine fix_charges
subroutine add_sp_H(i1, i, i2)
!
! Given three atoms, i1, i, and i2, put a hydrogen atom at the point 
! of a triangle
!
  use molkst_C, only: natoms, maxatoms
  use common_arrays_C, only : geo, nat, txtatm
  implicit none
  integer :: i1, i, i2
  logical :: first = .true.
  double precision :: scale
  natoms = natoms + 1
  if (natoms > maxatoms) then
    if (first) then
      write(*,*)" Too many changes. Re-run using the data set generated by this job"
      first = .false.
    end if
    natoms = natoms - 1
    return
  end if
  geo(:,natoms) = 2.d0*geo(:,i1) - 2.d0*geo(:,i) + geo(:,i2)
  scale = 1.1d0/sqrt( (geo(1,i1) - geo(1,natoms))**2 + &
                      (geo(2,i1) - geo(2,natoms))**2 + &
                      (geo(3,i1) - geo(3,natoms))**2 )
  geo(:,natoms) = geo(:,i1) + scale*(geo(:,natoms) - geo(:,i1))
  nat(natoms) = 1
  txtatm(natoms) = " "
end subroutine add_sp_H
subroutine add_sp2_H(i1, i, i2)
!
! Given three atoms, i1, i, and i2, put a hydrogen atom at the point 
! of a triangle
!
  use molkst_C, only: natoms, maxatoms
  use common_arrays_C, only : geo, nat, txtatm
  implicit none
  integer :: i1, i, i2
  logical :: first = .true.
  double precision :: scale
  natoms = natoms + 1
  if (natoms > maxatoms) then
    if (first) then
      write(*,*)" Too many changes. Re-run using the data set generated by this job"
      first = .false.
    end if
    natoms = natoms - 1
    return
  end if
  geo(:,natoms) = 3.d0*geo(:,i) - geo(:,i1) - geo(:,i2)
  scale = 1.1d0/sqrt( (geo(1,i) - geo(1,natoms))**2 + &
                      (geo(2,i) - geo(2,natoms))**2 + &
                      (geo(3,i) - geo(3,natoms))**2 )
  geo(:,natoms) = geo(:,i) + scale*(geo(:,natoms) - geo(:,i))
  nat(natoms) = 1
  txtatm(natoms) = " "
end subroutine add_sp2_H
subroutine add_sp3_H(i1, i, i2, i3)
!
! Given four atoms, i1, i, and i2, put a hydrogen atom at the point 
! of a tetrahedron
!
  use molkst_C, only: natoms, maxatoms
  use common_arrays_C, only : geo, nat, txtatm
  implicit none
  integer :: i1, i, i2, i3
  logical :: first = .true.
  double precision :: scale
  natoms = natoms + 1
  if (natoms > maxatoms) then
    if (first) then
      write(*,*)" Too many changes. Re-run using the data set generated by this job"
      first = .false.
    end if
    natoms = natoms - 1
    return
  end if
  geo(:,natoms) = 4.d0*geo(:,i) - geo(:,i1) - geo(:,i2) - geo(:,i3)
  scale = 1.1d0/sqrt( (geo(1,i) - geo(1,natoms))**2 + &
                      (geo(2,i) - geo(2,natoms))**2 + &
                      (geo(3,i) - geo(3,natoms))**2 )
  geo(:,natoms) = geo(:,i) + scale*(geo(:,natoms) - geo(:,i))
  nat(natoms) = 1
  txtatm(natoms) = " "
end subroutine add_sp3_H

subroutine compare_sequence(n_new)
!
!  Compare calculated residue with residue name from data-set
!  and print any differences
!
  use common_arrays_C, only : txtatm, txtatm1, nat, chains, breaks
  use molkst_C, only: numat, line, maxtxt
  use MOZYME_C, only : tyr, tyres
  use chanel_C, only: log
  implicit none
  integer :: n_new
!
! Local
!  
  integer :: i_atom, i, j, i_delta = 0, new_res, old_res, previous = -200, &
    mbreaks, loop
  character :: old*3, new*3
  double precision, external :: reada
  logical :: first = .true.
  character :: num
  first = .true.
  mbreaks = 1
  i_atom = 0
  do loop = 1, numat
!
!  Find first non-hydrogen atom
!
    do 
      i_atom = i_atom + 1
      if (i_atom > numat) exit
       if (i_atom == breaks(mbreaks)) mbreaks = mbreaks + 1
      if (nat(i_atom) /= 1) exit
    end do
    if (i_atom > numat) exit
    new_res = nint(reada(txtatm(i_atom),23)) 
    new = txtatm (i_atom)(18:20)
    if (maxtxt == 14) then  
      old_res = nint(reada(txtatm1(i_atom),12))
      old = txtatm1(i_atom)(8:10)
    else
      old_res = nint(reada(txtatm1(i_atom),23))
      old = txtatm1(i_atom)(18:20)
    end if
    if (new_res - i_delta /= old_res) then
      i_delta = new_res - old_res
    end if
    if (new_res /= previous .and. n_new == 0) then
      if (old /= new) then
        if (first) then
          write(*,'(/16x,a)')"Residue names that have changed"
          write(*,'(7x,a,/)')"Original residue name   Calculated residue name"
          if (log) then
        !    write(ilog,'(/16x,a)')"Residue names that have changed"
        !    write(ilog,'(7x,a,/)')"Original residue name   Calculated residue name"
          end if
          first = .false.
          line="XENO=("
        end if
        do i = 1, 20
          if (old == tyres(i)) exit
        end do
        num = tyr(i)
        if (i == 21) num = " "
        write(*,'(i14,a, 2x,a3,3x, 11x,i5,a,2x,a3)')old_res, " "//chains(mbreaks), old, &
          new_res, " "//chains(mbreaks), new
        !if (log) write(ilog,'(i14,a, 2x,a3,3x, a1, 11x,i5,a,2x,a3)')old_res, " "//chains(mbreaks), old, &
        !  num, new_res, " "//chains(mbreaks), new
        if (i < 21) then
          j = len_trim(line) + 1
          if (j > 200) exit
          if (old_res > -1) then
            num = char(ichar("1") + int(log10(old_res + 0.05)))
            write(line(j:),'(a1,i'//num//',a1,a1,a1)')chains(mbreaks), old_res, "=", tyr(i),","
          else
            num = char(ichar("1") + int(log10(-old_res + 0.05)) + 1)
            write(line(j:),'(a1,i'//num//',a1,a1,a1)')chains(mbreaks), old_res, "=", tyr(i),","
          end if
        end if
      end if
      previous = new_res
    end if
  end do
  if (first .and. n_new == 0) then
    write(*,'(/,a)')"        Calculated and original residue sequences agree perfectly"
    !if (log) write(ilog,'(/,a)')"        Calculated and original residue sequences agree perfectly"
  else
    j = len_trim(line)
    if (j > 6 .and. n_new == 0) then
      line(j:j) = ")"
      write(line,'(a)')"(Use the XENO keyword to re-define unrecognized residues.)"
      write(*,'(/2x,a)')trim(line)
      !if (log) write(ilog,'(/a)')trim(line)
    end if
  end if
  return
end subroutine compare_sequence
subroutine update_txtatm(output, sort)
!
!  If "output" is true. and keyword RESIDUES is not present, then use input text in TXTATM1
!
!  if "sort" is true, then use the geometry in coorda to work out which atom is which.
!
  use common_arrays_C, only : nat, txtatm, txtatm1, coorda, nbonds, ibonds, coord, &
    breaks
  USE molkst_C, ONLY: numat, maxtxt, numat_old, keywrd
  implicit none
  logical :: output, sort, pdb, update_chain
!
! Local
!
  logical :: L_all
  integer :: i, j, k, l, H_Z
  integer, allocatable :: n_H(:), nn_H(:)
  if (maxtxt /= 26) return
!
!   Do the hydrogen atom labels need to be re-done?
!  "Yes," if atoms are added or removed or shuffled.
!
  if (index(keywrd," RESID") + index(keywrd," ADD-H") + &
    index(keywrd," SITE=") + index(keywrd," RESEQ") /= 0 .or. output) then
    H_Z = 1
  else
    H_Z = 0
  end if
  pdb = (index(keywrd, " PDBOUT") /= 0)
  L_all = (output .and. index(keywrd," RESID") == 0) 
  if (L_all .and. .not. sort .and. numat == numat_old) then
!
!  Do nothing!
!
    do i = 1, numat
      if (txtatm1(i) == " ") then
        txtatm1(i) = txtatm(i)
      else
        txtatm(i) = txtatm1(i)
      end if
    end do    
    goto 99
  end if
!
!  Add text to TXTATM to label hydrogen atoms and to add chain letter
!
  call set_up_dentate()  
  call check_cvs(.false.)
  call check_H(i)
  update_chain = (index(keywrd, " CHAINS=(") == 0)
  do i = 1, numat
    if (sort) then
      if (nat(i) > H_Z) then
        if (L_all .or. update_chain) then
          do j = 1, numat_old          
            if (abs(coord(1,i) - coorda(1,j)) > 0.01d0) cycle
            if (abs(coord(2,i) - coorda(2,j)) > 0.01d0) cycle
            if (abs(coord(3,i) - coorda(3,j)) > 0.01d0) cycle
            if (L_all) then
              if (txtatm1(j) /= " ") txtatm(i) = txtatm1(j)
            else
              if (txtatm1(j) /= " ") txtatm(i)(22:22) = txtatm1(j)(22:22)
            end if
            exit
          end do 
        end if
      else     
        if (nbonds(i) > 0) txtatm(i) = txtatm(ibonds(1,i))  
      end if
    else 
      if (L_all) then
        if (nat(i) > 1) then
          if (txtatm1(i) /= " ") txtatm(i) = txtatm1(i)
        else if (nbonds(i) > 0) then
          txtatm(i) = txtatm(ibonds(1,i)) 
        end if
      else
        if (update_chain .and. txtatm1(i) /= " ") txtatm(i)(22:22) = txtatm1(i)(22:22)
      end if
    end if
  end do
  99 continue
!
!  Number hydrogen atoms, if more than one on a heavy atom
!
  allocate(n_H(numat), nn_H(numat))
  nn_H = 0 
  do i = 1, numat
    if (nat(i) > H_Z) then
      j = 0
      do k = 1, nbonds(i)
        if (nat(ibonds(k,i)) == 1) j = j + 1
      end do
      n_H(i) = j
    end if
  end do
  do i = 1, numat
    if (nat(i) <= H_Z .and. nbonds(i) > 0) then
      k = ibonds(1,i)
      txtatm(i) = txtatm(k)(:12)//" H"//txtatm(k)(15:)
!
!  If a hydrogen atom is attached to an unlabeled carbon atom, make the hydrogen atom a terminal hydrogen
!
      if (nat(k) == 6 .and. txtatm(i)(15:15) == " ") txtatm(i)(15:15) = "T"
      if (n_H(k) == 1) then
        txtatm(i)(13:14) = " H"
      else
        nn_H(k) = nn_H(k) + 1
        txtatm(i)(13:14) = char(nn_H(k) + ichar("0"))//"H"
      end if         
    end if
    if (numat == numat_old) then
      if (txtatm1(i) == " ") txtatm1(i) = txtatm(i)  
    end if
  end do
!
!  Check for duplicate hydrogen atom labels
!
  do i = 1, numat
    if (nat(i) == 1) then
      l = 1
      do j = i + 1, numat
        if (nat(j) == 1) then
          if (txtatm(i)(13:) == txtatm(j)(13:)) then
            do k = 13, 16
              if (txtatm(j)(k:k) == " ") then
                l = l + 1
                txtatm(j)(k:k) = char(l + ichar("0"))
                if (numat == numat_old) then
                  if (txtatm1(j) == " ") txtatm1(j) = txtatm(j)  
                end if
                exit
              end if
            end do
          end if
        end if
        if (l == 9) exit
      end do
      if (l > 1) then
        do k = 13, 16
            if (txtatm(i)(k:k) == " ") then
              txtatm(i)(k:k) = "1"
              if (numat == numat_old) then
                if (txtatm1(i) == " ") txtatm1(i) = txtatm(i)  
              end if
              exit
            end if
        end do
      end if
    end if
  end do 
  
!
!  Add atom numbering, using PDB format, i.e., a TER group has its own number
!  (TER groups are represented by BREAKS)
!
  j = 1
  do i = 1, numat
    write(txtatm(i),'(a6,i5,a15)')txtatm(i)(:6),i + j - 1,txtatm(i)(12:)  
    if (PDB .and. i == breaks(j)) j = j + 1     
  end do 
  return
  end subroutine update_txtatm
  subroutine rectify_sequence()
    use common_arrays_C, only : geo, coord, txtatm, txtatm1, nat, labels
    use molkst_C, only: numat, maxtxt
    implicit none
!
!  Search atom sequence for the following structure:
!  A break in a chain: same chain letter, but the residue sequence increases suddenly.
!  After the increase, atoms of a residue in the same chain that are out-of-sequence are found.
!
! If this situation occurs, move the out-of-sequence atoms to their correct place.
!
! This should only be done if RESEQ or similar keyword is present.
!
    integer :: i, i_lower, res2, res1, res3, res, ii, jj, num_min, num_max
    character :: chain1*1, chain2*1, chain*1, used(numat)*1
    double precision, external :: reada
    if (maxtxt /= 26) return
    txtatm1(:numat) = txtatm(:numat)
    do i = 1, numat
      used(i) = txtatm1(i)(22:22)
    end do
    used(1) = "a"
    labels(:numat) = nat(:numat)
    geo(:,:numat) = coord(:,:numat)
    res1 = nint(reada(txtatm(1), 23))
    chain1 = "a"
!
!  ii keeps track of the new numbering sequence
!
    ii = 1
    do i_lower = 2, numat
      chain2 = used(i_lower)
      if (chain2 == "a") cycle           ! Atom already used
!
!   res1: Residue number of previous residue
!   res2: Residue number of current residue
!
      res2 = nint(reada(txtatm1(i_lower), 23))
!
! Look for a fault at the end of the first bit of chain
! This will occur in atoms in the domain res1 to res1 + 1
!
      do res3 = res1, res1 + 1
        num_min = max(i_lower - 50, 1)
        num_max = min(i_lower + 50, numat)
        do 
          do jj = num_min, num_max
            res = nint(reada(txtatm1(jj), 23))
            chain = used(jj)
            if (res == res3 .and. chain == chain1 .and. used(jj) /= "a") exit
          end do
          if (jj > num_max) exit
          ii = ii + 1
          coord(:,ii) = geo(:,jj)
          txtatm(ii) = txtatm1(jj)
          nat(ii) = labels(jj)
          used(jj) = "a"   
        end do 
      end do
!
! Look for a fault at the start of the next bit of chain
! This will occur in atoms in the domain res2 - 1 to res2
!
      do res3 = res2 - 1, res2 
        num_min = max(i_lower - 50, 1)
        num_max = min(i_lower + 50, numat)
        do 
          do jj = num_min, num_max
            res = nint(reada(txtatm1(jj), 23))
            chain = used(jj)
            if (res == res3 .and. chain == chain1 .and. used(jj) /= "a") exit
          end do
          if (jj > num_max) exit
          ii = ii + 1
          coord(:,ii) = geo(:,jj)
          txtatm(ii) = txtatm1(jj)
          nat(ii) = labels(jj)
          used(jj) = "a"   
        end do 
      end do   
      chain1 = chain2
      res1 = res2
    end do  
  end subroutine rectify_sequence
  subroutine write_sequence
    use common_arrays_C, only : txtatm, nat
!
    use MOZYME_C, only : ions, allres, tyr, allr, tyres, maxres
!
!    
    use molkst_C, only: numat, line, numcal
    implicit none
    integer :: i, nfrag, jj, ii, ires, kl, ku, irold, l, j, k, charge, icalcn = -50, &
      ifrag
    logical, save :: prt = .true.
    character :: chain, ch*2
    double precision, external :: reada
    if ( .not. prt) return
    if (icalcn == numcal) then
      return
    else
      icalcn = numcal
      prt = .false.
    end if
!
!  Write out all residue names, use information in txtatm to define sequence
!
      ii = 1
      ifrag = 0
      do nfrag = 1, 100
        if (ii > numat) exit
        if (txtatm(ii) == " ") exit
        chain = txtatm(ii)(22:22)        
        ires = nint(reada(txtatm(ii), 23))        
        charge = ions(ii)
        irold = ires
        j = 0
        do ii = ii + 1, numat
          if (nat(ii) /= 1 .and. txtatm(ii)(22:22) /= chain) exit !  chain letter has changed.
          jj = nint(reada(txtatm(ii),23))
          if (nat(ii - 1) /= 1) j = ii - 1
          charge = charge + ions(ii)
          if (ires == jj) cycle              !  Residue is same
          charge = charge - ions(ii)
          if (txtatm(ii)(14:14) == "H") cycle!  Ignore hydrogen atoms
          if (ires + 1 /= jj) exit           !  Residue is not contiguous
!
!  Take residue name from the previous atom (the atom label has just changed)
!
          allres(ires) = txtatm(ii - 1)(18:20)
          if (charge == 1) then
            allres(ires)(4:4) = "+"
          else if (charge == -1) then
            allres(ires)(4:4) = "-"
          end if
          charge = 0
          ires = ires + 1
          if (ires > maxres) then
            write(line,'(a,i5)')" Maximum residue number allowed:", maxres
            call mopend(trim(line))
            return
          end if          
        end do
        if (j /= 0) then
          if (txtatm(j)(18:20) /= "   ") allres(ires) = txtatm(j)(18:20) 
        end if
        if (charge == 1) then
          allres(ires)(4:4) = "+"
        else if (charge == -1) then
          allres(ires)(4:4) = "-"
        end if
        do j = irold, ires
          do k = 1, 20
            if (allres(j) == tyres(k)) exit
          end do
          if (k < 21) exit          
        end do
        if (j > ires .or. (j == ires .and. j == 0)) cycle
        ifrag = ifrag + 1
        if (ifrag == 1) then
          write (*, "(/,16X,A,/)") "RESIDUE SEQUENCE IN PROTEIN Chain: "//chain
        else
          write (*, "(/,16X,A,I2, a/)") "RESIDUE SEQUENCE IN PROTEIN FRAGMENT:", ifrag, " Chain: "//chain
        end if
        if (irold < 0) then
          jj = 100
          i = mod(irold + jj, 10)
          if (i == 0) i = 10
          kl = irold
          ku = min(ires, kl - i + 10)         
          line = " "
          l = 6*(i - 1) + 1
          j = ((kl + jj)/10)*10 - jj
          if (i == 10) j = j - 10
          write(*,"(8x,10i6)")((k - 10), k = 1, 10)
          write(*,"(5x,i4,2x,a,10(a4,'  '))")j, line(:l), allres(kl:ku)      
          do
            kl = ku + 1
            if (kl > min(0, ires)) exit
            j = j + 10
            ku = min(ires, ku + 10)
            write(*,"(5x,i4,3x,10(a4,'  '))")j, allres(kl:ku)        
         end do
         write(*,'(a)')" "
        end if
        kl = max(irold, 1)
        i = mod(kl, 10)
        if (i == 0) i = 10        
        ku = min(ires, kl - i + 10)         
        line = " "
        l = 6*(i - 1) + 1
        j = (kl/10)*10
        if (i == 10) j = j - 10
        write(*,"(8x,10i6)")(k, k = 1, 10)
        write(*,"(5x,i4,2x,a,10(a4,'  '))")j, line(:l), allres(kl:ku)
        do
          kl = ku + 1
          if (kl > ires) exit
          j = j + 10
          ku = min(ires, ku + 10)
          write(*,"(5x,i4,3x,10(a4,'  '))")j, allres(kl:ku)        
        end do
      end do
      ii = 1
      ifrag = 0      
      do nfrag = 1, 100
        j = 0
        if (ii > numat) exit
        if (txtatm(ii) == " ") exit
        chain = txtatm(ii)(22:22)        
        ires = nint(reada(txtatm(ii), 23))
        irold = ires
        do ii = ii + 1, numat
          if (nat(ii) /= 1 .and. txtatm(ii)(22:22) /= chain) exit !  chain letter has changed.
          jj = nint(reada(txtatm(ii),23))
          if (nat(ii - 1) /= 1) j = ii - 1
          if (ires == jj) cycle              !  Residue is same
          if (txtatm(ii)(14:14) == "H") cycle!  Ignore hydrogen atoms
          if (ires + 1 /= jj) exit           !  Residue is not contiguous
          allres(ires) = txtatm(ii - 1)(18:20)  
          ires = ires + 1
        end do
        if (j /= 0) then
          if (txtatm(j)(18:20) /= "   ") allres(ires) = txtatm(j)(18:20) 
        end if 
        do j = irold, ires
          do k = 1, 20
            if (allres(j) == tyres(k)) exit
          end do
          if (k < 21) exit          
        end do
        if (j > ires .or. (j == ires .and. j == 0)) cycle
        do i = irold, ires
          do k = 1, 20
            if (tyres(k) == allres(i)) exit
          end do
          if (k > 0 .and. k < 21) then
            allr(i) = tyr(k)
          else if (allres(i) /= " ") then
            allr(i) = "X"
          else
            allr(i) = " "
          end if
        end do
        ifrag = ifrag + 1
        if (ifrag == 1) then
          write (*, "(/,16X,A,/)") "RESIDUE SEQUENCE IN PROTEIN Chain: "//chain
        else
          write (*, "(/,16X,A,I2, a/)") "RESIDUE SEQUENCE IN PROTEIN FRAGMENT:", ifrag, " Chain: "//chain
        end if   
        jj = 100
        i = mod(irold + jj, 10)
        if (irold < 0) i = i - 10
        if (i == 0) i = 10
        kl = irold
        if (i < 0) then
          ku = min(ires, kl - i + 40) 
        else
          ku = min(ires, kl - i + 50) 
        end if
        line = " "
        j = ((kl + jj)/10)*10 - jj
        if (i == 10) j = j - 10
        if (i == 1) then
          ch = "10"
        else
          ch(1:1) = " "
          if (i < 0) then           
             ch(2:2) = char(9 + i + ichar("0"))
             if (ch(2:2) == "0") then
               ch(1:1) = "1"
               i = 12
             else
               i = 2 - i
             end if
          else
             ch(2:2) = char(11 - i + ichar("0"))
          end if         
        end if
        if (i /= 0) write (*, "(5x,i4,a,"//ch//"a,x,10(10a1,x))")j, line(:i + 1), (allr(i), i=kl, ku)
        do 
          kl = ku + 1
          if (kl > ires) exit
          j = j + 50
          ku = min(ires, kl + 49)
          write (*, "(5X,I4,2X,5(10A1,1X))") j, (allr(i), i=kl, ku)
        end do
      end do
  end subroutine write_sequence
  
  subroutine find_salt_bridges(in_cat, in_ani, n_cat, n_ani)
!
!  Find Salt Bridges:  Identify all ionizable sites that can be used for making salt bridges
!
    use common_arrays_C, only : nat, txtatm, nbonds, ibonds
    use molkst_C, only : numat, line, keywrd, moperr
    implicit none
    integer, intent(in) :: n_cat, n_ani, in_cat(n_cat), in_ani(n_ani)    
    integer :: i, j, k, l, m, n, n_cations, n_anions, n_C, n_H, n_O, C, pairs(2,200), n_pairs, &
      salt_bridges(2,200), n_salt, i_length, ii, jj, nh_cat, nh_ani
    double precision :: r, cutoff, Rab(200), R_min, R_sorted(200)
    character :: bits*10
    double precision, external :: distance, reada
    integer, allocatable :: cations(:), anions(:)
    allocate (cations(numat), anions(numat))
    if (n_cat == 0) then
      i = index(keywrd,"SALT=")
      if (i /= 0) then
        cutoff = reada(keywrd, i + 4)
      else
        cutoff = 4.d0
      end if 
      nh_cat = 0
      nh_ani = 0
    else 
      cutoff = 8.d0
      nh_cat = 1
      nh_ani = -1
    end if    
!
!  Locate all Arg guanidine "C" atoms and all -C(R)H-NH2 "N" atoms
!
    n_cations = 0
    do i = 1, numat
      if (nat(i) == 6) then
        if (nbonds(i) == 3) then
          if (nat(ibonds(1,i)) == 7 .and. nat(ibonds(2,i)) == 7 .and. nat(ibonds(3,i)) == 7) then
!
!  Make sure it has the correct charge
!
            m = 0
            do j = 1, 3
              n = ibonds(j,i)
              do k = 1, nbonds(n)
                if (nat(ibonds(k,n)) == 1) m = m + 1
              end do
            end do
            if (m == 4 + nh_cat) then
              if (n_cat > 0) then
                do n = 1, n_cat
                  if (txtatm(in_cat(n))(18:26) == txtatm(i)(18:26)) then                   
                    exit
                  end if
                end do
              else
                n = 0              
              end if
              if (n <= n_cat .and. txtatm(i)(18:20) /= "UNK" &
                .and. txtatm(i)(22:22) >= "A" .and. txtatm(i)(22:22) <= "Z") then
                n_cations = n_cations + 1
                cations(n_cations) = i
              end if
            end if
          end if
        end if
      end if
    end do
    do i = 1, numat
      if (nat(i) == 7) then
        if (nbonds(i) == 3 + nh_cat) then
!
!  Test for -NH2, e.g. N-terminus and lysine 
!
          n_C = 0
          n_H = 0
          do j = 1, 3 + nh_cat
            if (nat(ibonds(j,i)) == 6  .and. txtatm(ibonds(j,i))(18:20) /= "UNK" &
                .and. txtatm(ibonds(j,i))(22:22) >= "A" .and. txtatm(ibonds(j,i))(22:22) <= "Z") then
              n_C = n_C + 1
              C = ibonds(j,i)
            end if
            if (nat(ibonds(j,i)) == 1) n_H = n_H + 1
          end do
!
!  Make sure it has the correct charge
!
          if (n_C == 1 .and. n_H == 2 + nh_cat) then
            do j = 1, nbonds(C)
              if (nat(ibonds(j,C)) == 8) exit
              if (nat(ibonds(j,C)) == 7 .and. ibonds(j,C) /= i) exit
            end do
            if (j > nbonds(C)) then
              if (n_cat > 0) then
                do n = 1, n_cat
                  if (txtatm(in_cat(n))(18:26) == txtatm(i)(18:26)) then                   
                    exit
                  end if
                end do
              else
                n = 0              
              end if
              if (n <= n_cat .and. txtatm(C)(18:20) /= "UNK" &
                .and. txtatm(C)(22:22) >= "A" .and. txtatm(C)(22:22) <= "Z") then
                n_cations = n_cations + 1
                cations(n_cations) = C
              end if
            end if
          end if
        else if (nbonds(i) == 2 + nh_cat .and. txtatm(i)(18:20) == "HIS") then
!
!  Check for Histidine
!
          do j = 1, nbonds(i)
            k = ibonds(j,i)
            if (nbonds(k) < 3) cycle
            do l = 1, nbonds(k)
              if (nat(ibonds(l,k)) == 7 .and. ibonds(l,k) /= i) exit
            end do
            if (l <= nbonds(k)) exit
          end do 
          if (j > nbonds(i)) cycle
          if (n_cat > 0) then
            do n = 1, n_cat
              if (txtatm(in_cat(n))(18:26) == txtatm(i)(18:26)) then                   
                exit
              end if
            end do
          else
            n = 0              
          end if
!
!  Avoid a "double count"
!
          do k = 1, n_cations
            if (txtatm(cations(k))(18:26) == txtatm(ibonds(j,i))(18:26)) exit
          end do
          if (n <= n_cat .and. k > n_cations) then
            n_cations = n_cations + 1
            cations(n_cations) = ibonds(j,i)
          end if
        end if
      end if
    end do     
!
!  Locate all -COOH "C" atoms
!
    n_anions = 0
    do i = 1, numat
      if (nat(i) == 6) then
      if (nbonds(i) == 3) then
        n_C = 0
        n_O = 0
        do j = 1,3
          if (nat(ibonds(j,i)) == 6) n_C = n_C + 1
          if (nat(ibonds(j,i)) == 8) n_O = n_O + 1
        end do
        if (n_C == 1 .and. n_O == 2) then
          m = 0
          ii = 0
          do j = 1, 3
            if (nat(ibonds(j,i)) == 8) then
              l = ibonds(j,i)
              ii = ii + nbonds(l)
              do k = 1, nbonds(l)
                if (nat(ibonds(k,l)) == 1) m = m + 1
              end do
            end if
          end do
          if (m == 1 + nh_ani .and. ii < 4) then
            if (n_ani > 0) then
              do n = 1, n_ani
                if (txtatm(in_ani(n))(18:26) == txtatm(i)(18:26)) then                   
                  exit
                end if
              end do
            else
              n = 0              
            end if
            if (n <= n_ani .and. txtatm(i)(18:20) /= "UNK" &
                .and. txtatm(i)(22:22) >= "A" .and. txtatm(i)(22:22) <= "Z") then
                n_anions = n_anions + 1
                anions(n_anions) = i
              end if
            end if
          end if
        end if
      end if
    end do
!
!  Find interatomic distances and weight them by type 
!
    n_pairs = 0
    do i = 1, n_cations
      do j = 1, n_anions
        r = distance(cations(i), anions(j))
        if (r < cutoff + 2.d0) then
          m = cations(i)
          n = anions(j)
          R_min = 100.d0
          do k = 1, nbonds(m)
            if (nat(ibonds(k,m)) /= 7) cycle
            do l = 1, nbonds(n)
              if (nat(ibonds(l,n)) /= 8) cycle
              r = distance(ibonds(k,m), ibonds(l,n))
              if (r < R_min) then
                R_min = r
                ii = ibonds(k,m)
                jj = ibonds(l,n)
              end if
            end do
          end do
          if (R_min < cutoff) then
            n_pairs = n_pairs + 1
            pairs(1,n_pairs) = ii
            pairs(2,n_pairs) = jj
            Rab(n_pairs) = R_min
          end if
        end if
      end do
    end do
!
!  Sort interatomic distances
!
    n_salt = 0
    do
      if (n_pairs == 0) exit
      r = 1.d10
      do i = 1, n_pairs
        if (r > Rab(i)) then
          r = Rab(i)
          j = i
        end if
      end do
      n_salt = n_salt + 1
      salt_bridges(:, n_salt) = pairs(:, j)
      R_sorted(n_salt) = r
      n_pairs = n_pairs - 1
      do k = j, n_pairs
        Rab(k) = Rab(k + 1)
        pairs(:, k) = pairs(:, k + 1)
      end do
    end do
    if (n_cat == 0) then
!
!  Eliminate duplicates
!
      do i = 1, n_salt
        if (salt_bridges(1,i) == 0) cycle
        do j = n_salt, i + 1, -1
          if (salt_bridges(1,j) == 0) cycle
          m = salt_bridges(1,i)
          n = salt_bridges(2,i)
          ii = salt_bridges(1,j)
          jj = salt_bridges(2,j)
          if (salt_bridges(1,j) == salt_bridges(1,i) .or. salt_bridges(2,j) == salt_bridges(2,i)) then
            salt_bridges(:,j) = 0
          else if (txtatm(m)(18:) == txtatm(jj)(18:)) then
            salt_bridges(:,j) = 0
          else
            do ii = 1,2
              m = salt_bridges(ii,i)
              n = salt_bridges(ii,j)
              do k = 1, nbonds(m)
                if (nat(ibonds(k,m)) /= 6) cycle
                do l = 1, nbonds(n)
                  if (nat(ibonds(l,n)) /= 6) cycle
                  if (ibonds(k,m) == ibonds(l,n))  salt_bridges(:,j) = 0
                end do
              end do
              if (salt_bridges(1,j) == 0) exit
            end do        
          end if        
        end do
      end do   
    end if
!
!  Build new SITE keyword
!
    line = " "
    if (n_salt > 0) then
      write(*,'(//19x,a)')"Salt Bridges Found"
      write(*,'(/4x,a,10x,a,28x,a,16x,a,/)')" No.","Cationic site","Anionic site","Dist. (Angstroms)"
     if (n_cat == 0) call update_txtatm(.true., .true.)
    else
       if (n_cat == 0) then
         call mopend("SALT option used in SITE command, but no Salt Bridges found")
         moperr = .false.
         return
       else
         write(*,'(//19x,a)')"No Salt Bridges Found"
       end if
    end if
    j = 0
    do i = 1, n_salt
      if (salt_bridges(1,i) == 0) cycle
      k = salt_bridges(1,i)
      do m = 26, 23, -1
        if (txtatm(k)(m:m) == " ") exit
      end do
      m = max(23,m + 1)
      l = salt_bridges(2,i)
      do n = 26, 23, -1
        if (txtatm(l)(n:n) == " ") exit
      end do
      n = max(23,n + 1)   
      i_length = len_trim(line)
      if (i_length > 900) exit
      write(line(i_length + 1:),'(2a)')","//txtatm(k)(22:22)//txtatm(k)(m:26)//"(+),", &
                                        txtatm(l)(22:22)//txtatm(l)(n:26)//"(-)"
      j = j + 1
      write(*,'(i7,2x,a,a9,4x,a,a9, f10.2)') &
      j, "("//txtatm(k)//")", txtatm(k)(22:22)//txtatm(k)(m:26)//"(+)", &
      "("//txtatm(l)//")", txtatm(l)(22:22)//txtatm(l)(n:26)//"(-)", R_sorted(i)
      if (txtatm(k)(22:22) < "A" .or. txtatm(k)(22:22) > "Z") then
        write(line,'(a)')"Chain letter for the cation is not in the range 'A' to 'Z'"      
        call mopend(trim(line))
        write(*,'(10x,a)')"Chain letter = """//txtatm(k)(22:22)//""""
        return
      end if
      if (txtatm(l)(22:22) < "A" .or. txtatm(l)(22:22) > "Z") then
        write(line,'(a)')"Chain letter for the anion is not in the range 'A' to 'Z'" 
        call mopend(trim(line))
        write(*,'(10x,a)')"Chain letter = """//txtatm(l)(22:22)//""""
        return
      end if
      
    end do
    if (n_cat > 0) return
    line = line(2:)
!
!  Delete the word "SALT"
!
    k = index(keywrd, "(SALT") + index(keywrd, ",SALT") + 1
    i = index(keywrd,"(SALT=") + index(keywrd, ",SALT=")
    if (i > 0) then
      do j = i, len_trim(keywrd)
        if (keywrd(j:j) == "," .or. keywrd(j:j) == ")") exit
      end do
      keywrd(i + 1:) = keywrd(j:)
      if (keywrd(i:i) == ",") keywrd(i:) = trim(keywrd(i + 1:))
    else
      i = index(keywrd, "(SALT") + index(keywrd, ",SALT") + 1
      keywrd(i:) = keywrd(i + 4:)
      if (keywrd(i:i) == ",") keywrd(i:) = trim(keywrd(i + 1:))
    end if  
!
!   Check for duplicate sites in keywrd
!
    m = index(keywrd," SITE")
    m = index(keywrd(m:), "(") + m
    j = index(keywrd(m:),") ") + m
    if (keywrd(m:m) == ",") m = m + 1
    do 
      k = index(keywrd(m:j), ",")
      if (k == 0) k = index(keywrd(m:j + 1), ") ")
      if (k == 0) exit
      k = k + m - 2
      bits = keywrd(m:k)
      if (k - m < 3) exit
      l = index(line, trim(bits))
      if (l > 0) then
!
!  Remove the comma
!
        if (line(l - 1:l - 1) == ",") then
          line = line(:l - 2)//line(l + k - m + 1:)
        else
          line = line(:l - 1)//line(l + k - m:)
        end if 
      end if
      m = k + 2
      if (keywrd(m:m) == ",") m = m + 1
    end do
!
!  Insert the new keyword
!   
    i = index(keywrd," SITE")
    if (keywrd(i + 7:i + 7) /= ")") then
      keywrd(i + 7:) = trim(line)//","//trim(keywrd(i + 8:))
    else
      keywrd(i + 7:) = trim(line)//trim(keywrd(i + 7:))
    end if
    return
  end subroutine find_salt_bridges



