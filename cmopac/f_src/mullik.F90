      subroutine mullik() 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!
      use molkst_C, only : numat, nelecs, nclose, nopen, fract,  &
        keywrd, norbs, id, verson, method_pm6, uhf, nalpha, nbeta, &
        numcal, escf, line
!
      use symmetry_C, only : jndex, namo
!
      use maps_C, only : rxn_coord
!
      use common_arrays_C, only : nfirst, nlast, nat, coord, &
      & c, h, pb, tvec, eigs, q, eigb, p, ifact, cb
!
      use parameters_C, only : zs, zp, zd, betas, betap, tore
!
      use chanel_C, only : igpt, gpt_fn, iw
!
      use mult_I
      Use mod_vars_cuda, only: lgpu
#if GPU     
      Use mod_vars_cuda, only: prec
#endif               
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, if, il, im1, k, ii, j, jf, jl, ij, icalcn=0, mo_l, mo_u
      character, dimension(:), allocatable :: namo_tmp*4  
      character :: paras*1
      real(double), dimension(norbs) :: eig 
      real(double), dimension(:), allocatable :: store, vecs, store_h
      real(double), dimension(:, :), allocatable :: work
      real(double) :: bi, bj, sum, summ, q2(numat)
      logical :: graph, graph_formatted, namo_ok
! GBR_new_addition
      integer :: nlower
!            
      save :: graph, graph_formatted, icalcn, mo_l, mo_u
!      
      
!-----------------------------------------------
!*********************************************************************
!
!  MULLIK DOES A MULLIKEN POPULATION ANALYSIS
! ON INPUT     C      =  SQUARE ARRAY OF EIGENVECTORS.
!              H      =  PACKED ARRAY OF ONE-ELECTRON MATRIX
!              VECS   =  WORKSTORE OF SIZE AT LEAST NORBS*NORBS
!              STORE  =  WORKSTORE OF SIZE AT LEAST (NORBS*(NORBS+1))/2
!
!*********************************************************************
 
     
     allocate(store((norbs*(norbs + 1))/2), store_h((norbs*(norbs + 1))/2), &
       vecs(norbs**2), stat = i)
     if (i /= 0) then
        write(iw,*)" Unable to allocate memory for eigenvector matrices in MULLIK"
        call mopend("Unable to allocate memory for eigenvector matrices in MULLIK")
        return
      end if
!*********************************************************************
!
!  FIRST, RE-CALCULATE THE OVERLAP MATRIX
!
!*********************************************************************
      if (icalcn /= numcal) then
        icalcn = numcal
        graph = index(keywrd,'GRAPH') /= 0 
        graph_formatted = index(keywrd,'GRAPHF') /= 0 
        if (allocated(ifact)) deallocate(ifact)
        allocate (ifact(norbs + 1))
        do i = 1, norbs 
          ifact(i) = (i*(i - 1))/2 
        end do 
        ifact(norbs+1) = (norbs*(norbs + 1))/2 
        mo_l = 1
        mo_u = norbs
      end if      
      call chrge (p, q2) 
      do i = 1, numat 
        q(i) = tore(nat(i) ) - q2(i)
      end do      
      store_h = h
      do i = 1, numat 
        if = nfirst(i) 
        il = nlast(i) 
        im1 = i - 1 
        bi = betas(nat(i)) 
        do k = if, il 
          ii = (k*(k - 1))/2 
          do j = 1, im1 
            jf = nfirst(j) 
            jl = nlast(j) 
            bj = betas(nat(j)) 
!  THE  +1.D-14 IS TO PREVENT POSSIBLE ERRORS IN THE DIAGONALIZATION.
            ij = ii + jf 
            h(ij) = 2.D0*h(ij)/(bi + bj) + 1.D-14 
            store(ij) = h(ij) 
            bj = betap(nat(j)) 
            bj = betap(nat(j)) 
            h(ii+jf+1:jl+ii) = 2.D0*h(ii+jf+1:jl+ii)/(bi + bj) + 1.D-14 
!  THE  +1.D-14 IS TO PREVENT POSSIBLE ERRORS IN THE DIAGONALIZATION.
            store(ii+jf+1:jl+ii) = h(ii+jf+1:jl+ii) 
          end do 
          store(ii+if:k+ii) = 0.D0 
          h(ii+if:k+ii) = 0.D0 
          bi = betap(nat(i)) 
        end do 
      end do 
      store(ifact(2:norbs+1)) = 1.D0 
      h(ifact(2:norbs+1)) = 1.D0                  
      call rsp (h, norbs, eig, vecs)   
      do i = 1, norbs 
        eig(i) = 1.D0/sqrt(abs(eig(i))) 
      end do 
      if (.not. allocated(work)) allocate(work(norbs,norbs))
      ij = 0 
      do i = 1, norbs 
        do j = 1, i 
          ij = ij + 1 
          sum = 0.D0 
          do k = 1, norbs 
            sum = sum + vecs(i+(k-1)*norbs)*eig(k)*vecs(j+(k-1)*norbs) 
          end do 
          work(i,j) = sum 
          work(j,i) = sum 
        end do 
      end do 
      call mult (c, work, vecs, norbs) 
      i = -1     
      nlower = (norbs*(norbs + 1))/2
      if (.not. lgpu) then
           call density_for_GPU (vecs, fract, nclose, nopen, 2.d0, nlower, norbs, 2, pb, 3)
      else
           call density_for_GPU (vecs, fract, nclose, nopen, 2.d0, nlower, norbs, 2, pb, 2)
      endif
!            
      pb = pb*store
      summ = 0.D0 
      do i = 1, norbs 
        sum = 0 
        do j = 1, i 
          sum = sum + pb(ifact(i)+j) 
        end do 
        do j = i + 1, norbs 
          sum = sum + pb(ifact(j)+i) 
        end do 
        summ = summ + sum 
        pb(ifact(i+1)) = sum 
      end do 
      h = store_h
      return  
      end subroutine mullik 
