    module radii_C
   double precision :: covalent_radii(118)
   data covalent_radii /0.37d0, 0.32d0, 1.34d0, 0.9d0, 0.82d0, 0.77d0, &
     0.75d0, 0.73d0, 0.71d0, 0.69d0, 1.54d0, 1.3d0, 1.18d0, 1.11d0, 1.06d0, 1.02d0, 0.99d0, 0.97d0, &
     1.96d0, 1.74d0, 1.44d0, 1.36d0, 1.25d0, 1.27d0, 1.39d0, 1.25d0, 1.26d0, 1.21d0, 1.38d0, 1.31d0, &
     1.26d0, 1.22d0, 1.19d0, 1.16d0, 1.14d0, 1.1d0, 2.11d0, 1.92d0, 1.62d0, 1.48d0, 1.37d0, 1.45d0, &
     1.56d0, 1.26d0, 1.35d0, 1.31d0, 1.53d0, 1.48d0, 1.44d0, 1.41d0, 1.38d0, 1.35d0, 1.33d0, 1.3d0, &
     2.25d0, 1.98d0, 1.69d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
     0.0d0, 0.0d0, 1.6d0, 1.5d0, 1.38d0, 1.46d0, 1.59d0, 1.28d0, 1.37d0, 1.28d0, 1.44d0, 1.49d0, 0.0d0, &
     0.0d0, 1.46d0, 0.0d0, 0.0d0, 1.45d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0/
   end module radii_C

  function energy_corr_hh_rep(l_grad, dxyz) 
    use common_arrays_C, only: nat, coord
    use molkst_C, only : numat, keywrd
    use chanel_C, only : iw
    use elemts_C, only: elemnt
    implicit none
    logical, intent (in) :: l_grad
    double precision, intent (inout) :: dxyz(3, numat)
    double precision :: grad_hh(3, numat)
    double precision :: e_corr, e_corr_sum, r, d_rad, g(3), energy_corr_hh_rep
    integer :: i, j ! iteration counters
    double precision, external :: distance, poly
    e_corr_sum = 0.d0
    grad_hh = 0.d0
!
! Iterate over H atoms twice
!
    do i = 1, numat
      if (nat(i) /= 1) cycle
      do j = 1, i - 1
        if (nat(j) /= 1) cycle
!
! Calculate distance
!
        r = distance(i,j)
        e_corr = poly(r, l_grad, d_rad)
        e_corr_sum = e_corr_sum + e_corr
        if (l_grad) then
!
! Cartesian components of the gradient
!
          g(:) = -(coord(:,i) - coord(:,j))/r*d_rad
!
! Add pair contribution to the global gradient
!  
          grad_hh(:,i) = grad_hh(:,i) - g(:)
          grad_hh(:,j) = grad_hh(:,j) + g(:)
        end if
      end do
    end do
    energy_corr_hh_rep = e_corr_sum
    if (l_grad) then
      dxyz = dxyz + grad_hh
      if (index(keywrd, " DERIV") > 0) then
        write (iw, '(/25X,a)')"HH REPULSION" 
        write (iw, '(3X,       ''NUMBER  ATOM '',5X,''X'',12X,''Y'',12X,''Z'',/)') 
        write (iw, '(I6,4x,a2,F13.6,2F13.6)') (i, elemnt(nat(i)), grad_hh(:,i), i = 1,numat) 
      end if 
    end if
    return
  end function energy_corr_hh_rep
  double precision function poly(r, l_grad, dpoly)
    implicit none
    double precision :: r, dpoly 
    logical :: l_grad
!
!  poly is the hydrogen-hydrogen correction for two atoms separated by "r" Angstroms
!  dpoly is the gradient of the internal coordinate
!
    if (r <= 1.d0) then
      poly = 25.46293603147693d0
      dpoly = 0.d0
    else if ( r > 1.d0 .and. r < 1.5d0) then      
      poly =  -2714.952351603469651d0 * r**5 &
             +17103.650110591705015d0 * r**4 &
             -42511.857982217959943d0 * r**3 &
             +52063.196799138342612d0 * r**2 &
             -31430.658335972289933d0 * r    &
              +7516.084696095140316d0
      if (l_grad) &
        dpoly = -2714.952351603469651d0*5.d0 * r**4 &
               +17103.650110591705015d0*4.d0 * r**3 &
               -42511.857982217959943d0*3.d0 * r**2 &
               +52063.196799138342612d0*2.d0 * r    &
               -31430.658335972289933d0 
    else
      poly = 118.7326d0*exp(-1.53965d0*(r**1.72905d0))  
      if (l_grad) &
      dpoly = -1.53965d0*1.72905d0*r**0.72905d0*118.7326d0*exp(-1.53965d0*(r**1.72905d0))
    end if      
    return
  end function poly

