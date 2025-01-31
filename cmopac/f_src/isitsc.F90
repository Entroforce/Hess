subroutine isitsc (escf, selcon, emin, iemin, iemax, okscf, niter, itrmax)
    use molkst_C, only: iscf
    use MOZYME_C, only : ovmax, energy_diff
    implicit none
    logical, intent (out) :: okscf
    integer, intent (in) :: itrmax, niter
    integer, intent (inout) :: iemax, iemin
    double precision, intent (in) :: emin, escf, selcon
    logical :: scf1 = .false.
    integer :: i, iemax1, iemin1
    double precision :: energy_test, fmo_test
    double precision, dimension (10) :: escf0
    data escf0 / 10 * 0.d0 /
   !
   ! Test the change in energy on successive iterations and the maximum
   ! element of the occ-vir block of the Fock matrix in the LMO basis
   !
    energy_test = selcon
    fmo_test = selcon * 5.0d0

    if (ovmax < fmo_test .and. Abs (energy_diff) < energy_test &
         & .and. scf1 .or. niter > itrmax) then
      okscf = .true.
      iscf = 2
      if (scf1) then
        iscf = 1
      end if
    else 
      if (ovmax < fmo_test .and. Abs (energy_diff) < energy_test) then
        scf1 = .true.
      else
        scf1 = .false.
      end if
      if (emin /= 0.d0) then
        !*****************************************************************
        !
        !  THE FOLLOWING TESTS ARE INTENDED TO ALLOW A FAST EXIT FROM 
        !  ITER IF THE RESULT IS 'GOOD ENOUGH' FOR THE CURRENT STEP IN 
        !  THE GEOMETRY OPTIMIZATION
        !
        if (escf < emin) then
          !
          !  THE ENERGY IS LOWER THAN THE PREVIOUS MINIMUM.  
          !  NOW CHECK THAT IT IS CONSISTENTLY LOWER.
          !
          iemax = 0
          iemin1 = iemin
          iemin = Min (5, iemin+1)
          if (iemin1 == 5) then
            do i = 2, 5
              escf0(i-1) = escf0(i)
            end do
          end if
          escf0(iemin) = escf
          !
          !  IS THE DIFFERENCE IN ENERGY BETWEEN TWO ITERATIONS LESS THAN 10%
          !  OF THE ENERGY GAIN FOR THIS GEOMETRY RELATIVE TO THE PREVIOUS
          !  MINIMUM.
          !
          if (iemin > 3) then
            do i = 2, iemin
              if (Abs (escf0(i)-escf0(i-1)) > 0.1d0*(emin-escf)) go to 1000
            end do
               !
               ! IS GOOD ENOUGH -- RAPID EXIT
               !
            okscf = .true.
            iscf = 1
            return
          end if
        else
            !
            !  THE ENERGY HAS RISEN ABOVE THAT OF THE PREVIOUS MINIMUM.
            !  WE NEED TO CHECK WHETHER THIS IS A FLUKE OR IS THIS REALLY
            !  A BAD GEOMETRY.
            !
          iemin = 0
          iemax1 = iemax
          iemax = Min (5, iemax+1)
          if (iemax1 == 5) then
            do i = 2, 5
              escf0(i-1) = escf0(i)
            end do
          end if
          escf0(iemax) = escf
          !
          !  IS THE DIFFERENCE IN ENERGY BETWEEN TWO ITERATIONS LESS THAN 10%
          !  OF THE ENERGY LOST FOR THIS GEOMETRY RELATIVE TO THE PREVIOUS
          !  MINIMUM.
          !
          if (iemax > 3) then
            do i = 2, iemax
              if (Abs (escf0(i)-escf0(i-1)) > 0.1d0*(escf-emin)) go to 1000
            end do
               !
               ! IS GOOD ENOUGH -- RAPID EXIT
               !
            okscf = .true.
            iscf = 1
            return
          end if
        end if
      end if
1000  okscf = .false.
    end if
end subroutine isitsc
