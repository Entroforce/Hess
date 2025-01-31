    subroutine calpar 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double  
      USE parameters_C, only : ios, iop, iod, qq, am, ad, aq, dd, &
      gpp, gp2, hsp, gss, gsp, &
      zs, zp, uss, upp, udd, eisol, &
      f0sd, g2sd, f0sd_store, g2sd_store
      USE funcon_C, only : ev
      USE molkst_C, only : keywrd, method_pm7, in_house_only
!...Translated by Pacific-Sierra Research 77to90  4.4G  08:24:24  03/10/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
     
      integer , dimension(107) :: nspqn 
      integer :: i, k, l, jmax, j
      character :: num*1, file_name*120
      real(double), dimension(107) :: gssc, gspc, hspc, gp2c, gppc 
      real(double) :: p, p4,  hpp, qn, gdd1, d1, d2, df, hsp1, hsp2, d3, &
        gqq, q1, q2, qf, hpp1, hpp2, q3
      save nspqn 
!----------------------------------------------- 
      data nspqn/ 2*1, 8*2, 8*3, 18*4, 18*5, 32*6, 21*0/  
!
! THE CONTINUATION LINES INDICATE THE PRINCIPAL QUANTUM NUMBER.
!
!
!     SET SCALING PARAMETER.
      p = 2.D0 
      p4 = p**4
      call sp_two_electron 
      am = 0.d0
      do i = 2, 97 
!
!  GSSC is the number of two-electron terms of type <SS|SS>
!
        gssc(i) = max(ios(i)-1,0) 
        k = iop(i) 
!
!  GSPC is the number of two-electron terms of type <SS|PP>
!
        gspc(i) = ios(i)*k 
        l = min(k,6 - k) 
!
!  GP2C is the number of two-electron terms of type <PP|PP>
!       plus 0.5 of the number of HPP integrals.
!  (HPP is not used; instead it is replaced by 0.5(GPP-GP2))
!
        gp2c(i) = (k*(k - 1))/2 + 0.5D0*(l*(l - 1))/2 
!
!  GPPC is minus 0.5 times the number of HPP integrals.
!
        gppc(i) = -0.5D0*(l*(l - 1))/2 
!
!  HSPC is the number of two-electron terms of type <SP|SP>.
!       (S and P must have the same spin.  In all cases, if
!  P is non-zero, there are two S electrons)
!
        hspc(i) = -k*ios(i)*0.5d0
!
!
        if (zp(i)<1.D-4 .and. zs(i)<1.D-4) cycle  
!*********************************************************************
!
!   CONSTRAINTS ON THE POSSIBLE VALUES OF PARAMETERS
!
!*********************************************************************
        zp(i) = dmax1(0.3D0,zp(i)) 
!  PUT IN ANY CONSTRAINTS AT THIS POINT
!*********************************************************************
        hpp = 0.5D0*(gpp(i)-gp2(i)) 
        hpp = max(0.1D0,hpp) 
        hsp(i) = max(1.D-7,hsp(i)) 
        eisol(i) = uss(i)*ios(i) + upp(i)*iop(i) + udd(i)*iod(i) + gss(i)*gssc(i) + &
          gpp(i)*gppc(i) + gsp(i)*gspc(i) + gp2(i)*gp2c(i) + hsp(i)*hspc(i) 
        qn = nspqn(i) 
        dd(i) = (2.D0*qn + 1)*(4.D0*zs(i)*zp(i))**(qn + 0.5D0)/(zs(i)+zp(i))**(2.D0*qn + 2)/sqrt(3.D0) 
        qq(i) = sqrt((4.D0*qn*qn + 6.D0*qn + 2.D0)/20.D0)/zp(i) 
!     CALCULATE ADDITIVE TERMS, IN ATOMIC UNITS.
        jmax = 5 
        gdd1 = (hsp(i)/(ev*dd(i)**2))**(1.D0/3.D0) 
        d1 = gdd1 
        d2 = gdd1 + 0.04D0 
        do j = 1, jmax 
          df = d2 - d1 
          hsp1 = 0.5D0*d1 - 0.5D0/sqrt(4.D0*dd(i)**2+1.D0/d1**2) 
          hsp2 = 0.5D0*d2 - 0.5D0/sqrt(4.D0*dd(i)**2+1.D0/d2**2) 
          if (abs(hsp2 - hsp1) < 1.D-25) exit  
          d3 = d1 + df*(hsp(i)/ev-hsp1)/(hsp2 - hsp1) 
          d1 = d2 
          d2 = d3 
        end do 
        gqq = (p4*hpp/(ev*48.D0*qq(i)**4))**0.2D0 
        q1 = gqq 
        q2 = gqq + 0.04D0 
        do j = 1, jmax 
          qf = q2 - q1 
          hpp1 = 0.25D0*q1 - 0.5D0/sqrt(4.D0*qq(i)**2+1.D0/q1**2) + 0.25D0/&
            sqrt(8.D0*qq(i)**2+1.D0/q1**2) 
          hpp2 = 0.25D0*q2 - 0.5D0/sqrt(4.D0*qq(i)**2+1.D0/q2**2) + 0.25D0/&
            sqrt(8.D0*qq(i)**2+1.D0/q2**2) 
          if (abs(hpp2 - hpp1) < 1.D-25) exit  
          q3 = q1 + qf*(hpp/ev - hpp1)/(hpp2 - hpp1) 
          q1 = q2 
          q2 = q3 
        end do 
        am(i) = gss(i)/ev 
        ad(i) = d2 
        aq(i) = q2 
      end do 
      do i = 1, 107 
        if (am(i) < 1.d-20) then
          if (gss(i) > 1.d-20) then
            am(i) = gss(i)/ev
          else
            am(i) = 1.d0
          end if
        end if
      end do 
      eisol(1) = uss(1) 
      am(1) = gss(1)/ev 
      ad(1) = am(1) 
      aq(1) = am(1)
      do i = 1, 100
        if (f0sd_store(i) < 1.d-20) f0sd(i) = 0.d0 ! Force f0sd and g2sd to zero if not already defined
        if (g2sd_store(i) < 1.d-20) g2sd(i) = 0.d0
      end do
      
      call inid        ! Calculate derived parameters for "d" orbital work
!
!   Atomic number 102 is the capped bond.  It should have a very
!   small AM to prevent division by zero in REPPD, and to avoid
!   spurious effects due to normal AM values.
!
      am(102) = 1.D-10 
!
!     DEBUG PRINTING.
!     THIS IS FORMATTED FOR DIRECT INSERTION INTO 'PARAM'
!
      if (index(keywrd,' DEP ') == 0) return  
      if ( .not. in_house_only) return      
      num = "6"
      if (method_pm7) num = "7"
   !   file_name = "parameters_for_PM"//num//"_C.F90"
      file_name = "parameters_for_PM7_TS_C.F90"
      call create_parameters_for_PMx_C(file_name, num)       
      return  
    end subroutine calpar 
!
!
!
!
  subroutine create_parameters_for_PMx_C(file_name, num)
!
    USE molkst_C, only : keywrd
!
    USE parameters_C, only : guess1, guess2, guess3, gpp, gp2, hsp, gss, gsp, betas, betap, betad, &
    zs, zp, zd, uss, upp, udd, zsn, zpn, zdn, pocord, alpb, xfac, &
    f0sd, g2sd, main_group, alp, polvol
!
    USE chanel_C, only : iw, param_out
!
    USE elemts_C, only : atom_names
!
    implicit none
    character :: file_name*120, num*1, num1*1
!
!  Local
!
    integer :: i, j,k
    logical :: new
    double precision :: abond, fff
    if (index(keywrd,"PKA") == 0) then
        write (iw, *)" A new file called '"//trim(file_name)//"' will be written"
        j = 0
   97   open (unit=param_out, file=trim(file_name), status="UNKNOWN", iostat=i)
        if ( i /= 0 .and. j < 10)then
!
!  The file can exist, but is not currently accessible
!
          write(iw,*) "Problem writing '"//trim(file_name)//"'"
          j = j + 1
          goto 97
        end if
        if (j < 10) then
          rewind (param_out)
          write (param_out,"(a)") "  module Parameters_for_PM7_TS_C"        
          write (param_out,"(a)") "    use vast_kind_param, ONLY:  double  "
          write (param_out,"(a)") &
          & "    real(double), dimension(107) :: uss7_TS, upp7_TS, udd7_TS, zs7_TS, zp7_TS, zd7_TS, betas7_TS, &"
          write (param_out,"(a)") "    betap7_TS, betad7_TS, gss7_TS, gsp7_TS, gpp7_TS, gp27_TS, hsp7_TS, polvo7_TS, poc_7_TS, &"
          write (param_out,"(a)") "    zsn7_TS, zpn7_TS, zdn7_TS, f0sd7_TS, g2sd7_TS, alp7_TS"
          write (param_out,"(a)") "    real(double), dimension(107,4) :: gues7_TS1, gues7_TS2, gues7_TS3" 
          do i = 1, 107 
            if (zs(i) < 1.d-20 .and. gss(i) < 1.d-20 .and. Abs(guess1(i,1)) &
            + Abs(alp(i)) < 1.d-20) cycle 
            write (param_out,"('!')") 
            write (param_out, '("!",20X,"Data for Element ",I3,5x,a)') i, atom_names(i)
            write (param_out,"('!')") 
            if (.false.) then
            if (uss(i) /= 0.D0) write (param_out, &
              '(6X,"data     uss'//num//'(",I3,")/",F17.6,"D0/")') i, uss(i) 
            if (upp(i) /= 0.D0) write (param_out, &
              '(6X,"data     upp'//num//'(",I3,")/",F17.6,"D0/")') i, upp(i) 
            if (udd(i) /= 0.D0) write (param_out, &
              '(6X,"data     udd'//num//'(",I3,")/",F17.6,"D0/")') i, udd(i) 
            if (betas(i) /= 0.D0) write (param_out, &
              '(6X,"data   betas'//num//'(",I3,")/",F17.6,"D0/")') i, betas(i) 
            if (betap(i) /= 0.D0) write (param_out, &
              '(6X,"data   betap'//num//'(",I3,")/",F17.6,"D0/")') i, betap(i) 
            if (betad(i) /= 0.D0) write (param_out, &
              '(6X,"data   betad'//num//'(",I3,")/",F17.6,"D0/")') i, betad(i) 
            if (zs(i) /= 0.D0) write (param_out, &
              '(6X,"data      zs'//num//'(",I3,")/",F17.6,"D0/")') i, zs(i) 
            if (zp(i) /= 0.D0) write (param_out, &
              '(6X,"data      zp'//num//'(",I3,")/",F17.6,"D0/")') i, zp(i) 
            if (zd(i) /= 0.D0) write (param_out, &
              '(6X,"data      zd'//num//'(",I3,")/",F17.6,"D0/")') i, zd(i) 
            if (zsn(i) /= 0.D0) write (param_out, &
              '(6X,"data     zsn'//num//'(",I3,")/",F17.6,"D0/")') i, zsn(i) 
            if (zpn(i) /= 0.D0) write (param_out, &
              '(6X,"data     zpn'//num//'(",I3,")/",F17.6,"D0/")') i, zpn(i) 
            if (zdn(i) /= 0.D0) write (param_out, &
              '(6X,"data     zdn'//num//'(",I3,")/",F17.6,"D0/")') i, zdn(i) 
            if (alp(i) /= 0.D0) write (param_out, &
              '(6X,"data     alp'//num//'(",I3,")/",F17.6,"D0/")') i, alp(i) 
            if (gss(i) /= 0.D0) write (param_out, &
              '(6X,"data     gss'//num//'(",I3,")/",F17.6,"D0/")') i, gss(i) 
            if (gsp(i) /= 0.D0) write (param_out, &
              '(6X,"data     gsp'//num//'(",I3,")/",F17.6,"D0/")') i, gsp(i) 
            if (gpp(i) /= 0.D0) write (param_out, &
              '(6X,"data     gpp'//num//'(",I3,")/",F17.6,"D0/")') i, gpp(i) 
            if (gp2(i) /= 0.D0) write (param_out, &
              '(6X,"data     gp2'//num//'(",I3,")/",F17.6,"D0/")') i, gp2(i) 
            if (hsp(i) /= 0.D0) write (param_out, &
              '(6X,"data     hsp'//num//'(",I3,")/",F17.6,"D0/")') i, hsp(i) 
            if (pocord(i) /= 0.D0) write (param_out, &
              '(6X,"data    poc_'//num//'(",I3,")/",F17.6,"D0/")') i, pocord(i) 
            if (polvol(i) /= 0.D0) write (param_out, &
              '(6X,"data   polvo'//num//'(",I3,")/",F17.6,"D0/")') i, polvol(i) 
            if (.not. main_group(i) .and. f0sd(i) /= 0.D0) write (param_out, &
              '(6X,"data    f0sd'//num//'(",I3,")/",F17.6,"D0/")') i, f0sd(i) 
            if (.not. main_group(i) .and. g2sd(i) /= 0.D0) write (param_out, &
              '(6X,"data    g2sd'//num//'(",I3,")/",F17.6,"D0/")') i, g2sd(i) 
            if (i /= 99 .and. i /= 98) then
              do j = 1, 4 
                if (guess1(i,j) /= 0.D0) &
                write (param_out, &
                '(6X,"data gues'//num//'1(",I3,",",I1,")/",          F17.6,"D0/")')&
                 i, j, guess1(i,j) 
                if (guess2(i,j) /= 0.D0) &
                write (param_out, &
                '(6X,"data gues'//num//'2(",I3,",",I1,")/",          F17.6,"D0/")')&
                 i, j, guess2(i,j) 
                if (guess3(i,j) /= 0.D0) &
                write (param_out, &
                '(6X,"data gues'//num//'3(",I3,",",I1,")/",          F17.6,"D0/")')&
                 i, j, guess3(i,j) 
              end do 
            else
              do j = 1, 4                 
                if (guess1(i,j) /= 0.D0) then
                  k = 1 + 3*(j - 1)
                  if (i == 98) k = k + 12
                  if (k < 10) then
                    num1 = "1"
                  else
                    num1 = "2"
                  end if
                  write (param_out, &
                    '(6X,"data gues'//num//'1(",I3,",",I1,")/",          F17.6,"D0/ ! = PAR",i'//num1//')')&
                     i, j, guess1(i,j),k
                 end if
                 if (guess2(i,j) /= 0.D0) then
                  k = 2 + 3*(j - 1)
                  if (i == 98) k = k + 12
                  if (k < 10) then
                    num1 = "1"
                  else
                    num1 = "2"
                  end if
                  write (param_out, &
                    '(6X,"data gues'//num//'2(",I3,",",I1,")/",          F17.6,"D0/ ! = PAR",i'//num1//')')&
                     i, j, guess2(i,j),k
                 end if
                 if (guess3(i,j) /= 0.D0) then
                  k = 3 + 3*(j - 1)
                  if (i == 98) k = k + 12
                  if (k < 10) then
                    num1 = "1"
                  else
                    num1 = "2"
                  end if
                  write (param_out, &
                    '(6X,"data gues'//num//'3(",I3,",",I1,")/",          F17.6,"D0/ ! = PAR",i'//num1//')')&
                     i, j, guess3(i,j),k
                 end if
               
              end do 
            end if
            else


            if (uss(i) /= 0.D0) write (param_out, &
              '(6X,"data     uss7_TS(",I3,")/",F17.6,"D0/")') i, uss(i) 
            if (upp(i) /= 0.D0) write (param_out, &
              '(6X,"data     upp7_TS(",I3,")/",F17.6,"D0/")') i, upp(i) 
            if (udd(i) /= 0.D0) write (param_out, &
              '(6X,"data     udd7_TS(",I3,")/",F17.6,"D0/")') i, udd(i) 
            if (betas(i) /= 0.D0) write (param_out, &
              '(6X,"data   betas7_TS(",I3,")/",F17.6,"D0/")') i, betas(i) 
            if (betap(i) /= 0.D0) write (param_out, &
              '(6X,"data   betap7_TS(",I3,")/",F17.6,"D0/")') i, betap(i) 
            if (betad(i) /= 0.D0) write (param_out, &
              '(6X,"data   betad7_TS(",I3,")/",F17.6,"D0/")') i, betad(i) 
            if (zs(i) /= 0.D0) write (param_out, &
              '(6X,"data      zs7_TS(",I3,")/",F17.6,"D0/")') i, zs(i) 
            if (zp(i) /= 0.D0) write (param_out, &
              '(6X,"data      zp7_TS(",I3,")/",F17.6,"D0/")') i, zp(i) 
            if (zd(i) /= 0.D0) write (param_out, &
              '(6X,"data      zd7_TS(",I3,")/",F17.6,"D0/")') i, zd(i) 
            if (zsn(i) /= 0.D0) write (param_out, &
              '(6X,"data     zsn7_TS(",I3,")/",F17.6,"D0/")') i, zsn(i) 
            if (zpn(i) /= 0.D0) write (param_out, &
              '(6X,"data     zpn7_TS(",I3,")/",F17.6,"D0/")') i, zpn(i) 
            if (zdn(i) /= 0.D0) write (param_out, &
              '(6X,"data     zdn7_TS(",I3,")/",F17.6,"D0/")') i, zdn(i) 
            if (alp(i) /= 0.D0) write (param_out, &
              '(6X,"data     alp7_TS(",I3,")/",F17.6,"D0/")') i, alp(i) 
            if (gss(i) /= 0.D0) write (param_out, &
              '(6X,"data     gss7_TS(",I3,")/",F17.6,"D0/")') i, gss(i) 
            if (gsp(i) /= 0.D0) write (param_out, &
              '(6X,"data     gsp7_TS(",I3,")/",F17.6,"D0/")') i, gsp(i) 
            if (gpp(i) /= 0.D0) write (param_out, &
              '(6X,"data     gpp7_TS(",I3,")/",F17.6,"D0/")') i, gpp(i) 
            if (gp2(i) /= 0.D0) write (param_out, &
              '(6X,"data     gp27_TS(",I3,")/",F17.6,"D0/")') i, gp2(i) 
            if (hsp(i) /= 0.D0) write (param_out, &
              '(6X,"data     hsp7_TS(",I3,")/",F17.6,"D0/")') i, hsp(i) 
            if (pocord(i) /= 0.D0) write (param_out, &
              '(6X,"data    poc_7_TS(",I3,")/",F17.6,"D0/")') i, pocord(i) 
            if (polvol(i) /= 0.D0) write (param_out, &
              '(6X,"data   polvo7_TS(",I3,")/",F17.6,"D0/")') i, polvol(i) 
            if (.not. main_group(i) .and. f0sd(i) /= 0.D0) write (param_out, &
              '(6X,"data    f0sd7_TS(",I3,")/",F17.6,"D0/")') i, f0sd(i) 
            if (.not. main_group(i) .and. g2sd(i) /= 0.D0) write (param_out, &
              '(6X,"data    g2sd7_TS(",I3,")/",F17.6,"D0/")') i, g2sd(i) 
            if (i /= 99 .and. i /= 98) then
              do j = 1, 4 
                if (guess1(i,j) /= 0.D0) &
                write (param_out, &
                '(6X,"data gues7_TS1(",I3,",",I1,")/",          F17.6,"D0/")')&
                 i, j, guess1(i,j) 
                if (guess2(i,j) /= 0.D0) &
                write (param_out, &
                '(6X,"data gues7_TS2(",I3,",",I1,")/",          F17.6,"D0/")')&
                 i, j, guess2(i,j) 
                if (guess3(i,j) /= 0.D0) &
                write (param_out, &
                '(6X,"data gues7_TS3(",I3,",",I1,")/",          F17.6,"D0/")')&
                 i, j, guess3(i,j) 
              end do 
            else
              do j = 1, 4                 
                if (guess1(i,j) /= 0.D0) then
                  k = 1 + 3*(j - 1)
                  if (i == 98) k = k + 12
                  if (k < 10) then
                    num1 = "1"
                  else
                    num1 = "2"
                  end if
                  write (param_out, &
                    '(6X,"data gues7_TS1(",I3,",",I1,")/",          F17.6,"D0/ ! = PAR",i'//num1//')')&
                     i, j, guess1(i,j),k
                 end if
                 if (guess2(i,j) /= 0.D0) then
                  k = 2 + 3*(j - 1)
                  if (i == 98) k = k + 12
                  if (k < 10) then
                    num1 = "1"
                  else
                    num1 = "2"
                  end if
                  write (param_out, &
                    '(6X,"data gues7_TS2(",I3,",",I1,")/",          F17.6,"D0/ ! = PAR",i'//num1//')')&
                     i, j, guess2(i,j),k
                 end if
                 if (guess3(i,j) /= 0.D0) then
                  k = 3 + 3*(j - 1)
                  if (i == 98) k = k + 12
                  if (k < 10) then
                    num1 = "1"
                  else
                    num1 = "2"
                  end if
                  write (param_out, &
                    '(6X,"data gues7_TS3(",I3,",",I1,")/",          F17.6,"D0/ ! = PAR",i'//num1//')')&
                     i, j, guess3(i,j),k
                 end if
               
              end do 
            end if


            end if
          end do 
          write (param_out,"(a)") '  contains'  
    
          write (param_out,"(a)") &
          & '  subroutine alpb_and_xfac_pm7_TS', &
          & '    use parameters_C, only : xfac, alpb'
    !
    !  Write out all the diatomic parameters
    !
          do i = 1, 100     
            new = .true.
            do j = 1, i
              abond = alpb(i,j)
              if (abond > 1.d-4) then
                fff = xfac(i,j)
                if ( new ) then
                new = .false.
                write (param_out,*)"!"
                end if
                 write(param_out,"(6x,a,i2,a,i2,a,f12.6,a)") &
                 "alpb(",i,",",j,") = ", abond,"d0 !"//atom_names(i)//" - "//atom_names(j)
                 write(param_out,"(6x,a,i2,a,i2,a,f12.6,a)") &
                 "xfac(",i,",",j,") = ", fff,"d0 !"//atom_names(i)//" - "//atom_names(j)
              end if
            end do
          end do
          write (param_out,"(a)") '    end subroutine alpb_and_xfac_pm7_TS'
          write (param_out,"(a)") "  end module Parameters_for_PM7_TS_C"   
        end if
      else
        if (Abs(uss(99)) < 0.1d0) return
        j = 0
        write (iw, *)" A new file called 'parameters_for_PKA.F90' will be written"
   98   open (unit=param_out, file="parameters_for_PKA.F90", status="UNKNOWN", iostat=i)
        if ( i /= 0 .and. j < 10)then
!
!  The file can exist, but is not currently accessible
!
          write(iw,*) "Problem with parameters_for_PKA.F90"
          j = j + 1
          goto 98
        end if
        if (j < 10) then
          rewind (param_out)
          write (param_out,"(a)") "  subroutine Parameters_for_PKA(c1, c2, c3)"        
          write (param_out,"(a)") "    use parameters_C, only : uss, upp, udd, zs, zp, zd, betas, &"
          write (param_out,"(a)") "    betap, betad, gss, gsp, gpp, gp2, hsp, zsn, zpn, zdn"
          write (param_out,"(a)") "    double precision, intent (out) :: c1, c2, c3"
          write (param_out,'(4x,"c1   =",F16.7,"d0")')  uss(99) 
          write (param_out,'(4x,"c2   =",F16.7,"d0")')  upp(99) 
          write (param_out,'(4x,"c3   =",F16.7,"d0")')  zs(99) 
          do i = 1, 100 
            select case(i)
            case (1, 6, 7, 8, 9, 17, 35, 53)      
              new = .true.
              if ( new ) then
                new = .false.
                write (param_out,"('!')") 
                write (param_out, '("!",20X,"Data for pKa for Element",I3,5x,a)') i, atom_names(i)
                write (param_out,"('!')") 
              end if
              if (uss(i) /= 0.d0)   write (param_out,'(4x,"uss  (",i2,")   =",F16.7,"d0")') i, uss(i) 
              if (upp(i) /= 0.d0)   write (param_out,'(4x,"upp  (",i2,")   =",F16.7,"d0")') i, upp(i) 
              if (udd(i) /= 0.d0)   write (param_out,'(4x,"udd  (",i2,")   =",F16.7,"d0")') i, udd(i) 
              if (betas(i) /= 0.d0) write (param_out,'(4x,"betas(",i2,")   =",F16.7,"d0")') i, betas(i) 
              if (betap(i) /= 0.d0) write (param_out,'(4x,"betap(",i2,")   =",F16.7,"d0")') i, betap(i) 
              if (betad(i) /= 0.d0) write (param_out,'(4x,"betad(",i2,")   =",F16.7,"d0")') i, betad(i) 
              if (zs(i) /= 0.d0)    write (param_out,'(4x,"zs   (",i2,")   =",F16.7,"d0")') i, zs(i)
              if (zp(i) /= 0.d0)    write (param_out,'(4x,"zp   (",i2,")   =",F16.7,"d0")') i, zp(i) 
              if (zd(i) /= 0.d0)    write (param_out,'(4x,"zd   (",i2,")   =",F16.7,"d0")') i, zd(i)  
              if (zsn(i) /= 0.d0)   write (param_out,'(4x,"zsn  (",i2,")   =",F16.7,"d0")') i, zsn(i)
              if (zpn(i) /= 0.d0)   write (param_out,'(4x,"zpn  (",i2,")   =",F16.7,"d0")') i, zpn(i) 
              if (zdn(i) /= 0.d0)   write (param_out,'(4x,"zdn  (",i2,")   =",F16.7,"d0")') i, zdn(i)  
              if (gss(i) /= 0.d0)   write (param_out,'(4x,"gss  (",i2,")   =",F16.7,"d0")') i, gss(i) 
              if (gsp(i) /= 0.d0)   write (param_out,'(4x,"gsp  (",i2,")   =",F16.7,"d0")') i, gsp(i) 
              if (gpp(i) /= 0.d0)   write (param_out,'(4x,"gpp  (",i2,")   =",F16.7,"d0")') i, gpp(i) 
              if (gp2(i) /= 0.d0)   write (param_out,'(4x,"gp2  (",i2,")   =",F16.7,"d0")') i, gp2(i)
              if (hsp(i) /= 0.d0)   write (param_out,'(4x,"hsp  (",i2,")   =",F16.7,"d0")') i, hsp(i)   
            end select
          end do
          write(param_out,"(a)")  "    return"
          write (param_out,"(a)") "  end subroutine Parameters_for_PKA"   
          close(param_out)
          call mopend("Job stopped because a pKa FORTRAN file has been created")
        end if
      end if
    return   
  end subroutine create_parameters_for_PMx_C
