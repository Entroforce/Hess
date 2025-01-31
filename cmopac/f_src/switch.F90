      subroutine switch 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE parameters_C, only : alp, guess1, guess2, guess3, &
      betas, betap, betad, uss, upp, udd, zs, zp, zd, zsn, zpn, zdn, &
      gss, gsp, gpp, gp2, hsp, f0sd, g2sd, f0sd_store, g2sd_store, &
      polvol, pocord, xfac, alpb
!
!    
      USE molkst_C, only : keywrd, method_mndo, method_pm3, &
      & method_mndod, method_pm6, method_rm1, method_pm7, method_PM7_ts
!
!
      USE parameters_for_PM7_C, only :  uss7, upp7, udd7, zs7, zp7, zd7, betas7, &
        betap7, betad7, gss7, gsp7, gpp7, gp27, hsp7, polvo7, gues71, f0sd7, g2sd7, &
        gues72, gues73, poc_7, zsn7, zpn7, zdn7, alpb_and_xfac_pm7, alp7
!
!
       USE parameters_for_PM7_TS_C, only : uss7_TS, upp7_TS, udd7_TS, zs7_TS, zp7_TS, zd7_TS, betas7_TS, &
         betap7_TS, betad7_TS, gss7_TS, gsp7_TS, gpp7_TS, gp27_TS, hsp7_TS, polvo7_TS, poc_7_TS, &
         zsn7_TS, zpn7_TS, zdn7_TS, f0sd7_TS, g2sd7_TS, alp7_TS, gues7_TS1, gues7_TS2, gues7_TS3, &
         alpb_and_xfac_pm7_ts
!
!
!
      Use Parameters_for_PM7_Sparkles_C, only : gss7sp, alp7sp, gues7sp1, gues7sp2, gues7sp3
!
!
      USE chanel_C, only : iw 
!***********************************************************************
!
!   SWITCH copies data from the reference modules into the
!          arrays used by MOPAC.  The choice of reference data used
!          is given by the keywords MNDO, AM1, PM3, or MNDOD
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:03:42  03/15/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use mopend_I 
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
        pocord = 0.d0
        if (method_pm7_ts) then 
!
!    SWITCH IN PM7_TS PARAMETERS
!
          guess1 = gues7_TS1 
          guess2 = gues7_TS2 
          guess3 = gues7_TS3 
          zs = zs7_TS 
          zp = zp7_TS 
          zd = zd7_TS 
          zsn = zsn7_TS 
          zpn = zpn7_TS 
          zdn = zdn7_TS 
          uss = uss7_TS 
          upp = upp7_TS 
          udd = udd7_TS 
          betas = betas7_TS 
          betap = betap7_TS 
          betad = betad7_TS
          gss = gss7_TS 
          gpp = gpp7_TS 
          gsp = gsp7_TS 
          gp2 = gp27_TS 
          hsp = hsp7_TS 
          f0sd = f0sd7_TS
          g2sd = g2sd7_TS
          alp = alp7_TS
          pocord = poc_7_TS
          polvol = polvo7_TS
          call alpb_and_xfac_pm7_TS
        else if (method_pm7) then
!
!    SWITCH IN PM7 PARAMETERS
!
          guess1 = gues71 
          guess2 = gues72 
          guess3 = gues73 
          zs = zs7 
          zp = zp7 
          zd = zd7 
          zsn = zsn7 
          zpn = zpn7 
          zdn = zdn7 
          uss = uss7 
          upp = upp7 
          udd = udd7 
          betas = betas7 
          betap = betap7 
          betad = betad7
          gss = gss7 
          gpp = gpp7 
          gsp = gsp7 
          gp2 = gp27 
          hsp = hsp7 
          f0sd = f0sd7
          g2sd = g2sd7
          alp = alp7
          pocord = poc_7
          polvol = polvo7
          call alpb_and_xfac_pm7
          if (index(keywrd, " SPARK") /= 0) then
            do i = 1, 102
              if (alp7sp(i) > 0.1d0) then
                zd(i) = 0.d0
                zp(i) = 0.d0
                zs(i) = 0.d0
                zsn(i) = 0.d0
                zpn(i) = 0.d0
                zdn(i) = 0.d0
                uss(i) = 0.d0
                upp(i) = 0.d0
                udd(i) = 0.d0
                betas(i) = 0.d0
                betap(i) = 0.d0
                betad(i) = 0.d0
                alp(i) = alp7sp(i)
                gss(i) = gss7sp(i)
                gpp(i) = 0.d0
                gp2(i) = 0.d0
                hsp(i) = 0.d0
                gsp(i) = 0.d0
                f0sd(i) = 0.d0
                g2sd(i) = 0.d0
                pocord(i) = 0.d0
                guess1(i,1) = gues7sp1(i,1)
                guess2(i,1) = gues7sp2(i,1)
                guess3(i,1) = gues7sp3(i,1)
                guess1(i,2) = gues7sp1(i,2)
                guess2(i,2) = gues7sp2(i,2)
                guess3(i,2) = gues7sp3(i,2)
                guess1(i,3:4) = 0.d0
                guess2(i,3:4) = 0.d0
                guess3(i,3:4) = 0.d0               
              end if
            end do  
            alpb(:100,57:71) = 0.d0           
            alpb(57:71,:100) = 0.d0  
            xfac(:100,57:71) = 0.d0
            xfac(57:71,:100) = 0.d0 
          end if
        end if 
!
!    
        f0sd_store = f0sd
        g2sd_store = g2sd 

!
!  Symmetrize the alpb and xfac arrays
!
      do i = 1, 100
        alpb(:i,i) = alpb(i,:i)
        xfac(:i,i) = xfac(i,:i) 
      end do
      call fractional_metal_ion
      if (index(keywrd,' EXTERNAL')/=0) return  
      if (uss(1) > (-1.D0)) then 
        write (iw, &
      '(''  THE HAMILTONIAN REQUESTED IS NOT AVAILABLE IN THIS PROGRAM'')') 
        call mopend (&
          'THE HAMILTONIAN REQUESTED IS NOT AVAILABLE IN THIS PROGRAM') 
        return  
      endif 
      return  
      end subroutine switch 
      subroutine fractional_metal_ion
      use parameters_C, only : tore, zs, zp, gss, gpp, gp2, gsp, hsp, &
      betas, betap, upp, alp
      integer :: i 
!
!  Pseudo - halide metal ion with core charge of -1/2 - for use with
!  transition metal ions only
    tore(85) = -0.5d0
!
!  Pseudo - alkali metal ion with core charge of 1/2 - for use with
!  transition metal ions only
      tore(87) = 0.5d0
      do i = 85, 87, 2
        upp(i) = 0.d0
        alp(i) = 3.0d0
        zs(i) = 0.d0
        zp(i) = 0.d0
        betas(i) = 0.d0
        betap(i) = 0.d0
        gss(i) = 10.d0
        gsp(i) =0.d0
        gpp(i) = 0.d0
        gp2(i) = 0.d0
        hsp(i) = 0.d0
      end do
      end subroutine fractional_metal_ion

