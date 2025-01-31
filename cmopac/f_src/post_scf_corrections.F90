
subroutine post_scf_corrections(correction, l_grad)  
!
!    Add dispersion, hydrogen bonding, extra H-H repulsion, and other energies to
!    improve intermolecular interaction geometries and energies.
!
  use molkst_C, only : keywrd, E_disp, E_hb, E_hh, method_pm7, method_PM6, P_Hbonds, &
    method_pm7_hh, method_pm7_minus
  use common_arrays_C, only: dxyz
  implicit none
  double precision, intent(out) ::  correction
  logical, intent (in) :: l_grad
!
! Local variables
!
  logical ::  prt
  double precision, external :: & !                  Original references
                          !
  energy_corr_hh_rep,   & ! Rezac J., Hobza P., "Advanced Corrections of Hydrogen Bonding and 
                          ! Dispersion for Semiempirical Quantum Mechanical Methods", J. Chem. 
                          ! Theory and Comp 8, 141-151 (2012).
                          !

                          !
  PM6_DH_Dispersion,    & ! S. Grimme, "Accurate description of van der Waals complexes by 
                          ! density functional theory including empirical corrections."
                          ! J Comput Chem. 2004 Sep;25(12):1463-73.
                          !
  disp_DnX,             & ! (DH2X) Rezac and Hobza's correction: "A halogen-bonding correction  
                          ! for the semiempirical PM6 method" Chem. Phys. Lett. 506 286-289 (2011)
                          !
                          ! (D3H4X) Rezac J., Hobza P., "Advanced Corrections of Hydrogen Bonding 
                          ! and Dispersion for Semiempirical Quantum Mechanical Methods", J. Chem. 
                          ! Theory and Comp 8, 141-151 (2012)
                          !
  PM6_DH_H_bond_corrections !(DH2) Korth M., Pitonak M., Rezac J., Hobza P., "A Transferable 
                          ! H-bonding Correction for Semiempirical Quantum-Chemical Methods", 
                          ! J. Chem. Theory and Computation 6, 344-352 (2010).
                          !
                          ! (DH+) Korth M., "Third-Generation Hydrogen-Bonding Corrections for 
                          ! Semiempirical QM Methods and Force Fields", J. Chem. Theory Comput.
                          ! 6, 3808-3816 (2010).
!
! PM6-D3H4X: Brahmkshatriya, P. S.; Dobes, P.; Fanfrlik, J.; Rezac, J.; Paruch, K.; Bronowska, 
! A.; Lepsik, M.; Hobza, P. "Quantum Mechanical Scoring: Structural and Energetic Insights into 
! Cyclin-Dependent Kinase 2 Inhibition by Pyrazolo[1,5-a]pyrimidines" Curr. Comput.-Aid. Drug. 
! 2013 , 9 (1), 118�129.
!
  prt = (index(keywrd," 0SCF ") + index(keywrd," PRT ") /= 0 .and. index(keywrd," DISP") /= 0)
  correction = 0.d0
  E_hb       = 0.d0
  E_hh       = 0.d0
  E_disp     = 0.d0
  P_Hbonds   = 0
  if (method_pm7_hh) then
    correction = correction + energy_corr_hh_rep(l_grad, dxyz)
    correction = correction + PM6_DH_Dispersion(l_grad, dxyz)
    correction = correction + PM6_DH_H_bond_corrections(l_grad, prt)
  else if (method_pm7_minus) then
    return
  else if (method_pm7) then
    correction = correction + PM6_DH_Dispersion(l_grad, dxyz)
    correction = correction + PM6_DH_H_bond_corrections(l_grad, prt)
  end if
  if (index(keywrd, " SILENT") == 0) then
    if (prt .and. P_Hbonds > 0) call print_post_scf_corrections
  end if
  return
end subroutine post_scf_corrections

