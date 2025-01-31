module molkst_C
  USE vast_kind_param, ONLY:  double 
!
!  This module contains all the scalars relating to the system being calculated. 
!  Every entry is unique, and has the same meaning in every subroutine.
!
!
!  This module can also be regarded as a dictionary of the scalars relating to the system.
!  Quantities which have only a one-line descriptionis are of no interest outside MOPAC.
!  The data-type for each quantity is given at the start of the statement
!
  integer, bind(C) :: &
  &  maxatoms, & !  Maximum number of atoms allowed = number of lines in the data set
                 !
  &  natoms,   & !  Term        Number of atoms, real plus dummy, etc., in the system
                 !  Units       Atoms
                 !  Min. value  0
                 !  Max. value  very large
                 !
  &  numat = 0,& !  Term        Number of real atoms in the system
                 !  Units       Atoms
                 !  Min. value  0
                 !  Max. value  very large
                 !
  &  numat_old,& !  Number of real atoms in the system supplied by the data set
  &  norbs,    & !  Term        Number of atomic orbitals in the system
                 !  Units       atomic orbitals
                 !  Min. value  0
                 !  Max. value  very large
                 !
  &  nelecs,   & !  Term        Number of electrons
                 !  Units       Electrons
                 !  Min. value  0
                 !  Max. value  2*norbs    
                 !             
  &  ispd,     & !  0 if a s-p basis set, 1 or more if any atoms have d-orbitals
  &  ndep = 0, & !  Number of dependent coordinates, related by symmetry
  &  nvar = 0, & !  Number of coordinate that are to be optimized
                 !
  &  nclose,   & !  Term        Number of doubly-occupied M.O.s
                 !  Units       Molecular orbitals
                 !  Min. value  0
                 !  Max. value  nelecs/2  
                 !
  &  nopen,    & !  Term        Upper bound of active space in ROHF - C.I.
                 !  Units       Molecular orbitals
                 !  Definition  nclose + number of M.O.s in the active space.
                 !  Min. value  0
                 !  Max. value  norbs     
                 !  Default     nclose
                 !  
  &  nalpha,   & !  Term        Number of alpha singly-occupied M.O.s
                 !  Units       Molecular orbitals
                 !  Min. value  0
                 !  Max. value  min(norbs, nelecs)
                 !  Default     0 (RHF) (nelecs + 1)/2 (UHF)
                 !
  &  nbeta,    & !  Term        Number of beta singly-occupied M.O.s
                 !  Units       Molecular orbitals
                 !  Min. value  0
                 !  Max. value  min(norbs, nalpha - nelecs)
                 !  Default     0 (RHF) nelecs/2 (UHF)
                 !
  &nalpha_open,& !  Term        Number of alpha fractionally-occupied M.O.s
                 !  Units       Molecular orbitals
                 !  Min. value  0
                 !  Max. value  norbs
                 !  Default     nalpha
                 !
  & nbeta_open,& !  Term        Number of beta fractionally-occupied M.O.s
                 !  Units       Molecular orbitals
                 !  Min. value  0
                 !  Max. value  norbs
                 !  Default     nbeta
                 !
  &  numcal=0, & !  Number of the calculation.  The first molecule calculated will have
                 !  numcal = 1.  The second molecule calculated will have numcal = 2, and so on.
                 !
   step_num=0, & !  The unique number of the stage or step within a calculation of a molecule.
                 !  The first stage of the first calculation will have step_num = 1. 
                 !  If that calculation involves a second stage, then step_num will be incremented.
                 !  When a new calculation starts, step_num will be incremented.
                 !  Each stage is limited to the same electronic structure.  Most calculations
                 !  will only have one stage, e.g. geometry optimization or force constants.
                 !
  &  mpack,    & !  Number of elements in a lower-half-triangle = (norbs*(norbs+1))/2
  &  n2elec,   & !  Number of two-electron integrals
  &  nscf,     & !  Number of SCF calculations done
  &  iscf,     & !  Index of message of how the SCF ended
  &  iflepo,   & !  Index of message of how the geometry operation ended
  &  maxtxt,   & !  Maximum number of characters in labeled atoms (e.g. C(on Si))
  &  last,     & !  
  &  na1,      & !  99 if coordinates must be Cartesian, 0 otherwise
  &  lm61,     & !
  &  nbreaks,  & !  Number of breaks in a PDB structure
                 !
  &  no_pKa,   & !  Number of ionizable hydrogen atoms used in pKa calculation
  &  id,       & !  Term        Dimensionality of system
                 !  Definition  Number of infinite dimensions of the system
                 !  Units       none
                 !  Min. value  0  (A discrete molecule or ion)
                 !  Max. value  3  (A regular solid)
                 !
  &  l1u,      & !  Number of unit cells in the first dimension used in solid-state work (0 for a molecule)
  &  l2u,      & !  Number of unit cells in the second dimension (0 for a molecule or polymer)
  &  l3u,      & !  Number of unit cells in the third dimension (0 for a molecule, polymer, or layer system)
  &  l123,     & !  Total number of unit cells (=1 for a molecule = (2*l1u + 1)*(2*l2u + 1)*(2*l3u + 1))
  &  l11,      & !  1 if l1u > 0, 0 otherwise
  &  l21,      & !  1 if l2u > 0, 0 otherwise
  &  l31,      & !  1 if l3u > 0, 0 otherwise
  &  msdel,    & !  Magnetic component of spin
  &  Run = 1,  & !  Run-number (=1 if a single job)
  &  ijulian   & !  Number of days that this program has left before it stops working
  &  = 100000, & !  Be generous
  &  ilim,     & !  Number of temperatures in a thermodynamics calculation
  &  site_no,  & !  The number of this site
  &  ncomments,& !  Number of lines of comment
  &  mers(3),  & !  Number of mers in each direction, in solids
  &  itemp_1,  & !  Used for very temporary transfer of information
  &  itemp_2,  & !  Used for very temporary transfer of information
  &  lpars,    & !  Number of parameters read from parameter files
  &  n_screen, & !  Used to pass information about type of output to "to_screen"
  &  nl_atoms, & !  number of atoms to be printed.
  &  old_chrge,& !  Net charge from the previous calculation, if it exists.
  &  N_Hbonds, & !  Number of hydrogen bonds, with energy < -1 kcal/mol, found in system
  &  P_Hbonds, & !  Number of hydrogen bonds found in system
  &  num_threads, & !  Number of threads to be used (in parallel work)
  &  num_bits, & !  Number of bits used in architecture (32 or 64)
  &  dummy
  real ::          &
    CPU_0,         & !  Start of CPU time (in seconds)
    CPU_1,         & !  Curresnt CPU time (in seconds)
    wall_clock_0,  & !  Start of wall-clock time (in seconds)
    wall_clock_1     !  Current wall-clock time (in seconds)
  real(double), bind(C)::  &
  & arc_hof_1,     & !  Heat of formation from an old ARC file
  & arc_hof_2,     & !  Heat of formation from another old ARC file
  &  escf,         & !  Term        Heat of formation at 298 K
                     !  Units       kcal/mole
                     !  Unit type   Heat of Formation in kcal/mol
                     !  Min. value  ~-20000
                     !  Max. value  ~+20000
                     !
  &  emin,         & !  Lowest heat of formation calculated.
                     !
  &  elect,        & !  Term        Electronic energy 
                     !  Units       eV
                     !  Definition  Total electronic energy of the system in electron volts
                     !  Min. value  Very negative
                     !  Max. value  0.0
                     !
  &  enuclr,       & !  Term        Nuclear energy
                     !  Units       eV
                     !  Definition  Total core-core repulsion energy in eV
                     !  Min. value  0.0
                     !  Max. value  Very large
                     ! 
  &  fract,        & !  Term        Fractional occupancy of M.O.s in active space
                     !  Units       Electrons
                     !  Definition  Number of electrons in active space divided by number of M.O.s in active space
                     !  Min. value  0.0
                     !  Max. value  2.0
                     !  Default     0.0
                     !
  &  gnorm,        & !  Term        Scalar of gradient vector
                     !  Units       kcal/mol/(Angstrom - Radians)
                     !  Definition  The square root of the sum of the squares of the gradient components
                     !
  &  mol_weight,   & !  Term        Molecular weight
                     !  Units       Atomic mass units (Hydrogen = 1.00794 on this scale)
                     !  Min. value  0.0
                     !  Max. value  Very large
                     !
  &  time0 = 0.d0, & !  The start of time, or when the job began (seconds)
  &  tleft,        & !  Number of seconds the job has left before running out of time
  &  tdump = 0.d0, & !  Time between checkpoint dumps (seconds)
  &  cosine,       & !  Angle between previous and current gradient vectors
  &  ux, uy, uz,   & !  Dipole components (?)
  &  step,         & !  "step" in SADDLE calculation
  &  rjkab1,       & !  "J" and "K" integrals for correction of doublet I.P.
  &  atheat,       & !  Sum of atomic heats of formation
  &  hpress,       & !  Energy due to pressure (solids only)
  &  nsp2_corr,    & !  PM6 >N- correction
  &  Si_O_H_corr,  & !  PM7 Si-O-H correction
  &  sum_dihed,    & !  Peptide twist correction (sets barrier in Me-CO-NH-ME to the correct value)
                     !
  &  efield(3),    & !  Term        Electric field
                     !  Definition  Applied external electric field components in X, Y, and Z.
                     !  Units       Atomic units per Bohr
                     !  Min. value  Very negative
                     !  Max. value  Very positive
                     !  Default     0.0, 0.0, 0.0
                     !
  &  cutofp,       & !  Cutoff distance for NDDO approximation in solid state 
                     !  Default     10**(10) for molecules, 30 for polymers, layers, solids
  &  clower,       & !  Lower bound for transition to point charges in solid-state 
                     !  Default     13 Angstroms
  &  cupper,       & !  Upper bound for transition to point charges in solid-state 
                     !  Default     cutofp
                     !
  &  pressure,     & !  Term        Pressure or stress for solid-state and polymer work
                     !  Definition  Applied isotropic pressure (for solids) or pull (for polymers)
                     !  Units       Pascals (Newtons per square meter) (solids) or Newtons per Mole (polymers)
                     !  Default     0.0
  &  zpe,          & !  Term        Zero point energy
                     !  Definition  Half the addition of the vibrational energies
                     !  Units       kcal/mol
  &  density,      & !  Term        Density of solid
                     !  Definition  Weight of cluster divided by its volume
                     !  Units       grams per cubic centimeter
  &  stress,       & !  Term        Energy due to distortion from reference geometry, used with GEO-REF=
                     !  Units       Kcal/mol
                     !
  &  press(3),     & ! Type         Pressure on crystal faces
                     ! Definition   Pressure required to stop crystal expanding
                     ! Units        Gigapascals
  &  E_disp,       & ! Dispersion term
  &  E_hb,         & ! Hydrogen-bond term
  &  E_hh,         & ! Hydrogen-hydrogen repulsion term (to correct for short H - H interactions)
  &  Rab,          & ! Distance between two atom, calculated in "connected"
  &  temp_1,       & !  Used for very temporary transfer of information
  &  temp_2,       & !  Used for very temporary transfer of information
  &  sz,           & ! Spin component
  &  ss2,          & ! Total spin
  param_constant,  & ! Used by PARAM only
   trunc_1,        &
   trunc_2,        &
  &  rdummy
!
   character ::     &
  & ltxt,          & !
  & formula*100,   & !  Type          Empirical formula
                     !  Definition    Type and number count of each element in the system
                     !  Units         Text
 &verson*7= "00.000X"!  Term         Version number
                     !  Definition   Version number for this copy of MOPAC
                     !  Pattern      "\d\.\d\d\d[X|W|L]"
                     !  Description  Year.Julian date. Operating System [X = placeholder]
  character ::     &
  & jobnam*240 = ' ', &!
  & line*1000         !
  character, bind(C) :: keywrd*1000, koment*200, title*200, refkey(6)*1000, geo_dat_name*241, &
    allkey*1000
  character, bind(C) :: errtxt*200, program_name*17="Standalone MOPAC ", dh*20
!
  logical, bind(C) ::            &
     moperr,            & !
     Academic = .true., & !  TRUE is this is an academic licence.
     uhf,               & !  Term        Should the Unrestricted Hartree Fock method be used?
                          !  Definition  True if the calculation is Unrestricted Hartree-Fock
                          !  Default    = .false.
     rhf,               & !  Term        Should the Restricted Hartree Fock method be used?
                          !  Definition  True if the calculation is Restricted Hartree-Fock
                          !  Default    = .false.
     isok,              & !
     mozyme,            & !  TRUE if any of the MOZYME functions are to be used 
     limscf,            & !  Convergence criterion for SCF: if TRUE, then exit the SCF
                          !  if the energy changes a lot (useful in geometry optimization)
                          !  if FALSE, then converge the SCF to the default criterion
     gui = .false.,      & !  By default, output information NOT for a Graphical User Interface
     in_house_only,     & !  TRUE only if run at Stewart Computational Chemistry. Used in debugging, etc.
     lxfac,             & !  TRUE if a diatomic is being used to define the values of XFAC and ALPB
     units,             & !  TRUE if units for input geometry are defined (Angstroms or A0), FALSE otherwise
     Angstroms,         & !  TRUE if units for input geometry must be in Angstroms, if FALSE then A0, see also units
     Sparkle,           & !  TRUE if basis set is missing and sparkles are present (any of elements 58:70)
     keep_res           !  TRUE if the original residue names are to be used
     logical::          &
     MM_corrections(20),& !  &
     N_3_present,       & !
     Si_O_H_present,    &
     dummy_present
     logical, bind(C):: &
     prt_coords,        & ! TRUE if coordinates are to be printed
     prt_gradients,     & ! TRUE if gradients are to be printed
     prt_cart,          & ! TRUE if Cartesian coordinates to be printed
     prt_charges,       & ! TRUE if atomic partial charges to be printed
     prt_pops,          & ! TRUE if atomic orbital populations to be printed
     prt_topo,          & ! TRUE if topography to be printed
     
     is_PARAM=.false.     !  This will be set "TRUE" in a PARAM run
  equivalence  &
    (MM_corrections(1), N_3_present),    & ! TRUE if the system contains at least one N bonded to three ligands
                                           ! and at least two are not hydrogen atoms
    (MM_corrections(2), Si_O_H_present), & ! TRUE if the system contains at least one Si-O-H group
    (MM_corrections(20), dummy_present)    ! If needed, then increase the size of the MM_corrections array.
!
!  Define names for all methods.  Adjust n_methods here and in the equivalence statement lower down.
!
  integer, parameter :: n_methods = 17
  logical::      &
       & methods(n_methods),   &
       & method_mndo,          &   !  1
       & method_am1,           &
       & method_pm3,           &
       & method_rm1,           &
       & method_mndod,         &   !  5
       & method_pm6,           &
       & method_pm6_dh_plus,   &   
       & method_pm6_dh2,       &
       & method_pm6_d3h4,      &
       & method_pm6_dh2x,      &   ! 10
       & method_pm6_d3h4x,     &
       & method_pm6_d3,        &   
       & method_pm6_d3_not_h4, &
       & method_pm7,           &  
       & method_pm7_ts,        &   ! 15
       & method_pm7_hh,        &   
       & method_pm7_minus,     &
       & method_x,             &  ! To be used some day 
       & method             !  Default method = PM6
  character :: methods_keys(n_methods)*11
  data methods_keys/ " MNDO ", " AM1 ", " PM3 ", " RM1 ", " MNDOD ", " PM6 ", "PM6-DH+ ", &
    & " PM6-DH2 ", " PM6-D3H4 ", " PM6-DH2X ", " PM6-D3H4X ", " PM6-D3 ", " PM6-D3(H4)",  &
    & " PM7 ", " PM7-TS ", " PM7-HH ", " PM7- "/
  equivalence (methods(1),  method_mndo), (methods(2),  method_am1), (methods(3),  method_pm3), &
    & (methods(4),  method_rm1), (methods(5),  method_mndod), (methods(6),  method_pm6), &
    & (methods(7),  method_pm6_dh_plus), (methods(8),  method_pm6_dh2), (methods(9),  method_pm6_d3h4), &
    & (methods(10),  method_pm6_dh2x), (methods(11),  method_pm6_d3h4x), (methods(12),method_pm6_d3), &
    & (methods(13),  method_pm6_d3_not_h4), (methods(14), method_pm7), (methods(15),  method_pm7_ts), &
    & (methods(16),  method_pm7_hh), (methods(17),  method_pm7_minus)

    contains
    
    subroutine setMethodsByIdx(i, val) bind(C, name="setMethodsByIdx")
        implicit none
        integer, intent (in) :: i
        logical, intent (in) :: val
        methods(i + 1) = val
      end subroutine
      
      subroutine getMethodsByIdx(i, val) bind(C, name="getMethodsByIdx")
        implicit none
        integer, intent (in) :: i
        logical, intent (out) :: val
        val = methods(i + 1)
      end subroutine
      
    subroutine setMM_correctionsByIdx(i, val) bind(C, name="setMM_correctionsByIdx")
        implicit none
        integer, intent (in) :: i
        logical, intent (in) :: val
        methods(i + 1) = val
      end subroutine
    
end module molkst_C
