      module parameters_C 
      use iso_C_binding
!
! This module holds all the parameters for the elements.  It is filled
! from the reference parameter sets 
!
      USE vast_kind_param, ONLY:  double 
!...Created by Pacific-Sierra Research 77to90  4.4G  11:04:16  03/09/06  
    logical, dimension(107) :: dorbs
    real(double), dimension(107) :: alp 
    real(double), dimension(107,4) :: guess1, guess2, guess3 
    real(double), dimension(6,107) :: ddp 
    real(double), dimension(9,107) :: po 
    real(double), dimension(107) :: betas, betap, betad 
    real(double), dimension(107) :: uss, upp, udd
    real(double), dimension(107) :: gpp, gp2, hsp, gss, gsp
    real(double), dimension(107) :: am, ad, aq, dd, qq
    real(double), dimension(107), bind(C) :: zs, zp, zd, zsn, zpn, zdn
    real(double), dimension(107), bind(C) :: eisol, eheat
    real(double), dimension(57:71) ::  eheat_sparkles
    real(double), dimension(107), target :: ams
    real(double), dimension (107), bind(C) :: tore, polvol, pocord, f0dd, f2dd, f4dd, &
    f0sd, g2sd, f0pd, f2pd, g1pd, g3pd 
    real(double), dimension (107) :: f0sd_store(107) ! Used by PARAM, not used in MOPAC
    real(double), dimension (107) :: g2sd_store(107) ! Used by PARAM, not used in MOPAC
    real(double), dimension (107) :: atom_radius_vdw, atom_radius_cosmo, &
    &  atom_radius_cosmo_oldcav
    real(double), dimension(100,100) :: xfac, alpb
    double precision :: dh2_a_parameters(6)
    integer, dimension (107,3) :: npq
    integer, dimension (107), bind(C) :: iod, iop, ios, natorb, ndelec
    logical, dimension (107) :: main_group
    double precision :: par1, par2, par3, par4, par5, par6, par7, par8, par9, par10, par11, par12, &
      par13, par14, par15, par16, par17, par18, par19, par20, par21
    character (len=5), dimension (37) :: partyp
    equivalence &
      (guess1(99,1), par1 ), (guess2(99,1), par2 ), (guess3(99,1), par3 ), &
      (guess1(99,2), par4 ), (guess2(99,2), par5 ), (guess3(99,2), par6 ), &
      (guess1(99,3), par7 ), (guess2(99,3), par8 ), (guess3(99,3), par9 ), &
      (guess1(99,4), par10), (guess2(99,4), par11), (guess3(99,4), par12), &
      (guess1(98,1), par13), (guess2(98,1), par14), (guess3(98,1), par15), &
      (guess1(98,2), par16), (guess2(98,2), par17), (guess3(98,2), par18), &
      (guess1(98,3), par19), (guess2(98,3), par20), (guess3(98,3), par21)

    data ndelec / 20 * 0, &
!       Sc Ti V  Cr Mn Fe Co Ni Cu  Zn
      & 0, 0, 2, 2, 4, 4, 6, 8, 10, 10, 8 * 0, &   !  First transition series
      & 0, 0, 2, 2, 4, 4, 6, 8, 10, 10, 22* 0, &   !  Second transition series
      & 0, 0, 2, 2, 4, 4, 6, 8, 10, 10, 27 * 0 /   !  Third transition series
    double precision, dimension (37) :: defmax, defmin  ! Used by PARAM, not used in MOPAC
   
    save
!              H           Initial "s" Orbital Occupancies                     He
!              Li Be                                            B  C  N  O  F  Ne
!              Na Mg                                            Al Si P  S  Cl Ar
!              K  Ca Sc            Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br Kr
!              Rb Sr Y             Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I  Xe
!              Cs Ba La Ce-Lu      Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At Rn
!              Fr Ra Ac Th Pa U    Np Pu Am Cm Bk Cf            Cb ++ +  -- -  Tv
!                                      "s" shell
    data ios &
        &/ 1,                                                                2, &!    2
        &  1, 2,                                              2, 2, 2, 2, 2, 0, &!   10
        &  1, 2,                                              2, 2, 2, 2, 2, 0, &!   18
        &  1, 2, 2,              2, 2, 1, 2, 2, 2, 2, 1, 2,   2, 2, 2, 2, 2, 0, &!   36
        &  1, 2, 2,              2, 1, 1, 2, 1, 1, 0, 1, 2,   2, 2, 2, 2, 2, 0, &!   54
        &  1, 2, 2, 5*0,3*2,6*2, 2, 2, 1, 2, 2, 2, 1, 1, 2,   2, 2, 2, 2, 2, 0, &!   86
        &  1, 1, 2, 4, 2, 2,     2, 2, 2, 2, 2, 1, 0, 3,-3,   1, 2, 1,-2,-1, 0 /
!                                  /
!
!              H           Initial "p" Orbital Occupancies                   He
!              Li Be                                          B  C  N  O  F  Ne
!              Na Mg                                          Al Si P  S  Cl Ar
!              K  Ca Sc          Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br Kr
!              Rb Sr Y           Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I  Xe
!              Cs Ba La Ce-Lu    Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At Rn
!              Fr Ra Ac Th Pa U  Np Pu Am Cm Bk Cf (The rest are reserved for MOPAC)
!                                      "p" shell
    data iop / 0 ,                                                           0, &!    2
            &  0, 0,                                          1, 2, 3, 4, 5, 6, &!   10
            &  0, 0,                                          1, 2, 3, 4, 5, 6, &!   18
            &  0, 0, 0,          0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, &!   36
            &  0, 0, 0,          0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, &!   54
            &  0, 0, 0,  14*0,   0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, &!   86
            &  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9*0                        /
!
!              H           Initial "d" Orbital Occupancies                   He
!              Li Be                                          B  C  N  O  F  Ne
!              Na Mg                                          Al Si P  S  Cl Ar
!              K  Ca Sc          Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br Kr
!              Rb Sr Y           Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I  Xe
!              Cs Ba La Ce-Lu    Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At Rn
!              Fr Ra Ac Th Pa U  Np Pu Am Cm Bk Cf (The rest are reserved for MOPAC)
!                                      "d" shell
    data iod / 0,                                                           0, &!    2
             & 0, 0,                                         0, 0, 0, 0, 0, 0, &!   10
             & 0, 0,                                         0, 0, 0, 0, 0, 0, &!   18
             & 0, 0, 1,          2, 3, 5, 5, 6, 7, 8, 10, 0, 0, 0, 0, 0, 0, 0, &!   36
             & 0, 0, 1,          2, 4, 5, 5, 7, 8,10, 10, 0, 0, 0, 0, 0, 0, 0, &!   54
             & 0, 0, 1,13*0,  1, 2, 3, 5, 5, 6, 7, 9, 10, 0, 0, 0, 0, 0, 0, 0, &!   86
             & 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9*0                     /
!
!                     Principal Quantum Numbers for all shells.
!
!              H                 "s"  shell                                  He
!              Li Be                                          B  C  N  O  F  Ne
!              Na Mg                                          Al Si P  S  Cl Ar
!              K  Ca Sc          Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br Kr
!              Rb Sr Y           Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I  Xe
!              Cs Ba La Ce-Lu    Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At Rn
!              Fr Ra Ac Th-Lr    ?? ?? ?? ??
!
data npq(1:107,1) / &
             & 1,                                                             1, &!  2
             & 2, 2,                                           2, 2, 2, 2, 2, 3, &! 10
             & 3, 3,                                           3, 3, 3, 3, 3, 4, &! 18
             & 4, 4,             4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, &! 36
             & 5, 5,             5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, &! 54
             & 6, 6, 14 * 6,     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, &! 86
             & 11 * 0, 1, 0, 0, 0,3, 5 * 0 /
!
!              H                "p"  shell                                   He
!              Li Be                                          B  C  N  O  F  Ne
!              Na Mg                                          Al Si P  S  Cl Ar
!              K  Ca Sc          Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br Kr
!              Rb Sr Y           Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I  Xe
!              Cs Ba La Ce-Lu    Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At Rn
!              Fr Ra Ac Th-Lr    ?? ?? ?? ??
!
data npq(1:107,2) / &
             & 1,                                                             2, &!  2
             & 2, 2,                                           2, 2, 2, 2, 2, 2, &! 10
             & 3, 3,                                           3, 3, 3, 3, 3, 3, &! 18
             & 4, 4,             4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &! 36
             & 5, 5,             5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &! 54
             & 6, 6, 14 * 6,     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &! 86
             & 21 * 0 /
!
!              H                 "d"  shell                                  He
!              Li Be                                          B  C  N  O  F  Ne
!              Na Mg                                          Al Si P  S  Cl Ar
!              K  Ca Sc          Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br Kr
!              Rb Sr Y           Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I  Xe
!              Cs Ba La Ce-Lu    Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At Rn
!              Fr Ra Ac Th-Lr    ?? ?? ?? ??
!
data npq(1:107,3) / &
             & 0,                                                             0, &!  2
             & 0, 0,                                           0, 0, 0, 0, 0, 0, &! 10
             & 3, 3,                                           3, 3, 3, 3, 3, 4, &! 18
             & 3, 3,             3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, &! 36
             & 4, 4,             4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, &! 54
             & 5, 5, 14 * 5,     5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, &! 86
             & 21 * 0 /
!  Main_group is .true. if the Gss, Gpp, etc., are independent of zsn, zpn, and zdn,
!  .false. otherwise.
!
!  For main-group elements, the Gss value is important, for the transition metals,
!  Gss is less important.
!
    data main_group /&
   &  2*.true.,                              & ! H  - He
   &  8*.true.,                              & ! Li - Ne
   &  8*.true.,                              & ! Na - Ar
   &  2*.true., 9*.false., 7*.true.,         & ! K  - Kr
   &  2*.true., 9*.false., 7*.true.,         & ! Rb - Xe
   &  2*.true.,23*.false., 7*.true.,         & ! Cs - Rn
   & 21*.true.                               / ! Fr - Tv
!
!     ENTHALPIES OF FORMATION OF GASEOUS ATOMS ARE TAKEN FROM \ANNUAL
!     REPORTS,1974,71B,P 117\  THERE ARE SOME SIGNIFICANT DIFFERENCES
!     BETWEEN THE VALUES REPORTED THERE AND THE VALUES PREVIOUSLY IN
!     THE BLOCK DATA OF THIS PROGRAM.  ONLY THE THIRD  ROW ELEMENTS
!     HAVE BEEN UPDATED.
!
! ALL THE OTHER ELEMENTS ARE TAKEN FROM CRC HANDBOOK 1981-1982
    data eheat(1) / 52.102d0 /
    data eheat(2) / 0.000d0 /
    data eheat(3) / 38.410d0 /
    data eheat(4) / 76.960d0 /
    data eheat(5) / 135.700d0 /
    data eheat(6) / 170.890d0 /
    data eheat(7) / 113.000d0 /
    data eheat(8) / 59.559d0 /
    data eheat(9) / 18.890d0 /
    data eheat(10) / 0.000d0 /
    data eheat(11) / 25.650d0 /
    data eheat(12) / 35.000d0 /
    data eheat(13) / 79.490d0 /
    data eheat(14) / 108.390d0 /
    data eheat(15) / 75.570d0 /
    data eheat(16) / 66.400d0 /
    data eheat(17) / 28.990d0 /
    data eheat(18) / 0.000d0 /
    data eheat(19) / 21.420d0 /
    data eheat(20) / 42.600d0 /
    data eheat(21) / 90.300d0 /
    data eheat(22) / 112.300d0 /
    data eheat(23) / 122.900d0 /
    data eheat(24) / 95.000d0 /
    data eheat(25) / 67.700d0 /
    data eheat(26) / 99.300d0 /
    data eheat(27) / 102.400d0 /
    data eheat(28) / 102.800d0 /
    data eheat(29) / 80.700d0 /
    data eheat(30) / 31.170d0 /
    data eheat(31) / 65.400d0 /
    data eheat(32) / 89.500d0 /
    data eheat(33) / 72.300d0 /
    data eheat(34) / 54.300d0 /
    data eheat(35) / 26.740d0 /
    data eheat(36) / 0.000d0 /
    data eheat(37) / 19.600d0 /
    data eheat(38) / 39.100d0 /
    data eheat(39) / 101.500d0 /
    data eheat(40) / 145.500d0 /
    data eheat(41) / 172.400d0 /
    data eheat(42) / 157.300d0 /
    data eheat(43) / 162.000d0 /
    data eheat(44) / 155.500d0 /
    data eheat(45) / 133.000d0 /
    data eheat(46) / 90.000d0 /
    data eheat(47) / 68.100d0 /
    data eheat(48) / 26.720d0 /
    data eheat(49) / 58.000d0 /
    data eheat(50) / 72.200d0 /
    data eheat(51) / 63.200d0 /
    data eheat(52) / 47.000d0 /
    data eheat(53) / 25.517d0 /
    data eheat(54) / 0.000d0 /
    data eheat(55) / 18.700d0 /
    data eheat(56) / 42.500d0 /
    data eheat(57) / 103.1D0/   
    data eheat(58) / 101.3D0 /  
  !  data eheat(59) / 952.9D0 /  
  !  data eheat(60) / 962.8D0 /  
  !  data eheat(61) / 976.9D0 /  
    data eheat(62) / 49.4d0 /  
  !  data eheat(63) / 1006.6d0/  
  !  data eheat(64) / 991.37d0/  
  !  data eheat(65) / 999.0d0/   
  !  data eheat(66) / 1001.3d0 / 
  !  data eheat(67) / 1009.6d0 / 
    data eheat(68) / 75.8d0 /
  !  data eheat(69) / 1022.06d0/ 
    data eheat(70) / 36.36d0/ 
    data eheat(71) / 102.1d0 / 
    data eheat(72) / 148.000d0 / 
    data eheat(73) / 186.900d0 / 
    data eheat(74) / 203.100d0 /
    data eheat(75) / 185.000d0 /
    data eheat(76) / 188.000d0 /
    data eheat(77) / 160.000d0 /
    data eheat(78) / 135.200d0 /
    data eheat(79) / 88.000d0 /
    data eheat(80) / 14.690d0 /
    data eheat(81) / 43.550d0 /
    data eheat(82) / 46.620d0 /
    data eheat(83) / 50.100d0 /
    data eheat(86) / 0.000d0 /
    data eheat(90) / 1674.64d0/
    data eheat(102) / 207.0d0 /
    data eheat(103) / 0.0d0 /
    data eheat(104) / 0.0d0 /
    data eheat(105) / 0.0d0 /
    data eheat(106) / 0.0d0 /
    data eheat(107) / 0.0d0 /
    data eheat_sparkles(57) / 928.9D0/   !  Represents La(+++)
    data eheat_sparkles(58) / 944.7D0 /  !  Represents Ce(+++)
    data eheat_sparkles(59) / 952.9D0 /  !  Represents Pr(+++)
    data eheat_sparkles(60) / 962.8D0 /  !  Represents Nd(+++)
    data eheat_sparkles(61) / 976.9D0 /  !  Represents Pm(+++)
    data eheat_sparkles(62) / 974.4d0 /  !  Represents Sm(+++)
    data eheat_sparkles(63) / 1006.6d0/  !  Represents Eu(+++)
    data eheat_sparkles(64) / 991.37d0/  !  Represents Gd(+++)
    data eheat_sparkles(65) / 999.0d0/   !  Represents Tb(+++)
    data eheat_sparkles(66) / 1001.3d0 / !  Represents Dy(+++)
    data eheat_sparkles(67) / 1009.6d0 / !  Represents Ho(+++)
    data eheat_sparkles(68) / 1016.15d0 /!  Represents Er(+++)
    data eheat_sparkles(69) / 1022.06d0/ !  Represents Tm(+++)
    data eheat_sparkles(70) / 1039.03d0/ !  Represents Ho(+++)
    data eheat_sparkles(71) / 1031.2d0 / !  Represents Lu(+++)
    data eisol / 107 * 0.0d0 /
    data ams / 1.00790d0, 4.00260d0, 6.94000d0, 9.01218d0, 10.81000d0, &
   & 12.01100d0, 14.00670d0, 15.99940d0, 18.99840d0, 20.17900d0, 22.98977d0, &
   & 24.30500d0, 26.98154d0, 28.08550d0, 30.97376d0, 32.06000d0, 35.45300d0, &
   & 39.94800d0, 39.09830d0, 40.08000d0, 44.95590d0, 47.90000d0, 50.94150d0, &
   & 51.99600d0, 54.93800d0, 55.84700d0, 58.93320d0, 58.71000d0, 63.54600d0, &
   & 65.38000d0, 69.73500d0, 72.59000d0, 74.92160d0, 78.96000d0, 79.90400d0, &
   & 83.80000d0, 85.46780d0, 87.62000d0, 88.90590d0, 91.22000d0, 92.90640d0, &
   & 95.94000d0, 98.90620d0, 101.0700d0, 102.9055d0, 106.4000d0, 107.8680d0, &
   & 112.4100d0, 114.8200d0, 118.6900d0, 121.7500d0, 127.6000d0, 126.9045d0, &
   & 131.3000d0, 132.9054d0, 137.3300d0, 138.9060d0, 140.1160d0, 140.9077d0, &
   & 144.2400d0, 145.0000d0, 150.3600d0, 151.9640d0, 157.2500d0, 158.9253d0, &
   & 162.5000d0, 164.9303d0, 167.2600d0, 168.9342d0, 173.0400d0, 174.9670d0, &
   & 178.4900d0, 180.9479d0, &
   & 183.8500d0, 186.2070d0, 190.2000d0, 192.2200d0, 195.0900d0, 196.9665d0, &
   & 200.5900d0, 204.3700d0, 207.2000d0, 208.9804d0, 209.0000d0, 210.0000d0, &
   & 222.0000d0, 223.0000d0, 226.0000d0, 227.0000d0, 232.0381d0, 231.0359d0, &
   & 238.0289d0, 5 * 0.000d0, 0.00005d0, 3*0.d0, 1.0079d0, &
   & 5 * 0.000d0 /
    data partyp (1),  defmin (1),  defmax (1) / "USS  ", - 200.d0, 60.d0 /
    data partyp (2),  defmin (2),  defmax (2) / "UPP  ", - 200.d0, 60.d0 /
    data partyp (3),  defmin (3),  defmax (3) / "UDD  ", - 200.d0, 60.d0 /
    data partyp (4),  defmin (4),  defmax (4) / "ZS   ", 0.6d0, 6.d0 /
    data partyp (5),  defmin (5),  defmax (5) / "ZP   ", 0.6d0, 6.d0 /
    data partyp (6),  defmin (6),  defmax (6) / "ZD   ", 0.6d0, 6.d0 /
    data partyp (7),  defmin (7),  defmax (7) / "BETAS", - 70.d0, 10.d0 /
    data partyp (8),  defmin (8),  defmax (8) / "BETAP", - 70.d0, 10.d0 /
    data partyp (9),  defmin (9),  defmax (9) / "BETAD", - 70.d0, 10.d0 /
    data partyp (10), defmin (10), defmax (10) / "GSS  ", 1.d0, 20.d0 /
    data partyp (11), defmin (11), defmax (11) / "GSP  ", 1.d0, 20.d0 /
    data partyp (12), defmin (12), defmax (12) / "GPP  ", 1.d0, 20.d0 /
    data partyp (13), defmin (13), defmax (13) / "GP2  ", 1.0d0, 19.d0 /
    data partyp (14), defmin (14), defmax (14) / "HSP  ", 1.d0, 5.d0 /
    data partyp (15), defmin (15), defmax (15) / "F0SD ", 1.d0,10.d0 /
    data partyp (16), defmin (16), defmax (16) / "G2SD ", 1.d0,10.d0 /
    data partyp (17), defmin (17), defmax (17) / "POC  ", 1.d0,10.d0 /
    data partyp (18), defmin (18), defmax (18) / "ALP  ", 0.5d0, 6.d0 /
    data partyp (19), defmin (19), defmax (19) / "ZSN  ", 0.5d0,25.d0 /
    data partyp (20), defmin (20), defmax (20) / "ZPN  ", 0.5d0,25.d0 /
    data partyp (21), defmin (21), defmax (21) / "ZDN  ", 0.5d0,25.d0 /
    data partyp (22), defmin (22), defmax (22) / "FN11 ", -1.d0, 1.d0 /
    data partyp (23), defmin (23), defmax (23) / "FN21 ", 0.5d0, 3.d0 /
    data partyp (24), defmin (24), defmax (24) / "FN31 ", 0.5d0, 3.d0 /
    data partyp (25), defmin (25), defmax (25) / "FN12 ", - 1.d0, 1.d0 /
    data partyp (26), defmin (26), defmax (26) / "FN22 ", 0.15d0, 3.d0 /
    data partyp (27), defmin (27), defmax (27) / "FN32 ", 0.5d0, 3.5d0 /
    data partyp (28), defmin (28), defmax (28) / "FN13 ", - 1.d0, 1.d0 /
    data partyp (29), defmin (29), defmax (29) / "FN23 ", 1.d0, 3.d0 /
    data partyp (30), defmin (30), defmax (30) / "FN33 ", 0.5d0, 4.d0 /
    data partyp (31), defmin (31), defmax (31) / "FN14 ", - 1.d0, 1.d0 /
    data partyp (32), defmin (32), defmax (32) / "FN24 ", 1.d0, 3.d0 /
    data partyp (33), defmin (33), defmax (33) / "FN34 ", 0.5d0, 3.d0 /
    data partyp (34), defmin (34), defmax (34) / "ALPB_", 0.9d0, 3.0d0 /
    data partyp (35), defmin (35), defmax (35) / "XFAC_", 0.5d0, 30.d0 /
    data partyp (36), defmin (36), defmax (36) / "FN33 ", 0.5d0, 6.d0 /
    data partyp (37), defmin (37), defmax (37) / "NORBS", 0.9d0, 9.1d0/
        
    contains
    
    subroutine getAms(ptr) bind(C, name="getAms")
        type(c_ptr), intent(out) :: ptr

        ptr = c_loc(ams(1))
    end subroutine
end module parameters_C 
