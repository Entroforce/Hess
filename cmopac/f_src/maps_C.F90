      module maps_C 
        use iso_C_binding
      double precision, dimension(:), allocatable  :: &
   &  surf,    &
   &  react

      integer :: &
  &  ijlp    , & !  
  &  ilp,      & !  
  &  jlp,      & !  
  &  jlp1,     & !  
  &  ione,     & !
  &  latom,    & !
  &  lparam,   & !
  &  kloop,    &
  &  latom1,   &
  &  lpara1,   &
  &  latom2,   &
  &  lpara2
      double precision, bind(C) :: &
  &  rxn_coord1,  &
  &  rxn_coord2,  &
  &  rxn_coord,   &
  &  rc_escf,     &  !  Reaction coordinate Heat of Formation
  &  ekin,        &  !  Kinetic energy
  &  dummy
      end module maps_C 
