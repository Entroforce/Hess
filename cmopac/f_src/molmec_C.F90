      module molmec_C 
      use iso_C_binding
      USE vast_kind_param, ONLY:  double 
!...Created by Pacific-Sierra Research 77to90  4.4G  11:05:01  03/09/06 
      integer, dimension(4,4000), bind(C) :: nhco
      integer, bind(C) :: nnhco      
      real(double), bind(C) :: htype 
      end module molmec_C 
