module mod_vars_cuda                                         
  use iso_C_binding
  integer, parameter :: nthreads_gpu = 256, nblocks_gpu = 256
  logical, bind(C) :: lgpu
  real, parameter :: real_cuda = selected_real_kind(8)
  integer, parameter :: prec = 8
  integer :: ngpus,gpu_id
  logical, parameter :: exe_gpu_kepler = .true.      
end module mod_vars_cuda  
      
