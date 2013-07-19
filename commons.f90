module params
  
  integer, parameter :: nvar=8 ! number of variables
  integer :: nx, ny, nz        ! number of MPI process in each direction
  integer :: xdim, ydim, zdim  ! size of the MPI sub-domains
  logical :: inline            ! if 'true', write contiguous data

end module params
