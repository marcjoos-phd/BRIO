module params
  
  integer, parameter :: nvar=8 ! number of variables
  integer :: nx, ny, nz        ! number of MPI process in each direction
  integer :: xdim, ydim, zdim  ! size of the MPI sub-domains
  integer :: ntoken            ! number of token for I/O
  logical :: inline            ! if 'true', write contiguous data
  logical :: xml               ! if 'true', use XML API for ADIOS I/O

end module params
