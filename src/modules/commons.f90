!===============================================================================
!======================== Benchmark for paRallel I/O ===========================
!===============================================================================
! Author: Marc B.R. Joos
!
! Created/last modified: sep 10, 2013/oct 24, 2013
!
! This file is distributed under GNU/GPL license, 
! see <http://www.gnu.org/licenses/>.
!===============================================================================
module params
  
  integer, parameter :: nvar=8 ! number of variables
  integer :: nx, ny, nz        ! number of MPI process in each direction
  integer :: xdim, ydim, zdim  ! size of the MPI sub-domains
  logical :: posix, potok      ! Selection for I/O strategy and library:
  logical :: pncdf, phdf5      ! if 'true', this library is selected
  logical :: adios, mpiio      ! 
  integer :: ntoken            ! number of token for I/O
  logical :: inline            ! if 'true', write contiguous data
  logical :: xml               ! if 'true', use XML API for ADIOS I/O

end module params
!===============================================================================
!===============================================================================
!===============================================================================
module filehandle

contains
!===============================================================================
  subroutine get_filename(myrank, prefix, filename)
    implicit none
    
    integer :: myrank
    character(LEN=*)  :: prefix
    character(LEN=6)  :: nrank
    character(LEN=80) :: filename, filedir
    
    call convtoasc(myrank, nrank)
    filedir = 'sequentialio/'
    filename = trim(filedir)//trim(prefix)//'.'//trim(nrank)
    
    return
  end subroutine get_filename
!===============================================================================
  subroutine convtoasc(number, nb_string)
    ! To convert an integer smaller than 999999 to a 6 characters string
    implicit none
    
    integer :: number, istring, num, nums10, i
    character(LEN=6) :: nb_string
    character(LEN=10), parameter :: nstring="0123456789"
    
    num    = 1000000
    nums10 = num/10
    do i = 1, 6
       istring        = 1 + mod(number,num)/nums10
       nb_string(i:i) = nstring(istring:istring)
       num            = num/10
       nums10         = nums10/10
    enddo
    
    return
  end subroutine convtoasc
!===============================================================================
end module filehandle
