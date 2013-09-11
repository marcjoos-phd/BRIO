program BRIO
  !=============================================================================
  !======================== Benchmark for paRallel I/O =========================
  !=============================================================================
  ! Author: Marc B.R. Joos
  !
  ! Created/last modified: jun 26, 2013/sep 10, 2013
  !
  ! This file is distributed under GNU/GPL licence, 
  ! see <http://www.gnu.org/licenses/>.
  !=============================================================================
  ! This program tests the Parallel I/O libraries in Fortran 90.
  !=============================================================================
  use params
  use writeIO
  use readIO
  use mpi
  implicit none

  integer :: nproc, myrank
  integer :: tag, first, last
  integer :: nparallel, token
  real(8), allocatable, dimension(:,:,:,:) :: data
  integer :: xpos, ypos, zpos
  integer :: i, j, k
  integer :: totalsize
  real(8) :: randval
  integer :: ierr
  real(8) :: tbegin_px, tend_px, tbegin_nc, tend_nc, tbegin_h5, tend_h5
  real(8) :: tbegin_ad, tend_ad, tbegin_mp, tend_mp

  ! MPI initialization
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD, nproc, ierr)

  ! Get parameters
  call read_params

  totalsize = xdim*nx*ydim*ny*zdim*nz*sizeof(randval)

  if(nx*ny*nz .ne. nproc) then
     if(myrank.eq.0) print '(/ "Your domain decomposition does not fit with the number of MPI processes! Exiting BRIO...",/ &
          & "ncpu = ", I10, " different from nx*ny*nz = ", I10 /)', nproc, nx*ny*nz
     call MPI_Finalize(ierr)
     call exit(0)
  endif
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if(myrank.eq.0) then
     print*, "=========================================================="
     print*, "         This is BRIO; Benchmark for paRallel I/O.        "
     print*, "=========================================================="
     print*, " Run parameters: "
     print*, "  - nproc: ", nx*ny*nz
     print*, "  - size of the sub-domains: ", xdim, ydim, zdim
     print*, "  - size of the whole domain: ", xdim*nx, ydim*ny, zdim*nz
#if POSIX == 0 && POTOK == 0
     if(inline)then
        print*, "  - with inline writing"
     else
        print*, "  - with naive writing"
     endif
#endif
#if POTOK == 1
     print*, "  - with ", ntoken, " token"
#endif
#if ADIOS == 1
     if(xml)then
        print*, "  - with XML ADIOS API"
     else
        print*, "  - with no-XML ADIOS API"
     endif
#endif
     print*, "=========================================================="
  endif

  allocate(data(1:xdim,1:ydim,1:zdim,1:nvar))

  ! Position of the thread in a 3D cartesian grid
  xpos = mod(myrank, nx)
  ypos = mod(myrank/(nx*nz), ny)
  zpos = mod(myrank/nx, nz)

  ! Fill data
  do k = 1, zdim
     do j = 1, ydim
        do i = 1, xdim
           call random_number(randval)
           data(i,j,k,1) = randval
           data(i,j,k,2) = sin(i - xdim/2.d0)
           data(i,j,k,3) = sin(j - ydim/2.d0)
           data(i,j,k,4) = sin(k - zdim/2.d0)
           data(i,j,k,5) = (i - xdim/2)
           data(i,j,k,6) = (j - ydim/2)
           data(i,j,k,7) = (k - zdim/2)
           data(i,j,k,8) = 0.
        enddo
     enddo
  enddo

  ! Write & read files
#if POSIX == 1 || POTOK == 1
  call create_dir(myrank)
#endif
#if POSIX == 1
  tbegin_px = MPI_Wtime()
  call write_posix(data, xpos, ypos, zpos, myrank)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  tend_px = MPI_Wtime()
  if(myrank==0) write(*,*) 'Sequential POSIX > writing speed: ' &
       , totalsize/(tend_px - tbegin_px)/1024**2, 'MiB/s'
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  tbegin_px = MPI_Wtime()
  call read_posix(data, xpos, ypos, zpos, myrank)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  tend_px = MPI_Wtime()
  if(myrank==0) write(*,*) 'Sequential POSIX > reading speed: ' &
       , totalsize/(tend_px - tbegin_px)/1024**2, 'MiB/s'
#endif
#if POTOK == 1
  tbegin_px = MPI_Wtime()
  nparallel = nproc/ntoken
  token     = 0
  tag       = 42
  first     = 0
  last      = nproc - 1
  if(myrank >= first .and. myrank <= last) then
     if(mod(myrank, nparallel) .ne. 0) then
        call MPI_Recv(token, 1, MPI_INT, myrank - 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
     endif
     call write_posix(data, xpos, ypos, zpos, myrank)
     if(mod((myrank + 1), nparallel) .ne. 0 .and. myrank < last) then
        call MPI_Send(token, 1, MPI_INT, myrank + 1, tag, MPI_COMM_WORLD, ierr)
     endif
  endif
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  tend_px = MPI_Wtime()
  if(myrank==0) write(*,*) 'Sequential POSIX w/ token > writing speed: ' &
       , totalsize/(tend_px - tbegin_px)/1024**2, 'MiB/s'
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  if(myrank >= first .and. myrank <= last) then
     if(mod(myrank, nparallel) .ne. 0) then
        call MPI_Recv(token, 1, MPI_INT, myrank - 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
     endif
     call write_posix(data, xpos, ypos, zpos, myrank)
     if(mod((myrank + 1), nparallel) .ne. 0 .and. myrank < last) then
        call MPI_Send(token, 1, MPI_INT, myrank + 1, tag, MPI_COMM_WORLD, ierr)
     endif
  endif
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  tend_px = MPI_Wtime()
  if(myrank==0) write(*,*) 'Sequential POSIX w/ token > reading speed: ' &
       , totalsize/(tend_px - tbegin_px)/1024**2, 'MiB/s'
#endif
#if PNCDF == 1
  tbegin_nc = MPI_Wtime()
  call write_pncdf(data, xpos, ypos, zpos, myrank)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  tend_nc = MPI_Wtime()
  if(myrank==0) write(*,*) 'Parallel NetCDF > writing speed: ' &
       , totalsize/(tend_nc - tbegin_nc)/1024**2, 'MiB/s'  
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  tbegin_nc = MPI_Wtime()
  call read_pncdf(data, xpos, ypos, zpos, myrank)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  tend_nc = MPI_Wtime()
  if(myrank==0) write(*,*) 'Parallel NetCDF > reading speed: ' &
       , totalsize/(tend_nc - tbegin_nc)/1024**2, 'MiB/s'  
#endif
#if PHDF5 == 1
  tbegin_h5 = MPI_Wtime()
  call write_phdf5(data, xpos, ypos, zpos, myrank)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  tend_h5 = MPI_Wtime()
  if(myrank==0) write(*,*) 'Parallel HDF5 > writing speed: ' &
       , totalsize/(tend_h5 - tbegin_h5)/1024**2, 'MiB/s'  
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  tbegin_h5 = MPI_Wtime()
  call read_phdf5(data, xpos, ypos, zpos, myrank)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  tend_h5 = MPI_Wtime()
  if(myrank==0) write(*,*) 'Parallel HDF5 > reading speed: ' &
       , totalsize/(tend_h5 - tbegin_h5)/1024**2, 'MiB/s'  
#endif
#if ADIOS == 1
  if(xml) then
     tbegin_ad = MPI_Wtime()
     call write_adios_XML(data, xpos, ypos, zpos, myrank)
     call MPI_Barrier(MPI_COMM_WORLD, ierr)
     tend_ad = MPI_Wtime()
     if(myrank==0) write(*,*) 'ADIOS > writing speed: ' &
       , totalsize/(tend_ad - tbegin_ad)/1024**2, 'MiB/s'  
  else
     tbegin_ad = MPI_Wtime()
     call write_adios_noXML(data, xpos, ypos, zpos, myrank)
     call MPI_Barrier(MPI_COMM_WORLD, ierr)
     tend_ad = MPI_Wtime()
     if(myrank==0) write(*,*) 'ADIOS-noXML > writing speed: ' &
       , totalsize/(tend_ad - tbegin_ad)/1024**2, 'MiB/s'  
  endif

  tbegin_ad = MPI_Wtime()
  call read_adios(data, xpos, ypos, zpos, myrank)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  tend_ad = MPI_Wtime()
  if(myrank==0) write(*,*) 'ADIOS > reading speed: ' &
       , totalsize/(tend_ad - tbegin_ad)/1024**2, 'MiB/s'  
#endif
#if MPIIO == 1
  tbegin_mp = MPI_Wtime()
  call write_mpiio(data, xpos, ypos, zpos, myrank)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  tend_mp = MPI_Wtime()
  if(myrank==0) write(*,*) 'MPI-IO > writing speed: ' &
       , totalsize/(tend_mp - tbegin_mp)/1024**2, 'MiB/s'  
#endif

  deallocate(data)
  call MPI_Finalize(ierr)

  contains
    !===========================================================================
    !=============================  SUBROUTINES  ===============================
    !===========================================================================
    subroutine read_params
      ! Read the parameters for the main program
      ! - interactively if INTERACTIVE is set to 1
      ! - from a namelist 'input' if INTERACTIVE is set to 0
      use params
      implicit none
      
#if INTERACTIVE==1
      integer       :: i,n
      integer       :: iargc
      character(len=7)   :: opt
      character(len=128) :: arg
  
      call default_params

      n = iargc()
      do i = 1,n,2
         call getarg(i,opt)
         if (i == n) then
            print '("option ",a2," has no argument")', opt
            stop 2
         end if
         call getarg(i+1,arg)
         select case (opt)
         case ('-xdim')
            read(arg,*) xdim
         case ('-ydim')
            read(arg,*) ydim
         case ('-zdim')
            read(arg,*) zdim
         case ('-nx')
            read(arg,*) nx
         case ('-ny')
            read(arg,*) ny
         case ('-nz')
            read(arg,*) nz
         case ('-inline')
            read(arg,*) inline
         case ('-ntoken')
            reand(arg,*) ntoken
         case ('-xml')
            read(arg,*) xml
         case default
            print '("unknown option ",a2," ignored")', opt
         end select
      end do
#else
      namelist /model_params/ xdim, ydim, zdim, nx, ny, nz
      namelist /output_params/ inline, ntoken
      namelist /adios_params/ xml

      call default_params
      
      open(unit=1,file='input_BRIO',status='old')
      read(1,model_params)
      read(1,output_params)
      read(1,adios_params)
      close(1)
#endif
      
      return
    end subroutine read_params
    !===========================================================================
    subroutine default_params
      use params
      implicit none

      ! Model parameters
      xdim = 100; nx = 1
      ydim = 100; ny = 1
      zdim = 100; nz = 1

      ! Output parameters
      inline = .false.
      ntoken = 2

      ! ADIOS parameters
      xml = .false.

      return
    end subroutine default_params
    !===========================================================================
    !===========================================================================
    !===========================================================================
end program BRIO
