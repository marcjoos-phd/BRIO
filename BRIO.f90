program BRIO
  !=============================================================================
  !======================== Benchmark for paRallel I/O =========================
  !=============================================================================
  ! Author: Marc B.R. Joos
  !
  ! Created/last modified: jun 26, 2013/jul 22, 2013
  !
  ! This file is distributed under GNU/GPL licence, 
  ! see <http://www.gnu.org/licenses/>.
  !=============================================================================
  ! This program tests the Parallel I/O libraries in Fortran 90.
  !=============================================================================
  use params
  implicit none

  include "mpif.h"

  integer :: nproc, myrank
  real(8), allocatable, dimension(:,:,:,:) :: data
  integer :: xpos, ypos, zpos
  integer :: i, j, k
  real(8) :: randval
  integer :: ierr
  real(8) :: tbegin_px, tend_px, tbegin_nc, tend_nc, tbegin_h5, tend_h5
  real(8) :: tbegin_ad, tend_ad

  ! MPI initialization
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD, nproc, ierr)

  ! Get parameters
  call read_params
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
     if(inline)then
        print*, "  - with inline writing"
     else
        print*, "  - with naive writing"
     endif
     if(xml)then
        print*, "  - with XML ADIOS API"
     else
        print*, "  - with no-XML ADIOS API"
     endif
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

  ! Write files
#if POSIX == 1
  tbegin_px = MPI_Wtime()
  call write_posix(data, xpos, ypos, zpos, myrank)
  tend_px = MPI_Wtime()
  if(myrank==0) write(*,*) 'Sequential POSIX > time elapsed: '&
       , tend_px - tbegin_px, 's'  
#endif
#if PNCDF == 1
  tbegin_nc = MPI_Wtime()
  call write_pncdf(data, xpos, ypos, zpos, myrank)
  tend_nc = MPI_Wtime()
  if(myrank==0) write(*,*) 'Parallel NetCDF > time elapsed: '&
       , tend_nc - tbegin_nc, 's'  
#endif
#if PHDF5 == 1
  tbegin_h5 = MPI_Wtime()
  call write_phdf5(data, xpos, ypos, zpos, myrank)
  tend_h5 = MPI_Wtime()
  if(myrank==0) write(*,*) 'Parallel HDF5 > time elapsed: '&
       , tend_h5 - tbegin_h5, 's'  
#endif
#if ADIOS == 1
  if(xml) then
     tbegin_ad = MPI_Wtime()
     call write_adios_XML(data, xpos, ypos, zpos, myrank)
     tend_ad = MPI_Wtime()
     if(myrank==0) write(*,*) 'ADIOS > time elapsed: '&
          , tend_ad - tbegin_ad, 's'  
  else
     tbegin_ad = MPI_Wtime()
     call write_adios_noXML(data, xpos, ypos, zpos, myrank)
     tend_ad = MPI_Wtime()
     if(myrank==0) write(*,*) 'ADIOS-noXML > time elapsed: '&
          , tend_ad - tbegin_ad, 's'  
  endif
#endif

  deallocate(data)
  call MPI_Finalize(ierr)

  contains
    !===========================================================================
    !=============================  SUBROUTINES  ===============================
    !===========================================================================
    subroutine write_posix(data, xpos, ypos, zpos, myrank)
      use params
      implicit none

      integer :: xpos, ypos, zpos, myrank
      real(8), dimension(xdim,ydim,zdim,nvar) :: data
      integer, dimension(3) :: boxsize, domdecomp
      character(LEN=80) :: filename
      
      call get_filename(myrank, 'posix', filename)

      ! Define and write some metadata
      boxsize   = (/ xdim, ydim, zdim /)
      domdecomp = (/ nx, ny, nz /)
      
      open(unit=10, file=filename, status='unknown', form='unformatted')
      write(10) boxsize
      write(10) domdecomp
      write(10) data
      close(10)

      return
    end subroutine write_posix
    !===========================================================================
    subroutine write_pncdf(data, xpos, ypos, zpos, myrank)
      use pnetcdf
      use params
      implicit none
      
      include "mpif.h"
      
      integer :: xpos, ypos, zpos, myrank
      real(8), dimension(xdim,ydim,zdim,nvar) :: data
      integer, dimension(3) :: boxsize, domdecomp
      character(LEN=13) :: filename

      ! PnetCDF variables
      integer(kind=MPI_OFFSET_KIND) :: nxtot, nytot, nztot
      integer :: nout, ncid, bdid, ddid, xdimid, ydimid, zdimid
      integer, dimension(3) :: sdimid
      integer :: bsid
      integer :: vid1, vid2, vid3, vid4, vid5, vid6, vid7, vid8
      integer(kind=MPI_OFFSET_KIND), dimension(3) :: dims, start, count
      integer :: diminline
      integer :: ierr
      
      ! Define the global dimension and the local dimension
      if(inline) then
         nxtot = xdim
         nytot = ydim
         nztot = nx*ny*nz*zdim
      else
         nxtot = nx*xdim
         nytot = ny*ydim
         nztot = nz*zdim
      endif

      dims = (/ xdim, ydim, zdim /)
      
      ! Create file & write metadata (by one thread only)
      filename = 'parallelio.nc'
      if(myrank.eq.0) then
         nout = nfmpi_create(MPI_COMM_SELF, filename, NF_CLOBBER, MPI_INFO_NULL&
              , ncid)
      
         ! Define and write some metadata
         boxsize   = (/ xdim, ydim, zdim /)
         domdecomp = (/ nx, ny, nz /)
         
         nout = nfmpi_def_dim(ncid, "boxdim", 3_8, bdid)
         nout = nfmpi_def_var(ncid, "boxsize", NF_INT, 1, (/bdid/), bsid)
         nout = nfmpi_def_var(ncid, "domdecomp", NF_INT, 1, (/bdid/), ddid)
         nout = nfmpi_enddef(ncid)
         
         ! Write the metadata & close the file
         nout = nfmpi_put_vara_int_all(ncid, bsid, (/1_8/), (/3_8/), boxsize)
         nout = nfmpi_put_vara_int_all(ncid, ddid, (/1_8/), (/3_8/), domdecomp)
         nout = nfmpi_close(ncid)
      end if
      
      ! Synchronize the threads
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      
      ! Reopen the file and allow the definition of new variables
      nout = nfmpi_open(MPI_COMM_WORLD, filename, NF_WRITE, MPI_INFO_NULL, ncid)
      nout = nfmpi_redef(ncid)
      
      ! Define dimensions
      nout = nfmpi_def_dim(ncid, "x", nxtot, xdimid)
      nout = nfmpi_def_dim(ncid, "y", nytot, ydimid)
      nout = nfmpi_def_dim(ncid, "z", nztot, zdimid)
      sdimid = (/ xdimid, ydimid, zdimid /)
      
      ! Create variable
      nout = nfmpi_def_var(ncid, "var1", NF_DOUBLE, 3, sdimid, vid1)
      nout = nfmpi_def_var(ncid, "var2", NF_DOUBLE, 3, sdimid, vid2)
      nout = nfmpi_def_var(ncid, "var3", NF_DOUBLE, 3, sdimid, vid3)
      nout = nfmpi_def_var(ncid, "var4", NF_DOUBLE, 3, sdimid, vid4)
      nout = nfmpi_def_var(ncid, "var5", NF_DOUBLE, 3, sdimid, vid5)
      nout = nfmpi_def_var(ncid, "var6", NF_DOUBLE, 3, sdimid, vid6)
      nout = nfmpi_def_var(ncid, "var7", NF_DOUBLE, 3, sdimid, vid7)
      nout = nfmpi_def_var(ncid, "var8", NF_DOUBLE, 3, sdimid, vid8)
      
      ! End of definitions
      nout = nfmpi_enddef(ncid)
      
      ! Indices start at 1!
      if(inline) then
         diminline = myrank*dims(3) + 1 
         start = (/ 1, 1, diminline /)
      else
         start = (/ xpos, ypos, zpos /)*dims+1
      endif
      count = dims
      
      ! Write data
      nout = nfmpi_put_vara_double_all(ncid, vid1, start, count, data(:,:,:,1))
      nout = nfmpi_put_vara_double_all(ncid, vid2, start, count, data(:,:,:,2))
      nout = nfmpi_put_vara_double_all(ncid, vid3, start, count, data(:,:,:,3))
      nout = nfmpi_put_vara_double_all(ncid, vid4, start, count, data(:,:,:,4))
      nout = nfmpi_put_vara_double_all(ncid, vid5, start, count, data(:,:,:,5))
      nout = nfmpi_put_vara_double_all(ncid, vid6, start, count, data(:,:,:,6))
      nout = nfmpi_put_vara_double_all(ncid, vid7, start, count, data(:,:,:,7))
      nout = nfmpi_put_vara_double_all(ncid, vid8, start, count, data(:,:,:,8))
      
      ! Close file
      nout = nfmpi_close(ncid)

      return
    end subroutine write_pncdf
    !===========================================================================
    subroutine write_phdf5(data, xpos, ypos, zpos, myrank)
      use hdf5
      use params
      implicit none
      
      include "mpif.h"
      
      integer :: xpos, ypos, zpos, myrank
      real(8), dimension(xdim,ydim,zdim,nvar) :: data
      integer, dimension(3) :: boxsize, domdecomp
      character(LEN=13) :: filename

      ! HDF5 variables
      integer :: ierr
      integer(HID_T) :: file_id, fapl_id, driver_id
      integer(HID_T) :: h5_dspace, h5_dset
      integer(HSIZE_T), dimension(3) :: dim_meta

      ! Initialize HDF5 interface
      call H5open_f(ierr)

      ! Create HDF5 property IDs for parallel file access
      filename = 'parallelio.h5'
      call H5Pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
      call H5Pset_fapl_mpio_f(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)
      call H5Pget_driver_f(fapl_id, driver_id, ierr)
      call H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr &
           , access_prp=fapl_id)

      ! Define and write some metadata
      boxsize   = (/ xdim, ydim, zdim /)
      domdecomp = (/ nx, ny, nz /)
      dim_meta  = 3
      call H5Screate_simple_f(1, dim_meta, h5_dspace, ierr)
      call H5Dcreate_f(file_id, "boxsize", H5T_NATIVE_INTEGER, &
           & h5_dspace, h5_dset, ierr)
      call H5Dwrite_f(h5_dset, H5T_NATIVE_INTEGER, boxsize, dim_meta, ierr)
      call H5Dclose_f(h5_dset, ierr)
      call H5Sclose_f(h5_dspace, ierr)
      call H5Screate_simple_f(1, dim_meta, h5_dspace, ierr)
      call H5Dcreate_f(file_id, "domdecomp", H5T_NATIVE_INTEGER, &
           & h5_dspace, h5_dset, ierr)
      call H5Dwrite_f(h5_dset, H5T_NATIVE_INTEGER, domdecomp, dim_meta, ierr)
      call H5Dclose_f(h5_dset, ierr)
      call H5Sclose_f(h5_dspace, ierr)

      ! Write data
      call write_dataset_h5(file_id, "var1", data(:,:,:,1), xpos, ypos, zpos&
           , myrank)
      call write_dataset_h5(file_id, "var2", data(:,:,:,2), xpos, ypos, zpos&
           , myrank)
      call write_dataset_h5(file_id, "var3", data(:,:,:,3), xpos, ypos, zpos&
           , myrank)
      call write_dataset_h5(file_id, "var4", data(:,:,:,4), xpos, ypos, zpos&
           , myrank)
      call write_dataset_h5(file_id, "var5", data(:,:,:,5), xpos, ypos, zpos&
           , myrank)
      call write_dataset_h5(file_id, "var6", data(:,:,:,6), xpos, ypos, zpos&
           , myrank)
      call write_dataset_h5(file_id, "var7", data(:,:,:,7), xpos, ypos, zpos&
           , myrank)
      call write_dataset_h5(file_id, "var8", data(:,:,:,8), xpos, ypos, zpos&
           , myrank)

      ! Close HDF5 property IDs and interface
      call H5Fclose_f(file_id, ierr)
      call H5Pclose_f(fapl_id, ierr)
      call H5close_f(ierr)

      return
    end subroutine write_phdf5
    !===========================================================================
    subroutine write_dataset_h5(loc_id, dsetname, data, xpos, ypos, zpos&
         , myrank)
      use hdf5
      use params
      implicit none

      character(LEN=*) :: dsetname
      real(8), dimension(xdim,ydim,zdim) :: data
      integer :: xpos, ypos, zpos, myrank

      ! HDF5 variables
      integer(HID_T) :: prop_data_chunk_id, dxpl_id
      integer(HID_T) :: loc_id, h5_dspace, h5_dset, h5_dspace_file
      integer(HSIZE_T), dimension(3) :: start, count, stride, blockSize
      integer(HSIZE_T), dimension(3) :: dims, dims_file, dims_chunk

      dims      = (/ xdim, ydim, zdim /)
      if(inline) then
         dims_file = (/ xdim, ydim, zdim*nx*ny*nz /)
      else
         dims_file = (/ xdim*nx, ydim*ny, zdim*nz /)
      endif

      call H5Screate_simple_f(3, dims, h5_dspace, ierr)
      call H5Screate_simple_f(3, dims_file, h5_dspace_file, ierr)
       
      ! Hyperslab for selecting data in h5_dspace
      start     = (/ 0, 0, 0 /)
      stride    = (/ 1, 1, 1 /)
      count     = dims
      blockSize = (/ 1, 1, 1 /)
      call H5Sselect_hyperslab_f(h5_dspace, H5S_SELECT_SET_F, start, count&
           , ierr, stride, blockSize)
       
      ! Hyperslab for selecting location in h5_dspace_file (to set the
      ! correct location in file where we want to put our piece of data)
      if(inline) then
         start(1) = 0
         start(2) = 0
         start(3) = dims(3)*myrank
      else
         start     = (/ xpos, ypos, zpos /)*dims
      endif
      stride    = (/ 1,1,1 /)
      count     = dims
      blockSize = (/ 1,1,1 /)
      call H5Sselect_hyperslab_f(h5_dspace_file, H5S_SELECT_SET_F, start, count&
           , ierr, stride, blockSize)
       
      ! create property id for our data chunk
      dims_chunk = dims
      call H5Pcreate_f(H5P_DATASET_CREATE_F, prop_data_chunk_id, ierr)
      call H5Pset_chunk_f(prop_data_chunk_id, 3, dims_chunk, ierr)
       
      ! enable parallel collective IO
      call H5Pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
      call H5Pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)
       
      ! create data set
      call H5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_DOUBLE&
           , h5_dspace_file, h5_dset, ierr, prop_data_chunk_id) 
       
      ! finally write data to file
      call H5Dwrite_f(h5_dset, H5T_NATIVE_DOUBLE, data, dims, ierr&
           , mem_space_id=h5_dspace, file_space_id=h5_dspace_file&
           , xfer_prp=dxpl_id)
       
      ! clean HDF5 IDs
      call H5Pclose_f(prop_data_chunk_id, ierr)
      call H5Pclose_f(dxpl_id, ierr)
      call H5Dclose_f(h5_dset, ierr)
      call H5Sclose_f(h5_dspace, ierr)
      call H5Sclose_f(h5_dspace_file, ierr)

      return
    end subroutine write_dataset_h5
    !===========================================================================
    subroutine write_adios_noXML(data, xpos, ypos, zpos, myrank)
      use adios_write_mod
      use params
      implicit none
      
      include "mpif.h"
      
      integer :: xpos, ypos, zpos, myrank
      real(8), dimension(xdim,ydim,zdim,nvar) :: data
      integer, dimension(3) :: boxsize, domdecomp
      character(LEN=13) :: filename

      ! MPI & ADIOS variables
      integer :: sizeMB=64
      integer(8) :: adios_group, adios_handle, varid
      integer(8) :: group_size, total_size, offset_x, offset_y, offset_z
      integer :: xdimglob, ydimglob, zdimglob
      integer :: ierr

      ! Init ADIOS
      call ADIOS_Init_noXML(MPI_COMM_WORLD, ierr)
      call ADIOS_Allocate_Buffer(sizeMB, ierr)

      ! Group declaration:
      !  - first argument is name of the group
      !  - second argument is time index, if we want to store some variables
      ! over time
      !  - third argument is a flag saying if we want statistics (with some
      ! overhead)
      filename = "parallelio.bp"
      call ADIOS_Declare_Group(adios_group, "dump", "", 0, ierr)

      ! Method selection:
      !  - first argument in 'POSIX', 'MPI', 'PHDF5' etc.
      !  - second argument is options relatives to the chosen method
      !  - third argument is root directory; default is current directory
      call ADIOS_Select_Method(adios_group, "MPI", "", "", ierr)

      ! Metadata definitions
      boxsize   = (/ xdim, ydim, zdim /)
      domdecomp = (/ nx, ny, nz /)
      call ADIOS_Define_Var(adios_group, "xdim", "", adios_integer, "", "" &
           & , "", varid)
      call ADIOS_Define_Var(adios_group, "ydim", "", adios_integer, "", "" &
           & , "", varid)
      call ADIOS_Define_Var(adios_group, "zdim", "", adios_integer, "", "" &
           & , "", varid)
      call ADIOS_Define_Var(adios_group, "xdimglob", "", adios_integer, "", "" &
           & , "", varid)
      call ADIOS_Define_Var(adios_group, "ydimglob", "", adios_integer, "", "" &
           & , "", varid)
      call ADIOS_Define_Var(adios_group, "zdimglob", "", adios_integer, "", "" &
           & , "", varid)
      call ADIOS_Define_Var(adios_group, "offset_x", "", adios_integer, "" &
           & , "", "", varid)
      call ADIOS_Define_Var(adios_group, "offset_y", "", adios_integer, "" &
           & , "" , "", varid)
      call ADIOS_Define_Var(adios_group, "offset_z", "", adios_integer, "" &
           & , "" , "", varid)
      call ADIOS_Define_Var(adios_group, "boxsize", "", adios_integer &
           & , "3", "", "", varid)
      call ADIOS_Define_Var(adios_group, "domdecomp", "", adios_integer &
           & , "3", "", "", varid)

      ! Data definitions
      call ADIOS_Define_Var(adios_group, "var1", "", adios_double &
        & , "xdim, ydim, zdim", "xdimglob, ydimglob, zdimglob" &
        & , "offset_x, offset_y, offset_z", varid)
      call ADIOS_Define_Var(adios_group, "var2", "", adios_double &
        & , "xdim, ydim, zdim", "xdimglob, ydimglob, zdimglob" &
        & , "offset_x, offset_y, offset_z", varid)
      call ADIOS_Define_Var(adios_group, "var3", "", adios_double &
        & , "xdim, ydim, zdim", "xdimglob, ydimglob, zdimglob" &
        & , "offset_x, offset_y, offset_z", varid)
      call ADIOS_Define_Var(adios_group, "var4", "", adios_double &
        & , "xdim, ydim, zdim", "xdimglob, ydimglob, zdimglob" &
        & , "offset_x, offset_y, offset_z", varid)
      call ADIOS_Define_Var(adios_group, "var5", "", adios_double &
        & , "xdim, ydim, zdim", "xdimglob, ydimglob, zdimglob" &
        & , "offset_x, offset_y, offset_z", varid)
      call ADIOS_Define_Var(adios_group, "var6", "", adios_double &
        & , "xdim, ydim, zdim", "xdimglob, ydimglob, zdimglob" &
        & , "offset_x, offset_y, offset_z", varid)
      call ADIOS_Define_Var(adios_group, "var7", "", adios_double &
        & , "xdim, ydim, zdim", "xdimglob, ydimglob, zdimglob" &
        & , "offset_x, offset_y, offset_z", varid)
      call ADIOS_Define_Var(adios_group, "var8", "", adios_double &
        & , "xdim, ydim, zdim", "xdimglob, ydimglob, zdimglob" &
        & , "offset_x, offset_y, offset_z", varid)

      ! Open file
      call ADIOS_Open(adios_handle, "dump", filename, "w", MPI_COMM_WORLD, ierr)

      ! Size of the group = 4*(# metadata) + 8*(data size)
      group_size = 4*20 + 8*xdim*ydim*zdim*8
      call ADIOS_Group_Size(adios_handle, group_size, total_size, ierr)

      ! Define offset and global dimensions
      if(inline) then
         offset_x = 0
         offset_y = 0
         offset_z = zdim*myrank
         xdimglob = xdim; ydimglob = ydim; zdimglob = zdim*nx*ny*nz
      else
         offset_x = xdim*xpos
         offset_y = ydim*ypos
         offset_z = zdim*zpos
         xdimglob = xdim*nx; ydimglob = ydim*ny; zdimglob = zdim*nz
      endif

      ! Write metadata
      call ADIOS_Write(adios_handle, "xdim", xdim, ierr)
      call ADIOS_Write(adios_handle, "ydim", ydim, ierr)
      call ADIOS_Write(adios_handle, "zdim", zdim, ierr)
      call ADIOS_Write(adios_handle, "xdimglob", xdimglob, ierr)
      call ADIOS_Write(adios_handle, "ydimglob", ydimglob, ierr)
      call ADIOS_Write(adios_handle, "zdimglob", zdimglob, ierr)
      call ADIOS_Write(adios_handle, "offset_x", offset_x, ierr)
      call ADIOS_Write(adios_handle, "offset_y", offset_y, ierr)
      call ADIOS_Write(adios_handle, "offset_z", offset_z, ierr)
      call ADIOS_Write(adios_handle, "boxsize", boxsize, ierr)
      call ADIOS_Write(adios_handle, "domdecomp", domdecomp, ierr)

      ! Write data
      call ADIOS_Write(adios_handle, "var1", data(:,:,:,1), ierr)
      call ADIOS_Write(adios_handle, "var2", data(:,:,:,2), ierr)
      call ADIOS_Write(adios_handle, "var3", data(:,:,:,3), ierr)
      call ADIOS_Write(adios_handle, "var4", data(:,:,:,4), ierr)
      call ADIOS_Write(adios_handle, "var5", data(:,:,:,5), ierr)
      call ADIOS_Write(adios_handle, "var6", data(:,:,:,6), ierr)
      call ADIOS_Write(adios_handle, "var7", data(:,:,:,7), ierr)
      call ADIOS_Write(adios_handle, "var8", data(:,:,:,8), ierr)

      ! Close file and ADIOS interface
      call ADIOS_Close(adios_handle, ierr)
      call ADIOS_Finalize(myrank, ierr)

      return
    end subroutine write_adios_noXML
    !===========================================================================
    subroutine write_adios_XML(data, xpos, ypos, zpos, myrank)
      use adios_write_mod
      use params
      implicit none
      
      include "mpif.h"
      
      integer :: xpos, ypos, zpos, myrank
      real(8), dimension(xdim,ydim,zdim,nvar) :: data
      integer, dimension(3) :: boxsize, domdecomp
      character(LEN=13) :: filename

      ! MPI & ADIOS variables
      integer :: sizeMB=64
      integer    :: adios_err
      integer(8) :: adios_groupsize, adios_totalsize
      integer(8) :: adios_group, adios_handle, varid
      integer(8) :: offset_x, offset_y, offset_z
      integer :: xdimglob, ydimglob, zdimglob
      integer :: ierr

      ! Init ADIOS
      call ADIOS_Init("adios_BRIO.xml", MPI_COMM_WORLD, ierr)

      ! Metadata definitions
      boxsize   = (/ xdim, ydim, zdim /)
      domdecomp = (/ nx, ny, nz /)

      ! Define offset and global dimensions
      if(inline) then
         offset_x = 0
         offset_y = 0
         offset_z = zdim*myrank
         xdimglob = xdim; ydimglob = ydim; zdimglob = zdim*nx*ny*nz
      else
         offset_x = xdim*xpos
         offset_y = ydim*ypos
         offset_z = zdim*zpos
         xdimglob = xdim*nx; ydimglob = ydim*ny; zdimglob = zdim*nz
      endif

      ! Open ADIOS file & write data
      call ADIOS_Open(adios_handle, "dump", "parallelio.bp", "w" &
           & , MPI_COMM_WORLD, ierr)
      ! Write I/O
#include "gwrite_dump.fh"

      ! Close ADIOS file and interface
      call ADIOS_Close(adios_handle, ierr)
      call ADIOS_Finalize(myrank, ierr)

      return
    end subroutine write_adios_XML
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
         case ('-xml')
            read(arg,*) xml
         case default
            print '("unknown option ",a2," ignored")', opt
         end select
      end do
#else
      namelist /model_params/ xdim, ydim, zdim, nx, ny, nz
      namelist /output_params/ inline
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

      ! ADIOS parameters
      xml = .false.

      return
    end subroutine default_params
    !===========================================================================
    subroutine get_filename(myrank, prefix, filename)
      implicit none

      include "mpif.h"

      integer :: myrank
      character(LEN=*)  :: prefix
      character(LEN=6)  :: nrank
      character(LEN=80) :: filename,filedir,filecmd
      integer :: ierr

      call convtoasc(myrank, nrank)
      filedir = 'sequentialio/'
      filecmd = 'mkdir -p '//trim(filedir)
#if BLUEGENE==1
      if (myrank==0) call mkdir(trim(filedir)//'\0', %val(511))
#else
      if (myrank==0) call system(filecmd)
#endif
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      filename = trim(filedir)//trim(prefix)//'.'//trim(nrank)

      return
    end subroutine get_filename
    !===========================================================================
    subroutine convtoasc(number, nb_string)
      ! To convert an integer smaller than 999999 to a 6 characters string
      implicit none

      integer :: number, istring, num, nums10, i
      character(LEN=6) :: nb_string
      character(LEN=10), parameter :: nstring="0123456789"
     
      num=1000000
      nums10=num/10
      do i = 1, 6
         istring        = 1 + mod(number,num)/nums10
         nb_string(i:i) = nstring(istring:istring)
         num            = num/10
         nums10         = nums10/10
      enddo
  
      return
    end subroutine convtoasc
    !===========================================================================
    !===========================================================================
    !===========================================================================
end program BRIO
