module readIO
  !=============================================================================
  !======================== Benchmark for paRallel I/O =========================
  !=============================================================================
  ! Author: Marc B.R. Joos
  !
  ! Created/last modified: sep 10, 2013/sep 10, 2013
  !
  ! This file is distributed under GNU/GPL licence, 
  ! see <http://www.gnu.org/licenses/>.
  !=============================================================================

contains
!===============================================================================
#if POSIX == 1 || POTOK == 1
  subroutine read_posix(data, xpos, ypos, zpos, myrank)
    use params
    use writeIO
    implicit none

    integer :: xpos, ypos, zpos, myrank
    real(8), dimension(xdim,ydim,zdim,nvar) :: data
    integer, dimension(3) :: boxsize, domdecomp
    character(LEN=80) :: filename
    
    call get_filename(myrank, 'posix', filename)

    open(unit=10, file=filename, status='unknown', form='unformatted')
    read(10) boxsize
    read(10) domdecomp
    read(10) data
    close(10)

    return
  end subroutine read_posix
#endif
!===============================================================================
#if PNCDF == 1
  subroutine read_pncdf(data, xpos, ypos, zpos, myrank)
    use pnetcdf
    use params
    implicit none

    include "mpif.h"

    integer :: xpos, ypos, zpos, myrank
    real(8), dimension(xdim,ydim,zdim,nvar) :: data
    integer, dimension(3) :: boxsize, domdecomp
    character(LEN=13) :: filename
  
    ! PnetCDF variables
    integer(kind=MPI_OFFSET_KIND), dimension(3) :: dims, start, count
    integer :: nout, ncid, ndims, nvars, natts, nulms
    character(Len=20) :: name
    integer :: i, id, diminline, type, ndimvar, dimid(3)
    integer :: ierr

    filename = 'parallelio.nc'

    ! Open netCDF file
    nout = nfmpi_open(MPI_COMM_WORLD, filename &
         , ior(NF_NOWRITE, NF_64BIT_OFFSET) &
         , MPI_INFO_NULL, ncid)

    ! Get basic informations (number of dimensions, variables...)
    nout = nfmpi_inq(ncid, ndims, nvars, natts, nulms)

    ! Loop over variables to retrieve metadata
    do i = 1, nvars
       nout = nfmpi_inq_var(ncid, i, name, type, ndimvar, dimid, natts)

       if(ndimvar.eq.1) then
          if(trim(name).eq."boxsize") then
             nout =  nfmpi_get_vara_int_all(ncid, i, (/1_8/), (/3_8/) &
                  , boxsize)
          endif
          if(trim(name).eq."domdecomp") then
             nout =  nfmpi_get_vara_int_all(ncid, i, (/1_8/), (/3_8/) &
                  , domdecomp)
          endif
       endif
    enddo

    ! Synchronize the threads
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! Loop over variables to retrieve data
    do i = 1, nvars
       nout = nfmpi_inq_var(ncid, i, name, type, ndimvar, dimid, natts)
       
       if(trim(name).eq."var1") id = 1
       if(trim(name).eq."var2") id = 2
       if(trim(name).eq."var3") id = 3
       if(trim(name).eq."var4") id = 4
       if(trim(name).eq."var5") id = 5
       if(trim(name).eq."var6") id = 6
       if(trim(name).eq."var7") id = 7
       if(trim(name).eq."var8") id = 8

       if(inline) then
          diminline = myrank*dims(3) + 1
          start = (/ 1, 1, diminline /)
       else
          start = (/ xpos, ypos, zpos /)*dims + 1
       endif
       count = dims
       
       if(ndimvar.eq.3) then
          nout = nfmpi_get_vara_double_all(ncid, i, start, count &
               , data(:,:,:,id))
       endif
    enddo

    ! Close the file
    nout = nfmpi_close(ncid)

    return
  end subroutine read_pncdf
#endif
!===============================================================================
#if PHDF5 == 1
  subroutine read_phdf5(data, xpos, ypos, zpos, myrank)
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
    integer(HID_T) :: file_id, fapl_id
    integer(HID_T) :: h5_dspace, h5_dset
    integer(HSIZE_T), dimension(3) :: dim_meta
  
    ! Initialize HDF5 interface
    call H5open_f(ierr)
  
    ! Create HDF5 property IDs for parallel file access
    filename = 'parallelio.h5'
    call H5Pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    call H5Pset_fapl_mpio_f(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)
    call H5Fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, ierr &
         , access_prp=fapl_id)

    ! Read metadata
    dim_meta = 3
    call H5Dopen_f(file_id, "boxsize", h5_dset, ierr)
    call H5Dread_f(h5_dset, H5T_NATIVE_INTEGER, boxsize, dim_meta, ierr)
    call H5Dclose_f(h5_dset, ierr)
    call H5Dopen_f(file_id, "domdecomp", h5_dset, ierr)
    call H5Dread_f(h5_dset, H5T_NATIVE_INTEGER, domdecomp, dim_meta, ierr)
    call H5Dclose_f(h5_dset, ierr)

    ! Read data
    call read_dataset_h5(file_id, "var1", data(:,:,:,1), xpos, ypos, zpos &
         , myrank)
    call read_dataset_h5(file_id, "var2", data(:,:,:,2), xpos, ypos, zpos &
         , myrank)
    call read_dataset_h5(file_id, "var3", data(:,:,:,3), xpos, ypos, zpos &
         , myrank)
    call read_dataset_h5(file_id, "var4", data(:,:,:,4), xpos, ypos, zpos &
         , myrank)
    call read_dataset_h5(file_id, "var5", data(:,:,:,5), xpos, ypos, zpos &
         , myrank)
    call read_dataset_h5(file_id, "var6", data(:,:,:,6), xpos, ypos, zpos &
         , myrank)
    call read_dataset_h5(file_id, "var7", data(:,:,:,7), xpos, ypos, zpos &
         , myrank)
    call read_dataset_h5(file_id, "var8", data(:,:,:,8), xpos, ypos, zpos &
         , myrank)

    call H5Fclose_f(file_id, ierr)
    call H5Pclose_f(fapl_id, ierr)
    call H5close_f(ierr)

    return
  end subroutine read_phdf5
!===============================================================================
  subroutine read_dataset_h5(loc_id, dsetname, data, xpos, ypos, zpos, myrank)
    use hdf5
    use params
    implicit none

    character(LEN=*) :: dsetname
    real(8), dimension(xdim,ydim,zdim) :: data
    integer :: xpos, ypos, zpos, myrank
  
    ! HDF5 variables
    integer(HID_T) :: loc_id, dxpl_id, h5_dspace, h5_dset, h5_dspace_file
    integer(HSIZE_T), dimension(3) :: start, count, stride, blockSize
    integer(HSIZE_T), dimension(3) :: dims, dims_file
    integer :: ierr
    
    dims      = (/ xdim, ydim, zdim /)
    if(inline) then
       dims_file = (/ xdim, ydim, zdim*nx*ny*nz /)
    else
       dims_file = (/ xdim*nx, ydim*ny, zdim*nz /)
    endif

    call H5Screate_simple_f(3, dims, h5_dspace, ierr)
    call H5Screate_simple_f(3, dims_file, h5_dspace_file, ierr)

    ! create data set
    call H5Dopen_f(loc_id, trim(dsetname), h5_dset, ierr)
     
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
     
    ! Enable parallel collective IO
    call H5Pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
    call H5Pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)
     
    ! Finally write data to file
    call H5Dread_f(h5_dset, H5T_NATIVE_DOUBLE, data, dims, ierr &
         , mem_space_id=h5_dspace, file_space_id=h5_dspace_file &
         , xfer_prp=dxpl_id)
     
    ! Clean HDF5 IDs
    call H5Pclose_f(dxpl_id, ierr)
    call H5Dclose_f(h5_dset, ierr)
    call H5Sclose_f(h5_dspace, ierr)
    call H5Sclose_f(h5_dspace_file, ierr)

    return
  end subroutine read_dataset_h5
#endif
!===============================================================================
#if ADIOS == 1
  subroutine read_adios(data, xpos, ypos, zpos, myrank)
    use adios_read_mod
    use params
    implicit none
    
    include "mpif.h"
    
    integer :: xpos, ypos, zpos, myrank
    real(8), dimension(xdim,ydim,zdim,nvar) :: data
    integer, dimension(3) :: boxsize, domdecomp
    character(LEN=17) :: filename
    integer :: i

    ! MPI & ADIOS variables
    integer(8) :: adios_handle, sel
    integer :: var_count, att_count, gr_count, tfirst, tlast
    integer :: vtype, vstep, vrank
    integer(8) :: offset_x, offset_y, offset_z
    integer(8), dimension(3) :: dims, start, count
    integer :: totalsize
    character(LEN=64), allocatable, dimension(:) :: vnames, anames, gnames
    integer :: ierr

    if(xml) then
       filename = "parallelio_XML.bp"
    else
       filename = "parallelio_noXML.bp"
    endif

    call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD &
         , "verbose=1", ierr)
    call adios_read_open_file(adios_handle, filename, 0, MPI_COMM_WORLD, ierr)
    call adios_inq_file(adios_handle, var_count, att_count, tfirst, tlast, ierr)

    allocate(vnames(var_count), anames(att_count))
    call adios_inq_varnames(adios_handle, vnames, ierr)
    call adios_inq_attrnames(adios_handle, anames, ierr)

    call adios_inq_ngroups(adios_handle, gr_count, ierr)
    allocate(gnames(gr_count))
    call adios_inq_groupnames(adios_handle, gnames, ierr)

    do i = 1, var_count
       call adios_inq_var(adios_handle, vnames(i), vtype, vstep, vrank, dims &
            , ierr)
       if(vrank .eq. 0) then
          if(vnames(i) .eq. "xdim") call adios_get_scalar(adios_handle &
               , vnames(i), xdim, ierr)
          if(vnames(i) .eq. "ydim") call adios_get_scalar(adios_handle &
               , vnames(i), ydim, ierr)
          if(vnames(i) .eq. "zdim") call adios_get_scalar(adios_handle &
               , vnames(i), zdim, ierr)
       endif
    enddo

    if(inline) then
       offset_x = 0
       offset_y = 0
       offset_z = zdim*myrank
    else
       offset_x = xdim*xpos
       offset_y = ydim*ypos
       offset_z = zdim*zpos
    endif
    start = (/ offset_x, offset_y, offset_z /)
    count = (/ xdim, ydim, zdim /)

    do i = 1, var_count
       call adios_inq_var(adios_handle, vnames(i), vtype, vstep, vrank, dims &
            , ierr)
       if(vrank .eq. 3) then
          call adios_selection_boundingbox(sel, vrank, start, count)
          if(vnames(i) .eq. "var1") call adios_schedule_read(adios_handle, sel &
               , vnames(i), 0, 1, data(:,:,:,1), ierr)
          if(vnames(i) .eq. "var2") call adios_schedule_read(adios_handle, sel &
               , vnames(i), 0, 1, data(:,:,:,2), ierr)
          if(vnames(i) .eq. "var3") call adios_schedule_read(adios_handle, sel &
               , vnames(i), 0, 1, data(:,:,:,3), ierr)
          if(vnames(i) .eq. "var4") call adios_schedule_read(adios_handle, sel &
               , vnames(i), 0, 1, data(:,:,:,4), ierr)
          if(vnames(i) .eq. "var5") call adios_schedule_read(adios_handle, sel &
               , vnames(i), 0, 1, data(:,:,:,5), ierr)
          if(vnames(i) .eq. "var6") call adios_schedule_read(adios_handle, sel &
               , vnames(i), 0, 1, data(:,:,:,6), ierr)
          if(vnames(i) .eq. "var7") call adios_schedule_read(adios_handle, sel &
               , vnames(i), 0, 1, data(:,:,:,7), ierr)
          if(vnames(i) .eq. "var8") call adios_schedule_read(adios_handle, sel &
               , vnames(i), 0, 1, data(:,:,:,8), ierr)
          call adios_perform_reads(adios_handle, ierr)
       endif
    enddo

    call adios_read_close(adios_handle, ierr)

    return
  end subroutine read_adios
#endif
!===============================================================================
end module readIO
