module writeIO
  !=============================================================================
  !======================== Benchmark for paRallel I/O =========================
  !=============================================================================
  ! Author: Marc B.R. Joos
  !
  ! Created/last modified: sep 10, 2013/oct 24, 2013
  !
  ! This file is distributed under GNU/GPL license, 
  ! see <http://www.gnu.org/licenses/>.
  !=============================================================================

contains
!===============================================================================
#if POSIX == 1 || POTOK == 1
  subroutine write_posix(data, xpos, ypos, zpos, myrank)
    use params
    use filehandle
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
#endif
!===============================================================================
#if PNCDF == 1
  subroutine write_pncdf(data, xpos, ypos, zpos, myrank)
    use pnetcdf
    use params
    use mpi
    implicit none
    
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
#endif
!===============================================================================
#if PHDF5 == 1
  subroutine write_phdf5(data, xpos, ypos, zpos, myrank)
    use hdf5
    use params
    use mpi
    implicit none
    
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
    if(myrank.eq.0) then
       call H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr)
       
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
       
       call H5Fclose_f(file_id, ierr)
    endif
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
  
    call H5Pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    call H5Pset_fapl_mpio_f(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)
    call H5Fopen_f(filename, H5F_ACC_RDWR_F, file_id, ierr &
         , access_prp=fapl_id)
  
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
!===============================================================================
  subroutine write_dataset_h5(loc_id, dsetname, data, xpos, ypos, zpos&
       , myrank)
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
     
    ! Create data set
    call H5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_DOUBLE &
         , h5_dspace_file, h5_dset, ierr, H5P_DEFAULT_F, H5P_DEFAULT_F &
         , H5P_DEFAULT_F)
  
    ! Finally write data to file
    call H5Dwrite_f(h5_dset, H5T_NATIVE_DOUBLE, data, dims, ierr &
         , mem_space_id=h5_dspace, file_space_id=h5_dspace_file &
         , xfer_prp=dxpl_id)
     
    ! Clean HDF5 IDs
    call H5Pclose_f(dxpl_id, ierr)
    call H5Dclose_f(h5_dset, ierr)
    call H5Sclose_f(h5_dspace, ierr)
    call H5Sclose_f(h5_dspace_file, ierr)
  
    return
  end subroutine write_dataset_h5
#endif
!===============================================================================
#if ADIOS == 1
  subroutine write_adios_noXML(data, xpos, ypos, zpos, myrank)
    use adios_write_mod
    use params
    use mpi
    implicit none
    
    integer :: xpos, ypos, zpos, myrank
    real(8), dimension(xdim,ydim,zdim,nvar) :: data
    integer, dimension(3) :: boxsize, domdecomp
    character(LEN=19) :: filename
  
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
    filename = "parallelio_noXML.bp"
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
!===============================================================================
  subroutine write_adios_XML(data, xpos, ypos, zpos, myrank)
    use adios_write_mod
    use params
    use mpi
    implicit none
    
    integer :: xpos, ypos, zpos, myrank
    real(8), dimension(xdim,ydim,zdim,nvar) :: data
    integer, dimension(3) :: boxsize, domdecomp
    character(LEN=17) :: filename
  
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
    call ADIOS_Open(adios_handle, "dump", "parallelio_XML.bp", "w" &
         & , MPI_COMM_WORLD, ierr)
    ! Write I/O
#include "gwrite_dump.fh"
  
    ! Close ADIOS file and interface
    call ADIOS_Close(adios_handle, ierr)
    call ADIOS_Finalize(myrank, ierr)
  
    return
  end subroutine write_adios_XML
#endif
!===============================================================================
#if MPIIO == 1
  subroutine write_mpiio(data, xpos, ypos, zpos, myrank)
    use params
    use mpi
    implicit none
    
    integer :: xpos, ypos, zpos, myrank, i
    real(8), dimension(xdim,ydim,zdim,nvar) :: data
    integer, dimension(3) :: boxsize, domdecomp
    character(LEN=13) :: filename
  
    ! MPI variables
    integer :: fhandle, ierr
    integer :: int_size, double_size
    integer(kind=MPI_OFFSET_KIND) :: buf_size
    integer :: written_arr
    integer, dimension(3) :: wa_size, wa_subsize, wa_start
  
    ! Metadata definitions
    boxsize   = (/ xdim, ydim, zdim /)
    domdecomp = (/ nx, ny, nz /)
  
    ! Create MPI array type
    wa_size    = (/ nx*xdim, ny*ydim, nz*zdim /)
    wa_subsize = (/ xdim, ydim, zdim /)
    wa_start   = (/ xpos, ypos, zpos /)*wa_subsize
    call MPI_Type_Create_Subarray(3, wa_size, wa_subsize, wa_start &
         & , MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, written_arr, ierr)
    call MPI_Type_Commit(written_arr, ierr)
  
    call MPI_Type_Size(MPI_INTEGER, int_size, ierr)
    call MPI_Type_Size(MPI_DOUBLE_PRECISION, double_size, ierr)
       
    filename = 'parallelio.mp'
    ! Write metadata
    if(myrank.eq.0) then
       call MPI_File_Open(MPI_COMM_SELF, trim(filename), &
            & MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, fhandle, ierr)
       buf_size = 0
       call MPI_File_Seek(fhandle, buf_size, MPI_SEEK_SET, ierr)
       call MPI_File_Write(fhandle, boxsize, 3, MPI_INTEGER &
            & , MPI_STATUS_IGNORE, ierr)
       call MPI_File_Write(fhandle, domdecomp, 3, MPI_INTEGER &
            & , MPI_STATUS_IGNORE, ierr)
       call MPI_File_Close(fhandle, ierr)
    endif
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
  
    ! Open file
    call MPI_File_Open(MPI_COMM_WORLD, trim(filename) &
         & , MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fhandle, ierr)
  
    ! Write data
    do i = 1, 8
       buf_size = 6*int_size + xdim*ydim*zdim*double_size*myrank &
            & + (i-1)*xdim*ydim*zdim*nx*ny*nz*double_size
       call MPI_File_Seek(fhandle, buf_size, MPI_SEEK_SET, ierr)
       call MPI_File_Write_All(fhandle, data(:,:,:,i), xdim*ydim*zdim &
            & , MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
    enddo
  
    ! Close file
    call MPI_File_Close(fhandle, ierr)
  
    return
  end subroutine write_mpiio
#endif
!===============================================================================
  subroutine create_dir(myrank)
    use mpi
    implicit none
    
    integer :: myrank
    character(LEN=80) :: filedir, filecmd
    integer :: ierr
    
    filedir = 'sequentialio/'
    filecmd = 'mkdir -p '//trim(filedir)
#if BLUEGENE==1
    if (myrank==0) call mkdir(trim(filedir)//'\0', %val(511))
#else
    if (myrank==0) call system(filecmd)
#endif
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    
    return
  end subroutine create_dir
!===============================================================================
end module writeIO
