! Distributed HDF5 hyperslab I/O adapted from QuOp_MPI, native/io/io.f90.
module parallel_io
    use, intrinsic :: iso_c_binding, only: c_int8_t, c_loc
    use, intrinsic :: iso_fortran_env, only: int64
    use hdf5
    use mpi
    implicit none
    private
    public :: open_file, write_dataset, close_file

contains

    logical function failed(error, comm)
        integer, intent(in) :: error, comm
        integer :: ierr
        call MPI_Allreduce(error < 0, failed, 1, MPI_LOGICAL, MPI_LOR, comm, ierr)
    end function

    subroutine open_file(name, comm, file)
        character(len=*), intent(in) :: name
        integer, intent(in) :: comm
        integer(hid_t), intent(out) :: file
        integer(hid_t) :: access
        integer :: error, ierr

        access = -1
        file = -1
        call h5open_f(error)
        if (error == 0) call h5pcreate_f(H5P_FILE_ACCESS_F, access, error)
        if (.not. failed(error, comm)) then
            call h5pset_fapl_mpio_f(access, comm, MPI_INFO_NULL, error)
            if (.not. failed(error, comm)) call h5fopen_f(name, H5F_ACC_RDWR_F, file, error, access)
        end if
        if (access >= 0) call h5pclose_f(access, ierr)
        if (failed(error, comm)) file = -1
    end subroutine

    subroutine write_dataset(file, name, comm, ndim, nblocks, offsets, counts, nbytes, buffer, status)
        integer(hid_t), intent(in) :: file
        character(len=*), intent(in) :: name
        integer, intent(in) :: comm, ndim, nblocks
        integer(int64), intent(in) :: offsets(ndim,nblocks), counts(ndim,nblocks), nbytes
        integer(c_int8_t), intent(in), target :: buffer(nbytes)
        integer, intent(out) :: status
        integer(hid_t) :: dataset, file_space, mem_space, disk_type, mem_type, xfer
        integer(hsize_t) :: count(ndim), offset(ndim), points
        integer(size_t) :: itemsize
        integer(int64) :: first, address
        integer :: i, rank, error, ierr

        dataset = -1
        file_space = -1
        disk_type = -1
        mem_type = -1
        xfer = -1
        first = 1
        call h5dopen_f(file, name, dataset, error)
        if (error == 0) call h5dget_space_f(dataset, file_space, error)
        if (error == 0) then
            call h5sget_simple_extent_ndims_f(file_space, rank, error)
            if (rank /= ndim .or. nbytes < 1) error = -1
        end if
        if (error == 0) call h5dget_type_f(dataset, disk_type, error)
        if (error == 0) call h5tget_native_type_f(disk_type, H5T_DIR_ASCEND_F, mem_type, error)
        if (error == 0) call h5tget_size_f(mem_type, itemsize, error)
        if (error == 0) call h5pcreate_f(H5P_DATASET_XFER_F, xfer, error)
        if (error == 0) call h5pset_dxpl_mpio_f(xfer, H5FD_MPIO_COLLECTIVE_F, error)

        do i = 1, nblocks
            ! HDF5's Fortran interface takes the fastest varying dimension first.
            count = counts(ndim:1:-1,i)
            offset = offsets(ndim:1:-1,i)
            points = product(count)
            mem_space = -1
            address = 1
            if (error == 0) then
                if (points > (nbytes - first + 1) / itemsize) error = -1
            end if
            if (error == 0) call h5screate_simple_f(1, [max(points, 1_hsize_t)], mem_space, error)
            if (points > 0) then
                address = first
                if (error == 0) call h5sselect_hyperslab_f(file_space, H5S_SELECT_SET_F, offset, count, error)
            else
                if (error == 0) call h5sselect_none_f(file_space, error)
                if (error == 0) call h5sselect_none_f(mem_space, error)
            end if
            if (failed(error, comm)) then
                error = -1
            else
                call h5dwrite_f(dataset, mem_type, c_loc(buffer(address)), error, &
                                mem_space_id=mem_space, file_space_id=file_space, xfer_prp=xfer)
            end if
            if (mem_space >= 0) call h5sclose_f(mem_space, ierr)
            if (failed(error, comm)) exit
            first = first + points * itemsize
        end do

        if (xfer >= 0) call h5pclose_f(xfer, ierr)
        if (mem_type >= 0) call h5tclose_f(mem_type, ierr)
        if (disk_type >= 0) call h5tclose_f(disk_type, ierr)
        if (file_space >= 0) call h5sclose_f(file_space, ierr)
        if (dataset >= 0) then
            call h5dclose_f(dataset, ierr)
            error = min(error, ierr)
        end if
        status = 0
        if (failed(error, comm)) status = -1
    end subroutine

    subroutine close_file(file, status)
        integer(hid_t), intent(in) :: file
        integer, intent(out) :: status
        call h5fclose_f(file, status)
    end subroutine
end module
