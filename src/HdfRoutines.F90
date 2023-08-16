!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: HdfRoutines.F90                                !
!    CONTAINS: subroutine hdf_read_serial_1d              !
!                                                         ! 
!    PURPOSE: I/O routines.                               !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine HdfStart
    use hdf5
    implicit none
    integer :: hdf_error
    call h5open_f(hdf_error)
end subroutine HdfStart
  
subroutine HdfClose
    use hdf5
    implicit none
    integer :: hdf_error
    call h5close_f(hdf_error)
end subroutine HdfClose

subroutine HdfCreateBlankFile(filename)

    use hdf5

    implicit none

    character*50,intent(in) :: filename
    integer(HID_T) :: file_id
    integer :: hdf_error

    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf_error)
    call h5fclose_f(file_id,hdf_error)

end subroutine HdfCreateBlankFile

subroutine HdfParallelCreateBlankFile(filename,comm)

    use mpih
    use hdf5

    implicit none

    character*50,intent(in) :: filename
    integer,intent(in)      :: comm
    integer(HID_T)          :: file_id,plist_id
    integer                 :: info,hdf_error

    info = MPI_INFO_NULL
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
    call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf_error,access_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    call h5fclose_f(file_id,hdf_error)

end subroutine HdfParallelCreateBlankFile

subroutine HdfCreateGroup(linkname,filename)

    use hdf5

    implicit none

    character*50,intent(in) :: filename,linkname
    integer(HID_T)          :: file_id,group_id
    integer                 :: hdf_error
    integer                 :: info
    logical                 :: lexist

    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf_error)
    call h5lexists_f(file_id,linkname,lexist,hdf_error)
    if (.not.lexist) then
        call h5gcreate_f(file_id,linkname,group_id,hdf_error)
        call h5gclose_f(group_id,hdf_error)
    end if
    call h5fclose_f(file_id,hdf_error)

end subroutine HdfCreateGroup

subroutine HdfParallelCreateGroup(linkname,filename,comm)

    use mpih
    use hdf5

    implicit none

    character*50, intent(in)        :: filename,linkname
    integer, intent(in), optional   :: comm
    integer(HID_T)                  :: file_id,group_id,plist_id
    integer                         :: hdf_error
    integer                         :: info
    logical                         :: lexist

    info = MPI_INFO_NULL
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
    call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf_error,access_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    call h5lexists_f(file_id,linkname,lexist,hdf_error)
    if (.not.lexist) then
        call h5gcreate_f(file_id,linkname,group_id,hdf_error)
        call h5gclose_f(group_id,hdf_error)
    end if
    call h5fclose_f(file_id,hdf_error)

end subroutine HdfParallelCreateGroup

! If MPI communicator comm is given, it is assumed that file is opened for the appropriate communicator
! If MPI communicator comm is not given, it is assumed that the subroutine is called only on root node 
subroutine HdfCreatePath(dsetname,filename,comm)

    use hdf5
    use mpih
    use decomp_2d, only:nrank

    implicit none

    character*50, intent(in) :: filename,dsetname
    integer, intent(in), optional :: comm
    
    integer(HID_T) :: file_id,group_id
    integer(HID_T) :: plist_id
    integer :: hdf_error,ierror,ndims
    integer :: info
    integer :: slashpos
    logical :: fexist,lexist
    character*50 :: linkname,grupname,chckname

    if (present(comm)) then
        ! write (6,*) 'Communicator is present' 
        info = MPI_INFO_NULL
        call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
        call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
        inquire(file=filename,exist=fexist)
        if (fexist) then
            ! write (6,*) 'File ',trim(filename),' exists' 
            call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf_error,access_prp=plist_id)
            ! write (6,*) 'Opened file ',trim(filename) 
        else
            ! write (6,*) 'File ',trim(filename),' does not exist' 
            call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf_error,access_prp=plist_id)
            ! write (6,*) 'Created file ',trim(filename) 
        end if
        call h5pclose_f(plist_id,hdf_error)
    else
        ! write (6,*) 'Communicator is absent' 
        inquire(file=filename,exist=fexist)
        if (fexist) then
            ! write (6,*) 'File ',trim(filename),' exists' 
            call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf_error)
            ! write (6,*) 'Opened file ',filename 
        else
            ! write (6,*) 'File ',trim(filename),' does not exist' 
            call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf_error)
            ! write (6,*) 'Created file ',trim(filename) 
        end if
    end if

    slashpos = 1
    linkname = dsetname
    chckname = trim('')
    do while (slashpos.gt.0)
        slashpos = index(linkname,'/')
        if (slashpos.ne.0) then
            grupname = trim(linkname(:slashpos-1))
            linkname = trim(linkname(slashpos+1:))
            chckname = trim(chckname)//'/'//trim(grupname)
            ! write (6,*) 'Checking if group ',trim(chckname),' is present' 
            call h5lexists_f(file_id,trim(chckname),lexist,hdf_error)
            if (.not.lexist) then
                ! write (6,*) 'Group ',trim(chckname),' does not exist' 
                call h5gcreate_f(file_id,trim(chckname),group_id,hdf_error)
                call h5gclose_f(group_id,hdf_error)
                ! write (6,*) 'Created group ',trim(chckname)
            ! else
            !     call h5gopen_f(file_id,trim(chckname),group_id,hdf_error)
            !     call h5gclose_f(group_id,hdf_error)
            end if
        end if
    end do

    call h5fclose_f(file_id,hdf_error)

    ! call MPI_BARRIER(comm,ierror)

    ! write (6,*) 'Created path by ',nrank

end subroutine HdfCreatePath

!==================================================================================================

subroutine HdfSerialWriteRealScalar(dsetname,filename,n)

    use hdf5

    implicit none

    character*50,intent(in) :: dsetname,filename
    real,intent(in) :: n
    integer(HID_T) :: file_id
    integer(HID_T) :: dset,filespace
    integer :: hdf_error
    integer(HSIZE_T) :: dims(1)
    logical :: fileexists

    dims(1)=1

    inquire(file=filename,exist=fileexists)
    if (fileexists) then
        call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf_error)
    else
        call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf_error)
    end if
    call h5lexists_f(file_id,dsetname,fileexists,hdf_error)
    if (fileexists) call h5ldelete_f(file_id,dsetname,hdf_error)
    call h5screate_simple_f(1,dims,filespace,hdf_error)
    call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset,hdf_error)
    call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,n,dims,hdf_error)
    call h5dclose_f(dset,hdf_error)
    call h5sclose_f(filespace,hdf_error)
    call h5fclose_f(file_id,hdf_error)
      
end subroutine HdfSerialWriteRealScalar

subroutine HdfSerialWriteIntScalar(dsetname,filename,n)

    use hdf5

    implicit none
    
    character*50,intent(in) :: dsetname,filename
    integer,intent(in) :: n
    integer(HID_T) :: file_id
    integer(HID_T) :: dset,filespace
    integer :: hdf_error
    integer(HSIZE_T) :: dims(1)
    logical :: fileexists

    dims(1)=1

    inquire(file=filename,exist=fileexists)
    if (fileexists) then
        call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf_error)
    else
        call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf_error)
    end if
    call h5lexists_f(file_id,dsetname,fileexists,hdf_error)
    if (fileexists) call h5ldelete_f(file_id,dsetname,hdf_error)
    call h5screate_simple_f(1,dims,filespace,hdf_error)
    call h5dcreate_f(file_id,dsetname,H5T_NATIVE_INTEGER,filespace,dset,hdf_error)
    call h5dwrite_f(dset,H5T_NATIVE_INTEGER,n,dims,hdf_error)
    call h5dclose_f(dset,hdf_error)
    call h5sclose_f(filespace,hdf_error)
    call h5fclose_f(file_id,hdf_error)
      
end subroutine HdfSerialWriteIntScalar

subroutine HdfSerialWriteReal1D(dsetname,filename,var,sz)

    use hdf5

    implicit none

    character*50,intent(in) :: dsetname,filename
    integer,intent(in) :: sz
    real,dimension(sz),intent(in) :: var
    integer(HID_T) :: file_id
    integer(HID_T) :: dset,filespace
    integer :: hdf_error
    integer(HSIZE_T) :: dims(1)
    logical :: fileexists

    dims(1)=sz

    inquire(file=filename,exist=fileexists)
    if (fileexists) then
        call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf_error)
    else
        call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf_error)
    end if
    call h5screate_simple_f(1,dims,filespace,hdf_error)
    call h5lexists_f(file_id,dsetname,fileexists,hdf_error)
    if (fileexists) call h5ldelete_f(file_id,dsetname,hdf_error)
    call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset,hdf_error)
    call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,var(1:sz),dims,hdf_error)
    call h5dclose_f(dset,hdf_error)
    call h5sclose_f(filespace,hdf_error)
    call h5fclose_f(file_id,hdf_error)
    
end subroutine HdfSerialWriteReal1D

subroutine HdfWriteReal2D_X(dsetname,filename,var)

    use mpih
    use param
    use hdf5
    use decomp_2d,only: xstart,xend

    implicit none

    character*50,intent(in)             :: dsetname,filename
    real,intent(in)                     :: var(xstart(2):xend(2),xstart(3):xend(3))
    integer(HID_T)                      :: file_id,group_id,plist_id
    integer(HID_T)                      :: filespace
    integer(HID_T)                      :: memspace
    integer(HID_T)                      :: dset
    integer(HSIZE_T)                    :: dims(2)
    integer(HSIZE_T),dimension(2)       :: data_count  
    integer(HSIZE_T),dimension(2)       :: data_offset 
    integer                             :: hdf_error,ndims
    integer                             :: comm,info
    logical                             :: dexist

    !RO   Sort out MPI definitions

    comm = comm_xcut
    info = MPI_INFO_NULL

    ndims = 2
    dims(1)=nym
    dims(2)=nzm

    data_count(1) = xend(2)-xstart(2)+1
    data_count(2) = xend(3)-xstart(3)+1

    data_offset(1) = xstart(2)-1
    data_offset(2) = xstart(3)-1
    
    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
    call h5pset_fapl_mpio_f(plist_id,comm,info,hdf_error)
    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf_error,access_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    call h5screate_simple_f(ndims,dims,filespace,hdf_error)
    call h5screate_simple_f(ndims,data_count,memspace,hdf_error)
    call h5lexists_f(file_id,dsetname,dexist,hdf_error)
    if (dexist) call h5ldelete_f(file_id,dsetname,hdf_error)
    call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset,hdf_error)
    call h5dget_space_f(dset,filespace,hdf_error)
    call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
    call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,var(xstart(2):xend(2),xstart(3):xend(3)),data_count,hdf_error,file_space_id=filespace,mem_space_id=memspace,xfer_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    call h5dclose_f(dset,hdf_error)
    call h5sclose_f(memspace,hdf_error)
    call h5sclose_f(filespace,hdf_error)
    call h5fclose_f(file_id,hdf_error)
      
end subroutine HdfWriteReal2D_X

! Ensure that only certain pencils (MPI ranks)
! containing the required plane use this subroutine
subroutine HdfWriteReal2D_Y(dsetname,filename,var)

    use mpih
    use param
    use hdf5
    use decomp_2d,only: xstart,xend

    implicit none

    character*50,intent(in)             :: dsetname,filename
    real,intent(in)                     :: var(1:nx,xstart(3):xend(3))
    integer(HID_T)                      :: file_id,group_id,plist_id
    integer(HID_T)                      :: filespace
    integer(HID_T)                      :: memspace
    integer(HID_T)                      :: dset
    integer(HSIZE_T)                    :: dims(2)
    integer(HSIZE_T),dimension(2)       :: data_count  
    integer(HSIZE_T),dimension(2)       :: data_offset 
    integer                             :: hdf_error,ndims
    integer                             :: comm,info
    logical                             :: dexist

    !RO   Sort out MPI definitions

    comm = comm_ycut
    info = MPI_INFO_NULL

    ndims = 2
    dims(1)=nx
    dims(2)=nzm

    data_count(1) = nx
    data_count(2) = xend(3)-xstart(3)+1

    data_offset(1) = 0
    data_offset(2) = xstart(3)-1

    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
    call h5pset_fapl_mpio_f(plist_id,comm,info,hdf_error)
    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf_error,access_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    call h5screate_simple_f(ndims,dims,filespace,hdf_error)
    call h5screate_simple_f(ndims,data_count,memspace,hdf_error) 
    call h5lexists_f(file_id,dsetname,dexist,hdf_error)
    if (dexist) call h5ldelete_f(file_id,dsetname,hdf_error)
    call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset,hdf_error)
    call h5dget_space_f(dset,filespace,hdf_error)
    call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
    call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,var(1:nx,xstart(3):xend(3)),data_count,hdf_error,file_space_id=filespace,mem_space_id=memspace,xfer_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    call h5dclose_f(dset,hdf_error)
    call h5sclose_f(memspace,hdf_error)
    call h5sclose_f(filespace,hdf_error)
    call h5fclose_f(file_id,hdf_error)
      
end subroutine HdfWriteReal2D_Y

! Ensure that only certain pencils (MPI ranks)
! containing the required plane use this subroutine
subroutine HdfWriteReal2D_Z(dsetname,filename,var)

    use mpih
    use param
    use hdf5
    use decomp_2d,only: xstart,xend

    implicit none

    character*50,intent(in)             :: dsetname,filename
    real,intent(in)                     :: var(1:nx,xstart(2):xend(2))
    integer(HID_T)                      :: file_id,group_id,plist_id
    integer(HID_T)                      :: filespace
    integer(HID_T)                      :: memspace
    integer(HID_T)                      :: dset
    integer(HSIZE_T)                    :: dims(2)
    integer(HSIZE_T),dimension(2)       :: data_count  
    integer(HSIZE_T),dimension(2)       :: data_offset 
    integer                             :: hdf_error,ndims
    integer                             :: comm,info
    logical                             :: dexist

    !RO   Sort out MPI definitions

    comm = comm_zcut
    info = MPI_INFO_NULL

    ndims = 2
    dims(1)=nx
    dims(2)=nym

    data_count(1) = nx
    data_count(2) = xend(2)-xstart(2)+1

    data_offset(1) = 0
    data_offset(2) = xstart(2)-1

    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
    call h5pset_fapl_mpio_f(plist_id,comm,info,hdf_error)
    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf_error,access_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    call h5screate_simple_f(ndims,dims,filespace,hdf_error)
    call h5screate_simple_f(ndims,data_count,memspace,hdf_error)
    call h5lexists_f(file_id,dsetname,dexist,hdf_error)
    if (dexist) call h5ldelete_f(file_id,dsetname,hdf_error)
    call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset,hdf_error)
    call h5dget_space_f(dset,filespace,hdf_error)
    call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
    call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,var(1:nx,xstart(2):xend(2)),data_count,hdf_error,file_space_id=filespace,mem_space_id=memspace,xfer_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    call h5dclose_f(dset,hdf_error)
    call h5sclose_f(memspace,hdf_error)
    call h5sclose_f(filespace,hdf_error)
    call h5fclose_f(file_id,hdf_error)
      
end subroutine HdfWriteReal2D_Z

subroutine HdfWriteRealHalo3D(filename,qua)

    use param
    use mpih
    use hdf5
    use decomp_2d,only: xstart,xend
      
    implicit none

    integer :: hdf_error
    integer(HID_T) :: file_id
    integer(HID_T) :: filespace
    integer(HID_T) :: slabspace
    integer(HID_T) :: memspace
    integer(HID_T) :: dset
    integer(HSIZE_T) :: dims(3)
    integer(HID_T) :: plist_id
    integer(HSIZE_T),dimension(3) :: data_count  
    integer(HSSIZE_T),dimension(3) :: data_offset 
    integer :: comm,info
    integer :: ndims
    real,intent(in),dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: qua
    character*50,intent(in) :: filename

    !RO   Sort out MPI definitions

    comm = MPI_COMM_WORLD
    info = MPI_INFO_NULL

    !RO   Set offsets and element counts
   
    ndims = 3

    dims(1)=nx
    dims(2)=nym
    dims(3)=nzm

    data_count(1) = nx
    data_count(2) = xend(2)-xstart(2)+1
    data_count(3) = xend(3)-xstart(3)+1

    data_offset(1) = 0
    data_offset(2) = xstart(2)-1
    data_offset(3) = xstart(3)-1

    call h5screate_simple_f(ndims,dims,filespace,hdf_error)
    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
    call h5pset_fapl_mpio_f(plist_id,comm,info,hdf_error)
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf_error,access_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    call h5dcreate_f(file_id,'var',H5T_NATIVE_DOUBLE,filespace,dset,hdf_error)
    call h5screate_simple_f(ndims,data_count,memspace,hdf_error) 
    call h5dget_space_f(dset,slabspace,hdf_error)
    call h5sselect_hyperslab_f (slabspace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
    call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,qua(1:nx,xstart(2):xend(2),xstart(3):xend(3)),dims,hdf_error,file_space_id=slabspace,mem_space_id=memspace,xfer_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    call h5dclose_f(dset,hdf_error)
    call h5sclose_f(memspace,hdf_error)
    call h5fclose_f(file_id,hdf_error)

end subroutine HdfWriteRealHalo3D

!==================================================================================================

subroutine HdfSerialReadRealScalar(dsetname,filename,n)

    use hdf5

    implicit none

    character*50,intent(in) :: dsetname,filename
    real,intent(out) :: n
    integer(HID_T) :: file_id
    integer(HID_T) :: dset
    integer :: hdf_error
    integer(HSIZE_T) :: dims(1)

    dims(1)=1
    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf_error)
    call h5dopen_f(file_id,dsetname,dset,hdf_error)
    call h5dread_f(dset,H5T_NATIVE_DOUBLE,n,dims,hdf_error)
    call h5dclose_f(dset,hdf_error)
    call h5fclose_f(file_id,hdf_error)
      
end subroutine HdfSerialReadRealScalar

subroutine HdfSerialReadIntScalar(dsetname,filename,n)

    use hdf5

    implicit none
    
    character*50,intent(in) :: dsetname,filename
    integer,intent(out) :: n
    integer(HID_T) :: file_id
    integer(HID_T) :: dset
    integer :: hdf_error
    integer(HSIZE_T) :: dims(1)

    dims(1)=1

    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf_error)
    call h5dopen_f(file_id,dsetname,dset,hdf_error)
    call h5dread_f(dset,H5T_NATIVE_INTEGER,n,dims,hdf_error)
    call h5dclose_f(dset,hdf_error)
    call h5fclose_f(file_id,hdf_error)
    
end subroutine HdfSerialReadIntScalar

subroutine HdfSerialReadReal1D(dsetname,filename,var,sz)

    use hdf5

    implicit none

    character*50,intent(in) :: dsetname,filename
    integer,intent(in) :: sz
    real,dimension(sz),intent(out) :: var
    integer(HID_T) :: file_id
    integer(HID_T) :: dset
    integer :: hdf_error
    integer(HSIZE_T) :: dims(1)

    dims(1)=sz

    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf_error)
    call h5dopen_f(file_id,dsetname,dset,hdf_error)
    call h5dread_f(dset,H5T_NATIVE_DOUBLE,var,dims,hdf_error)
    call h5dclose_f(dset,hdf_error)
    call h5fclose_f(file_id,hdf_error)
      
end subroutine HdfSerialReadReal1D

subroutine HdfReadReal2D_X(filename,dsetname,qua,st2,en2,st3,en3,n2o,n3o)
    
    use param
    use mpih
    use hdf5
    
    implicit none

    integer,intent(in) :: n2o,n3o,st2,en2,st3,en3
    real,intent(out),dimension(st2:en2,st3:en3) :: qua
    integer :: hdf_error
    integer(HID_T) :: file_id
    integer(HID_T) :: slabspace
    integer(HID_T) :: memspace
    integer(HID_T) :: dset_qua
    integer(HSIZE_T) :: dims(2)
    integer(HID_T) :: plist_id
    integer(HSIZE_T),dimension(2) :: data_count  
    integer(HSSIZE_T),dimension(2) :: data_offset 
    integer :: comm,info
    integer :: ndims
    character*50,intent(in) :: filename,dsetname

    !RO   Sort out MPI definitions

    comm = comm_xcut
    info = MPI_INFO_NULL

    !RO   Set offsets and element counts

    ndims = 2

    dims(1)=n2o-1
    dims(2)=n3o-1

    data_count(1) = en2-st2+1
    data_count(2) = en3-st3+1

    data_offset(1) = st2-1
    data_offset(2) = st3-1

    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
    call h5pset_fapl_mpio_f(plist_id,comm,info,hdf_error)
    call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf_error,access_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    call h5dopen_f(file_id,dsetname,dset_qua,hdf_error)
    call h5screate_simple_f(ndims,data_count,memspace,hdf_error) 
    call h5dget_space_f(dset_qua,slabspace,hdf_error)
    call h5sselect_hyperslab_f (slabspace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
    call h5dread_f(dset_qua,H5T_NATIVE_DOUBLE,qua(st2:en2,st3:en3),dims,hdf_error,file_space_id=slabspace,mem_space_id=memspace,xfer_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    call h5dclose_f(dset_qua,hdf_error)
    call h5sclose_f(memspace,hdf_error)
    call h5fclose_f(file_id,hdf_error)

end subroutine HdfReadReal2D_X

! Ensure that only certain pencils (MPI ranks)
! containing the required plane use this subroutine
subroutine HdfReadReal2D_Y(filename,dsetname,qua,st3,en3,n1o,n3o)
    
    use param
    use mpih
    use hdf5
    
    implicit none

    integer,intent(in) :: n1o,n3o,st3,en3
    real,intent(out),dimension(n1o,st3:en3) :: qua
    integer :: hdf_error,ierror
    integer(HID_T) :: file_id
    integer(HID_T) :: slabspace
    integer(HID_T) :: memspace
    integer(HID_T) :: dset_qua
    integer(HSIZE_T) :: dims(2)
    integer(HID_T) :: plist_id
    integer(HSIZE_T),dimension(2) :: data_count  
    integer(HSSIZE_T),dimension(2) :: data_offset 
    integer :: comm,info
    integer :: ndims
    character*50,intent(in) :: filename,dsetname

    !RO   Sort out MPI definitions

    comm = comm_ycut
    info = MPI_INFO_NULL

    !RO   Set offsets and element counts

    ndims = 2

    dims(1)=n1o
    dims(2)=n3o-1

    data_count(1) = n1o
    data_count(2) = en3-st3+1

    data_offset(1) = 0
    data_offset(2) = st3-1

    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
    call h5pset_fapl_mpio_f(plist_id,comm,info,hdf_error)
    call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf_error,access_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    call h5dopen_f(file_id,dsetname,dset_qua,hdf_error)
    call h5screate_simple_f(ndims,data_count,memspace,hdf_error) 
    call h5dget_space_f(dset_qua,slabspace,hdf_error)
    call h5sselect_hyperslab_f (slabspace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
    call h5dread_f(dset_qua,H5T_NATIVE_DOUBLE,qua(1:nx,st3:en3),dims,hdf_error,file_space_id=slabspace,mem_space_id=memspace,xfer_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    call h5dclose_f(dset_qua,hdf_error)
    call h5sclose_f(memspace,hdf_error)
    call h5fclose_f(file_id,hdf_error)

end subroutine HdfReadReal2D_Y

! Ensure that only certain pencils (MPI ranks)
! containing the required plane use this subroutine
subroutine HdfReadReal2D_Z(filename,dsetname,qua,st2,en2,n1o,n2o)
    
    use param
    use mpih
    use hdf5
    
    implicit none

    integer,intent(in) :: n1o,n2o,st2,en2
    real,intent(out),dimension(n1o,st2:en2) :: qua
    integer :: hdf_error,ierror
    integer(HID_T) :: file_id
    integer(HID_T) :: slabspace
    integer(HID_T) :: memspace
    integer(HID_T) :: dset_qua
    integer(HSIZE_T) :: dims(2)
    integer(HID_T) :: plist_id
    integer(HSIZE_T),dimension(2) :: data_count  
    integer(HSSIZE_T),dimension(2) :: data_offset 
    integer :: comm,info
    integer :: ndims
    character*50,intent(in) :: filename,dsetname

    !RO   Sort out MPI definitions

    comm = comm_zcut
    info = MPI_INFO_NULL

    !RO   Set offsets and element counts

    ndims = 2

    dims(1)=n1o
    dims(2)=n2o-1

    data_count(1) = n1o
    data_count(2) = en2-st2+1

    data_offset(1) = 0
    data_offset(2) = st2-1

    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
    call h5pset_fapl_mpio_f(plist_id,comm,info,hdf_error)
    call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf_error,access_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    call h5dopen_f(file_id,dsetname,dset_qua,hdf_error)
    call h5screate_simple_f(ndims,data_count,memspace,hdf_error) 
    call h5dget_space_f(dset_qua,slabspace,hdf_error)
    call h5sselect_hyperslab_f (slabspace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
    call h5dread_f(dset_qua,H5T_NATIVE_DOUBLE,qua(1:nx,st2:en2),dims,hdf_error,file_space_id=slabspace,mem_space_id=memspace,xfer_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    call h5dclose_f(dset_qua,hdf_error)
    call h5sclose_f(memspace,hdf_error)
    call h5fclose_f(file_id,hdf_error)

end subroutine HdfReadReal2D_Z

subroutine HdfReadRealHalo3D(filename,qua,st2,en2,st3,en3,n1o,n2o,n3o)
    
    use param
    use mpih
    use hdf5
    
    implicit none

    integer,intent(in) :: n1o,n2o,n3o,st2,en2,st3,en3
    real,intent(out),dimension(1:n1o,st2-lvlhalo:en2+lvlhalo,st3-lvlhalo:en3+lvlhalo) :: qua
    integer :: hdf_error
    integer(HID_T) :: file_id
    integer(HID_T) :: slabspace
    integer(HID_T) :: memspace
    integer(HID_T) :: dset_qua
    integer(HSIZE_T) :: dims(3)
    integer(HID_T) :: plist_id
    integer(HSIZE_T),dimension(3) :: data_count  
    integer(HSSIZE_T),dimension(3) :: data_offset 
    integer :: comm,info
    integer :: ndims
    character*50,intent(in) :: filename

    !RO   Sort out MPI definitions

    comm = MPI_COMM_WORLD
    info = MPI_INFO_NULL

    !RO   Set offsets and element counts

    ndims = 3

    dims(1)=n1o
    dims(2)=n2o-1
    dims(3)=n3o-1

    data_count(1) = n1o
    data_count(2) = en2-st2+1
    data_count(3) = en3-st3+1

    data_offset(1) = 0
    data_offset(2) = st2-1
    data_offset(3) = st3-1

    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
    call h5pset_fapl_mpio_f(plist_id,comm,info,hdf_error)
    call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf_error,access_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    call h5dopen_f(file_id,'var',dset_qua,hdf_error)
    call h5screate_simple_f(ndims,data_count,memspace,hdf_error) 
    call h5dget_space_f(dset_qua,slabspace,hdf_error)
    call h5sselect_hyperslab_f (slabspace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
    call h5dread_f(dset_qua,H5T_NATIVE_DOUBLE,qua(1:n1o,st2:en2,st3:en3),dims,hdf_error,file_space_id=slabspace,mem_space_id=memspace,xfer_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    call h5dclose_f(dset_qua,hdf_error)
    call h5sclose_f(memspace,hdf_error)
    call h5fclose_f(file_id,hdf_error)

end subroutine HdfReadRealHalo3D
