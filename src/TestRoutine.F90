subroutine TestRoutine

    use mpi
    use param
    use local_arrays, only: rhsx,rhsy,rhsz
    use implicit_decomp

    implicit none

    integer :: k,j,i
    real    :: test_array(dxst(1):dxen(1),dxst(2):dxen(2),dxst(3):dxen(3))
    real    :: count,res_dummy
    character*50 :: filename
    TYPE(DECOMP_INFO) :: testinfo

    call decomp_info_init(nx,nym,nzm,testinfo)

    ! count = 0.0d0
    test_array(:,:,:) = 0.0d0
    do i=dxst(3),dxen(3)
        do j=dxst(2),dxen(2)
            do k=dxst(1),dxen(1)
                test_array(k,j,i) = 1.0d0!count
                rhsx(k,j,i) = 1.0d0!count
                ! count = count + 1.0d0
            end do
        end do
    end do

    filename = trim('Debug/Pre.h5')
    call HdfDebugWriteRealHalo3D(filename,rhsx)
    
    call transpose_x_to_y(rhsx,rhsy,testinfo)
    ! call transpose_y_to_z(rhsy,rhsz,decomp_diff)
    ! call transpose_z_to_y(rhsz,rhsy,decomp_diff)
    call transpose_y_to_x(rhsy,rhsx,testinfo)

    ! call transpose_x_to_y(rhsx,rhsy)
    ! call transpose_y_to_z(rhsy,rhsz)
    ! call transpose_z_to_y(rhsz,rhsy)
    ! call transpose_y_to_x(rhsy,rhsx)

    filename = trim('Debug/Post.h5')
    call HdfDebugWriteRealHalo3D(filename,rhsx)

    count = 0.0d0
    do i=dxst(3),dxen(3)
        do j=dxst(2),dxen(2)
            do k=dxst(1),dxen(1)
                ! k = nx
                ! if ((abs(rhsx(k,j,i)-test_array(k,j,i))).gt.0.1d0) count = count + 1.0d0
                if (abs(rhsx(k,j,i)-1.0d0).gt.1E-6) count = count + 1.0d0
            end do
        end do
    end do

    call MpiAllSumRealScalar(count,res_dummy)
    count = res_dummy

    if (ismaster) then
        if (count.gt.0.0d0) write (*,*) 'TRANSFORMS FAILED FOR ', count, ' ELEMENTS!'
    end if
    
end subroutine TestRoutine