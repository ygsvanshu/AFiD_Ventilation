!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolvePressureCorrection.F90                    !
!    CONTAINS: subroutine SolvePressureCorrection,        !
!     CreateFFTTmpArrays, DestroyFFTTmpArrays             ! 
!                                                         ! 
!    PURPOSE: Compute the pressure correction by solving  !
!     a Poisson equation                                  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SolvePressureCorrection

    use, intrinsic :: iso_c_binding
    use param
    use fftw_params
    use local_arrays, only: dph
    use decomp_2d
    use decomp_2d_fft
    use mpih

    implicit none

    integer :: i,j,k,info
    real :: acphT_b
    real :: appph(nxm-2)
    real, dimension(nxm) :: acphT,apph,amph
    integer :: phpiv(nxm)

    integer(C_FFTW_R2R_KIND), dimension(1):: kind_forw, kind_back
    type(fftw_iodim),dimension(1) :: iodim
    type(fftw_iodim),dimension(2) :: iodim_howmany

    integer mpisize,mpirank,n
    logical fexist
    integer, parameter:: funit = 71

    !RO   Allocate variables for FFT transform

    call CreateFFTTmpArrays

    kind_forw(1) = FFTW_REDFT10
    kind_back(1) = FFTW_REDFT01

    ! warning "initialisation for matrix transpose test"
    ! do i=xstart(3),xend(3)
    !     do j=xstart(2),xend(2)
    !         do k=1,nxm
    !             dph(k,j,i) = j
    !         enddo
    !     enddo
    ! enddo

    call transpose_x_to_y(dph,ry1,ph)

    !RO   Plan FFT transforms if not planned previously

    if (.not.planned) then
        iodim(1)%n=nzm
        iodim(1)%is=(ph%zen(1)-ph%zst(1)+1)*(ph%zen(2)-ph%zst(2)+1)
        iodim(1)%os=(ph%zen(1)-ph%zst(1)+1)*(ph%zen(2)-ph%zst(2)+1)
        iodim_howmany(1)%n=(ph%zen(1)-ph%zst(1)+1)
        iodim_howmany(1)%is=1
        iodim_howmany(1)%os=1
        iodim_howmany(2)%n=(ph%zen(2)-ph%zst(2)+1)
        iodim_howmany(2)%is=(ph%zen(1)-ph%zst(1)+1)
        iodim_howmany(2)%os=(ph%zen(1)-ph%zst(1)+1)
        fwd_guruplan_z=fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,rz2,rz2,kind_forw,FFTW_ESTIMATE)
        iodim(1)%n=nzm
        bwd_guruplan_z=fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,rz2,rz2,kind_back,FFTW_ESTIMATE)

        if (.not.c_associated(bwd_guruplan_z)) then
            if (ismaster) print*,'Failed to create guru plan. You should'
            if (ismaster) print*,'link with FFTW3 before MKL'
            if (ismaster) print*,'Please check linking order.'
            call MPI_Abort(MPI_COMM_WORLD,1,info)
        endif

        iodim(1)%n=nym
        iodim(1)%is=ph%yen(1)-ph%yst(1)+1
        iodim(1)%os=ph%yen(1)-ph%yst(1)+1
        iodim_howmany(1)%n=(ph%yen(1)-ph%yst(1)+1)
        iodim_howmany(1)%is=1
        iodim_howmany(1)%os=1
        iodim_howmany(2)%n=(ph%yen(3)-ph%yst(3)+1)
        iodim_howmany(2)%is=(ph%yen(1)-ph%yst(1)+1)*(ph%yen(2)-ph%yst(2)+1)
        iodim_howmany(2)%os=(ph%yen(1)-ph%yst(1)+1)*(ph%yen(2)-ph%yst(2)+1)
        fwd_guruplan_y=fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,ry1,ry2,kind_forw,FFTW_ESTIMATE)

        iodim(1)%n=nym
        iodim(1)%is=ph%yen(1)-ph%yst(1)+1
        iodim(1)%os=ph%yen(1)-ph%yst(1)+1
        iodim_howmany(1)%n=(ph%yen(1)-ph%yst(1)+1)
        iodim_howmany(1)%is=1
        iodim_howmany(1)%os=1
        iodim_howmany(2)%n=(ph%yen(3)-ph%yst(3)+1)
        iodim_howmany(2)%is=(ph%yen(1)-ph%yst(1)+1)*(ph%yen(2)-ph%yst(2)+1)
        iodim_howmany(2)%os=(ph%yen(1)-ph%yst(1)+1)*(ph%yen(2)-ph%yst(2)+1)
        bwd_guruplan_y=fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,ry2,ry1,kind_back,FFTW_ESTIMATE)
        planned=.true.
    endif

    call dfftw_execute_r2r(fwd_guruplan_y,ry1,ry2)
    call transpose_y_to_z(ry2,rz2,ph)
    call dfftw_execute_r2r(fwd_guruplan_z,rz2,rz2)

    !EP   Normalize. FFT does not do this
    rz2 = 0.25d0*rz2 / (nzm*nym)

    call transpose_z_to_x(rz2,dphr,ph)

    !RO   Solve the tridiagonal matrix with complex coefficients

    !$OMP   PARALLEL DO COLLAPSE(2)                                         &
    !$OMP   DEFAULT(none)                                                   &
    !$OMP   SHARED(sp,nxm)                                                  &
    !$OMP   SHARED(acphk,ak2,ak1,dphc,apphk,amphk)                          &
    !$OMP   PRIVATE(apph,amph,acphT,acphT_b)                                &
    !$OMP   PRIVATE(phpiv,info,appph)

    do i=ph%xst(3),ph%xen(3)
        do j=ph%xst(2),ph%xen(2)
            do k = 1,nxm
                acphT_b=1.0/(acphk(k)-ak2(j)-ak1(i))
                dphr(k,j,i)=dphr(k,j,i)*acphT_b
                apph(k)=apphk(k)*acphT_b
                amph(k)=amphk(k)*acphT_b
                acphT(k)=1.0d0!+ 1.0d-15  ! Small perturbation needed to prevent singular matrix
            enddo                         ! when using uniform grid
        call dgttrf(nxm,amph(2),acphT,apph(1),appph,phpiv,info)
        call dgttrs('N',nxm,1,amph(2),acphT,apph(1),appph,phpiv,dphr(1,j,i),nxm,info)
        enddo
    enddo
    
    !$OMP END PARALLEL DO

    call transpose_x_to_z(dphr,rz2,ph)
    call dfftw_execute_r2r(bwd_guruplan_z,rz2,rz2)
    call transpose_z_to_y(rz2,ry2,ph)
    call dfftw_execute_r2r(bwd_guruplan_y,ry2,ry1)
    call transpose_y_to_x(ry1,dph,ph)
    
    
    ! call mpi_comm_size(mpi_comm_world, mpisize, ierr)
    ! call mpi_comm_rank(mpi_comm_world, mpirank, ierr)
    ! do n = 0, mpisize
    !     if(mpirank .eq. n) then
    !         inquire(file="check.dat", exist=fexist)
    !         if(fexist) then
    !             open(funit, file="check.dat", status="old", position="append", action="write")
    !         else
    !             open(funit, file="check.dat", status="new", action="write")
    !         endif
    !         do i = xstart(3), xend(3)
    !             do j = xstart(2), xend(2)
    !                 do k = 1, nxm
    !                     write(funit,*) k, j, i, dph(k, j, i)
    !                     call flush(6)
    !                 enddo
    !             enddo
    !         enddo
    !         close(funit)
    !     endif
    !     call mpi_barrier(MPI_COMM_WORLD, ierr)
    ! enddo
    ! write(*,*),dph(16,16,:)

    call DestroyFFTTmpArrays

    return

end subroutine SolvePressureCorrection

!=============================================================================
subroutine CreateFFTTmpArrays

    use fftw_params
    use decomp_2d
    use decomp_2d_fft

    implicit none
  
    allocate(ry1(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2),ph%yst(3):ph%yen(3)))
    allocate(rz1(ph%zst(1):ph%zen(1),ph%zst(2):ph%zen(2),ph%zst(3):ph%zen(3)))
    allocate(cy1(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2),sp%yst(3):sp%yen(3)))
    allocate(cz1(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)))
    allocate(dphc(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2),sp%xst(3):sp%xen(3)))

    allocate(ry2(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2),ph%yst(3):ph%yen(3)))
    allocate(rz2(ph%zst(1):ph%zen(1),ph%zst(2):ph%zen(2),ph%zst(3):ph%zen(3)))
    allocate(dphr(ph%xst(1):ph%xen(1),ph%xst(2):ph%xen(2),ph%xst(3):ph%xen(3)))

    return
  
end subroutine CreateFFTTmpArrays
  
!=============================================================================
subroutine DestroyFFTTmpArrays

    use fftw_params

    implicit none

    if(allocated(dphc)) deallocate(dphc)
    if(allocated(rz1)) deallocate(rz1)
    if(allocated(cz1)) deallocate(cz1)
    if(allocated(ry1)) deallocate(ry1)
    if(allocated(cy1)) deallocate(cy1)

    if(allocated(rz2)) deallocate(rz2)
    if(allocated(ry2)) deallocate(ry2)
    if(allocated(dphr)) deallocate(dphr)

    return

end subroutine DestroyFFTTmpArrays