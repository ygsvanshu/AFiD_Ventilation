!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: InitPressureSolver.F90                         !
!    CONTAINS: subroutine InitPressureSolver              !
!                                                         ! 
!    PURPOSE: Initialization routines. Compute the metric !
!     terms and modified wavenumbers for the pressure     !
!     correction                                          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InitPressureSolver

    use param
    use decomp_2d_fft

    implicit none

    integer :: nymh,nymp,j,i,nzmh,nzmp
    integer  :: kc,km,kp
    real :: ugmmm,a33icc,a33icp
    
    ! Initialize wave number definitions

    nymh=nym/2+1
    nymp=nymh+1

    nzmh=nzm/2+1
    nzmp=nzmh+1

    !****** This section was for DFT ******!
    ! do i=1,nzmh
    !     ao(i)=(i-1)*2.d0*pi
    ! enddo
    ! do i=nzmp,nzm
    !     ao(i)=-(nzm-i+1)*2.d0*pi
    ! enddo
    ! do i=1,nzm
    !     ak1(i)=2.d0*(1.d0-dcos(ao(i)/nzm))*(float(nzm)/zlen)**2
    ! enddo

    ! do j=1,nymh
    !     ap(j)=(j-1)*2.d0*pi
    ! enddo
    ! do j=nymp,nym
    !     ap(j)=-(nym-j+1)*2.d0*pi
    ! enddo
    ! do j=1,nym
    !     ak2(j)=2.d0*(1.d0-dcos(ap(j)/nym))*(float(nym)/ylen)**2
    ! enddo
    !****** This section was for DFT ******!

    do i = 1,nzm
        ao(i) = (i-1)*pi
        ak1(i) = 2.d0*(1.d0-dcos(ao(i)/nzm))*(float(nzm)/zlen)**2
    end do

    do j = 1,nym
        ap(j) = (j-1)*pi
        ak2(j) = 2.d0*(1.d0-dcos(ap(j)/nym))*(float(nym)/ylen)**2
    end do

    !RO! Initialize Tridiagonal matrices for Poisson solver

    do kc=1,nxm
        km=kmv(kc)
        kp=kpv(kc)
        a33icc=kmc(kc)*dxq/g3rc(kc)
        a33icp=kpc(kc)*dxq/g3rc(kp)
        ! =ModR10=Robert=2020-11-02=====================================================
        ugmmm=1.0d0/dx3c(kc) ! Nomeclature g3rm->dx3c Stick to old discretization
        ! =End=of=ModR10================================================================
        amphk(kc)=a33icc*ugmmm
        apphk(kc)=a33icp*ugmmm
        acphk(kc)=-(amphk(kc)+apphk(kc))
    enddo

    !RO! Initialize Pencil transposes for pressure solver

    call decomp_2d_fft_init

    return

    end