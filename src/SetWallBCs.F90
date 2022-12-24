!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SetWallBCs.F90                                 !
!    CONTAINS: subroutines SetWallBCs, SetInletBC         !
!              and SetOutletBC                            ! 
!                                                         !
!    PURPOSE: Sets the wall boundary conditions including ! 
!             the inlet and outlet boundary conditions    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetWallBCs

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays
    use mpih

    implicit none

    vx(1,:,:) = 0.d0
    vx(nx,:,:) = 0.d0

    ! Just plain side walls, no vents
    if (xstart(2).eq.1) then
        vx(:,0,:) = -vx(:,1,:)                  ! Dirchlet condition
        vy(:,1,:) = 0.d0                        ! Dirchlet condition
        vz(:,0,:) = -vz(:,1,:)                  ! Dirchlet condition
        temp(:,0,:) = temp(:,1,:)               ! Adiabatic
        co2(:,0,:) = co2(:,1,:)                 ! Zero flux
        h2o(:,0,:) = h2o(:,1,:)                 ! Zero flux
    end if

    if (xend(2).eq.nym) then
        vx(:,ny,:) = -vx(:,nym,:)               ! Dirchlet condition
        vy(:,ny,:) = 0.d0                       ! Dirchlet condition
        vz(:,ny,:) = -vz(:,nym,:)               ! Dirchlet condition
        temp(:,ny,:) = temp(:,nym,:)            ! Adiabatic
        co2(:,ny,:) = co2(:,nym,:)              ! Zero flux
        h2o(:,ny,:) = h2o(:,nym,:)              ! Zero flux
    end if

    return
    
end subroutine SetWallBCs

subroutine SetInletBC

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays
    use mpih

    implicit none

    integer :: k,j
    real    :: cell_ratio,ivel_slow

    ! This ensures that the inlet velocity slowly increases from zero using a smooth third order spline.
    if (time.le.tvel) then
        ivel_slow = ivel*((3.0d0*((time/tvel)**2)) - (2.0d0*((time/tvel)**3)))
    else 
        ivel_slow = ivel
    end if

    ! Wall with inlet vent
    if (xstart(3).eq.1) then
        vx(:,:,0) = -vx(:,:,1)                      ! Dirchlet condition
        vy(:,:,0) = -vy(:,:,1)                      ! Dirchlet condition
        vz(:,:,1) = 0.0d0                           ! Dirchlet condition
        temp(:,:,0) = temp(:,:,1)                   ! Adiabatic
        co2(:,:,0) = co2(:,:,1)                     ! Zero flux
        h2o(:,:,0) = h2o(:,:,1)                     ! Zero flux
    
        do k=1,nxm
            !Checking for all cells overlapping with the inlet
            if ((xc(k+1).gt.(iheight-(0.5d0*ilen))).and.(xc(k).lt.(iheight+(0.5d0*ilen)))) then
                cell_ratio = (min(xc(k+1),(iheight+(0.5d0*ilen))) - max(xc(k),(iheight-(0.5d0*ilen))))*dx/dx3c(k)
                vz(k,:,1) = ivel_slow*cell_ratio    ! Dirchlet condition
            end if
        end do
    end if

end subroutine SetInletBC

subroutine SetOutletBC

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vx,vy,vz,temp,co2,h2o
    use mpih

    implicit none

    integer :: i,j,k
    
    ! In the slab code developed by ChongShen Ng, Steven Chong, Rui Yang and Naoki Hori, the outlet
    ! boundary condition is an radiative/advective/non-reflecting boundary condition. The idea is to
    ! have a speed that is somehow related to the Courant number and determines the speed of 
    ! the Courant waves. Any disturbance at the outlet which travels in waves slower than this speed 
    ! gets advected out. They solve d(q)/d(t) + C d(q)/d(z) = 0 at the outlet.
    !
    ! We tried this boundary condition and many more including imposing a fixed velocity
    ! at the outlet. This always caused a numerical instability in the form waves entering the domain.
    ! We found through trial and error that this is an instability related to the diffusive terms
    ! since the diffusive terms in Z and Y direction are handled explicitly as compared to the slab code
    ! where they are handled implicitly. Through trial and error, we found that imposing the double derivative
    ! normal to the outflow as zero solves this issue (i.e. d^2(q)/dz^2 = 0). In a way, the solution at outflow 
    ! is linearly interpolated from the interior. 
    !
    ! We still do not fully understand why it solves the issue though. However, the outflow is always going 
    ! to be aphysical in nature and we do not care too much about the physical accuracy of the solution 
    ! close to the outlet for the expected application of this code.
    ! 
    ! - Vanshu and Chris

    ! Wall with outlet vent
    if (xend(3).eq.nzm) then
        vx(:,:,nz) = -vx(:,:,nzm)               ! Dirchlet condition
        vy(:,:,nz) = -vy(:,:,nzm)               ! Dirchlet condition
        vz(:,:,nz) = 0.0d0                      ! Dirchlet condition
        temp(:,:,nz) = temp(:,:,nzm)            ! Adiabatic 
        co2(:,:,nz) = co2(:,:,nzm)              ! Zero flux 
        h2o(:,:,nz) = h2o(:,:,nzm)              ! Zero flux 

        do k=1,nx
            ! Check if the node lies within the outlet dimensions
            if ((xc(k).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
                vx(k,:,nz)      = 2.0d0*vx(k,:,nzm)     - vx(k,:,nzm-1)
                temp(k,:,nz)    = 2.0d0*temp(k,:,nzm)   - temp(k,:,nzm-1)
                co2(k,:,nz)     = 2.0d0*co2(k,:,nzm)    - co2(k,:,nzm-1)
                h2o(k,:,nz)     = 2.0d0*h2o(k,:,nzm)    - h2o(k,:,nzm-1)
            end if
        end do

        do k=1,nxm
            ! Check for any overlap of grid cell with outlet
            if ((xc(k+1).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
                vy(k,:,nz)      = 2.0d0*vy(k,:,nzm)     - vy(k,:,nzm-1)
                vz(k,:,nz)      = 2.0d0*vz(k,:,nzm)     - vz(k,:,nzm-1)
            end if
        end do

    end if
            
end subroutine SetOutletBC

subroutine CorrectOutletFlux

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays
    use mpih

    implicit none

    integer :: j,k
    real    :: iflux,iarea,oflux,oarea
    real    :: cell_ratio,res_dummy

    !Get the outlet flux

    iarea = 0.0d0
    iflux = 0.0d0
    if (xstart(3).eq.1) then
        do k=1,nxm
            do j=xstart(2),xend(2)
                if ((xc(k+1).gt.(iheight-(0.5d0*ilen))).and.(xc(k).lt.(iheight+(0.5d0*ilen)))) then
                    iflux = iflux + (vz(k,j,1)*dx3c(k)/(dx*dy))         ! Compute vz inlet flux
                    iarea = iarea + (dx3c(k)/(dx*dy))                   ! Compute inlet area
                end if
            end do
        end do
    end if

    call MpiAllSumRealScalar(iflux,res_dummy)
    iflux = res_dummy

    call MpiAllSumRealScalar(iarea,res_dummy)
    iarea = res_dummy

    oarea = 0.0d0
    oflux = 0.0d0
    if (xend(3).eq.nzm) then
        do k=1,nxm
            do j=xstart(2),xend(2)
                if ((xc(k+1).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
                    oflux = oflux + (vz(k,j,nz)*dx3c(k)/(dx*dy))        ! Compute vz outlet flux
                    oarea = oarea + (dx3c(k)/(dx*dy))                   ! Compute outlet area
                end if
            end do
        end do
    end if

    call MpiAllSumRealScalar(oflux,res_dummy)
    oflux = res_dummy

    call MpiAllSumRealScalar(oarea,res_dummy)
    oarea = res_dummy

    ! Correct the outlet flux
    if (xend(3).eq.nzm) then
        do k=1,nxm
            ! if ((xm(k).gt.(oheight-(0.5d0*olen))).and.(xm(k).lt.(oheight+(0.5d0*olen)))) then
            if ((xc(k+1).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
                vz(k,:,nz) = vz(k,:,nz) + ((iflux - oflux)/oarea)    ! Correct vz outlet flux
                ! vz(k,:,nz) = (iflux)/oarea    ! Correct vz outlet flux
            end if
        end do
    end if

    ! DEBUG PURPOSES ==============================================================================
    ! if (ismaster) then
    !     open(100,file='Results/debug_flux.out',status='unknown',position='append',access='sequential')
    !     if ((ntime.eq.0).and.(.not. readflow)) then
    !         write(100,'(7(A14,X))') &
    !         'Time','Step','InFlux','InArea','OutFlux','OutArea','Correction'
    !     end if
    !     write(100,'((E14.6,X),(I14,X),5(E14.6,X))') &
    !     time,ntime,iflux,iarea,oflux,oarea,((iflux - oflux)/oarea)
    !     close(100)
    ! end if
    ! DEBUG PURPOSES ==============================================================================

end subroutine CorrectOutletFlux

subroutine SetPressureBC

    ! This subroutine is required to ensure that the outlet B.C.s are not
    ! replaced after a halo copy at the end of the time marcher scheme

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: dphhalo
    use mpih

    implicit none

    if (xstart(2).eq.1) then
        dphhalo(:,0,:)=dphhalo(:,1,:)
    end if

    if (xstart(3).eq.1) then
        dphhalo(:,:,0)=dphhalo(:,:,1)
    end if

    if (xend(2).eq.nym) then
        dphhalo(:,ny,:)=dphhalo(:,nym,:)
    end if

    if (xend(3).eq.nzm) then
        dphhalo(:,:,nz)=dphhalo(:,:,nzm)
    end if

end subroutine SetPressureBC

! ********** DEBUG ROUTINES **********
subroutine SetDebugWallBCs

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays
    use mpih

    implicit none

    integer :: i,j,k

    vx(1,:,:) = 0.0d0
    vx(nx,:,:) = 0.0d0

    ! Just plain side walls, no vents
    if (xstart(2).eq.1) then
        vx(:,0,:) = -vx(:,1,:)                  ! Dirchlet condition
        vy(:,1,:) = 0.0d0                       ! Dirchlet condition
        vz(:,0,:) = -vz(:,1,:)                  ! Dirchlet condition
        temp(:,0,:) = temp(:,1,:)               ! Adiabatic
        co2(:,0,:) = co2(:,1,:)                 ! Zero flux
        h2o(:,0,:) = h2o(:,1,:)                 ! Zero flux
    end if

    if (xend(2).eq.nym) then
        vx(:,ny,:) = -vx(:,nym,:)               ! Dirchlet condition
        vy(:,ny,:) = 0.0d0                       ! Dirchlet condition
        vz(:,ny,:) = -vz(:,nym,:)               ! Dirchlet condition
        temp(:,ny,:) = temp(:,nym,:)            ! Adiabatic
        co2(:,ny,:) = co2(:,nym,:)              ! Zero flux
        h2o(:,ny,:) = h2o(:,nym,:)              ! Zero flux
    end if

    ! Wall with inlet vent
    if (xstart(3).eq.1) then
        vx(:,:,0) = -vx(:,:,1)                  ! Dirchlet condition
        vy(:,:,0) = -vy(:,:,1)                  ! Dirchlet condition
        vz(:,:,1) = 0.0d0                       ! Dirchlet condition
        temp(:,:,0) = temp(:,:,1)               ! Adiabatic
        co2(:,:,0) = co2(:,:,1)                 ! Zero flux
        h2o(:,:,0) = h2o(:,:,1)                 ! Zero flux
    end if

    ! Wall with outlet vent
    if (xend(3).eq.nzm) then
        vx(:,:,nz) = -vx(:,:,nzm)               ! Dirchlet condition
        vy(:,:,nz) = -vy(:,:,nzm)               ! Dirchlet condition
        vz(:,:,nz) = 0.0d0                      ! Dirchlet condition
        temp(:,:,nz) = temp(:,:,nzm)            ! Adiabatic 
        co2(:,:,nz) = co2(:,:,nzm)              ! Zero flux 
        h2o(:,:,nz) = h2o(:,:,nzm)              ! Zero flux 
    end if

    return
    
end subroutine SetDebugWallBCs
! ********** DEBUG ROUTINES **********