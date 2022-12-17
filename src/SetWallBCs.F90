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
    use ventilation_arrays
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
    use ventilation_arrays
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

subroutine CalcOutletBC

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vx,vy,vz,temp,co2,h2o
    use ventilation_arrays
    use mpih

    implicit none

    integer :: i,j,k

    ! This is the first part of the upwind Crank Nicolson based update of outlet nodes
    ! Here we store the derivative for the current timestep at the wall with outlet vent
    if (xend(3).eq.nzm) then
        do k=1,nx
            ! Check if the node lies within the outlet dimensions
            if ((xc(k).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
                ! Handy hack to avoid halo update for out quantities
                dzoutvx(k,:) = (vx(k,:,nz) - vx(k,:,nzm))
                dzouttemp(k,:) = (temp(k,:,nz) - temp(k,:,nzm))
                dzoutco2(k,:) = (co2(k,:,nz) - co2(k,:,nzm))
                dzouth2o(k,:) = (h2o(k,:,nz) - h2o(k,:,nzm))
            end if
        end do

        do k=1,nxm
            ! Check for any overlap of grid cell with outlet
            if ((xc(k+1).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
                ! Handy hack to avoid halo update for out quantities
                dzoutvy(k,:) = (vy(k,:,nz) - vy(k,:,nzm))
                dzoutvz(k,:) = (vz(k,:,nz) - vz(k,:,nzm)) 
            end if
        end do
    end if

    return
    
end subroutine CalcOutletBC

subroutine SetOutletBC

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vx,vy,vz,temp,co2,h2o
    use ventilation_arrays
    use mpih

    implicit none

    integer :: i,j,k
    real    :: cou,bet,res_dummy

    ! This is the second part of the upwind Crank Nicolson based update of outlet nodes
    
    ! Apparently this is somehow related to the Courant number and determines the speed of 
    ! the Courant waves. Not sure why this is a fixed number instead of parameter*(delta_z/delta_t)
    ! Any disturbance at the outlet which travels in waves slower than this speed gets advected out
    ! Also called radiative/advective/non-reflecting boundary condition
    cou = limitCFL
    bet = 0.5*cou*al*dt*dz

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
                vx(k,:,nz) = (vx(k,:,nz) + (bet*(vx(k,:,nzm) - dzoutvx(k,:))))/(1.0d0 + bet)
                temp(k,:,nz) = (temp(k,:,nz) + (bet*(temp(k,:,nzm) - dzouttemp(k,:))))/(1.0d0 + bet)
                co2(k,:,nz) = (co2(k,:,nz) + (bet*(co2(k,:,nzm) - dzoutco2(k,:))))/(1.0d0 + bet)
                h2o(k,:,nz) = (h2o(k,:,nz) + (bet*(h2o(k,:,nzm) - dzouth2o(k,:))))/(1.0d0 + bet)
            end if
        end do

        do k=1,nxm
            ! Check for any overlap of grid cell with outlet
            if ((xc(k+1).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
                vy(k,:,nz) = (vy(k,:,nz) + (bet*(vy(k,:,nzm) - dzoutvy(k,:))))/(1.0d0 + bet)
                vz(k,:,nz) = (vz(k,:,nz) + (bet*(vz(k,:,nzm) - dzoutvz(k,:))))/(1.0d0 + bet)
            end if
        end do

    end if
            
end subroutine SetOutletBC

subroutine CorrectOutletFlux

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays
    use ventilation_arrays
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

subroutine CopyOutletBC

    ! This subroutine is required to ensure that the outlet B.C.s are not
    ! replaced after a halo copy at the end of the time marcher scheme

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays
    use ventilation_arrays
    use mpih

    implicit none

    if (xend(3).eq.nzm) then
        outvx(:,:)   = vx(:,:,nz)
        outvy(:,:)   = vy(:,:,nz)
        outvz(:,:)   = vz(:,:,nz)
        outtemp(:,:) = temp(:,:,nz)
        outco2(:,:)  = co2(:,:,nz)
        outh2o(:,:)  = h2o(:,:,nz)
    end if

end subroutine CopyOutletBC

subroutine PasteOutletBC

    ! This subroutine is required to ensure that the outlet B.C.s are not
    ! replaced after a halo copy at the end of the time marcher scheme

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays
    use ventilation_arrays
    use mpih

    implicit none

    if (xend(3).eq.nzm) then
        vx(:,:,nz)   = outvx(:,:)
        vy(:,:,nz)   = outvy(:,:)
        vz(:,:,nz)   = outvz(:,:)
        temp(:,:,nz) = outtemp(:,:)
        co2(:,:,nz)  = outco2(:,:)
        h2o(:,:,nz)  = outh2o(:,:)
    end if

end subroutine PasteOutletBC

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
        vy(:,1,:) = 0.0d0                        ! Dirchlet condition
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