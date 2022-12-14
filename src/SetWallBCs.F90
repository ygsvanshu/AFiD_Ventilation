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

subroutine CalcOutletBC_RK3

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vx,vy,vz,temp,co2,h2o
    use ventilation_arrays
    use mpih

    implicit none

    integer :: i,j,k
    real    :: cou,dqdz,cll,crr

    cll = 0.5d0*real(odesc)
    crr = 0.5d0*real(2-odesc)

    ! Apparently this is somehow related to the Courant number and determines the speed of 
    ! the Courant waves. Not sure why this is a fixed number instead of parameter*(delta_z/delta_t)
    ! Any disturbance at the outlet which travels in waves slower than this speed gets advected out
    ! Also called radiative/advective/non-reflecting boundary condition
    cou = 0.3

    ! Wall with outlet vent
    if (xend(3).eq.nzm) then
        do k=1,nx
            ! Check if the node lies within the outlet dimensions
            if ((xc(k).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
                ! Handy hack to avoid halo update for out quantities
                do j=xstart(2)-lvlhalo,xend(2)+lvlhalo
                    ! Radiative outflow boundary condition for vx
                    dqdz = (vx(k,j,nz) - vx(k,j,nzm))*dz
                    outvx(k,j) = ((crr*vx(k,j,nz)) + (cll*vx(k,j,nzm)))
                    outvx(k,j) = outvx(k,j) - (ga*dqdz+ro*dzoutvx(k,j))*cou*al*dt 
                    dzoutvx(k,j) = dqdz
                    ! Radiative outflow boundary condition for temp
                    dqdz = (temp(k,j,nz) - temp(k,j,nzm))*dz
                    outtemp(k,j) = ((crr*temp(k,j,nz)) + (cll*temp(k,j,nzm)))
                    outtemp(k,j) = outtemp(k,j) - (ga*dqdz+ro*dzouttemp(k,j))*cou*al*dt
                    dzouttemp(k,j) = dqdz
                    ! Radiative outflow boundary condition for co2
                    dqdz = (co2(k,j,nz) - co2(k,j,nzm))*dz
                    outco2(k,j) = ((crr*co2(k,j,nz)) + (cll*co2(k,j,nzm)))
                    outco2(k,j) = outco2(k,j) - (ga*dqdz+ro*dzoutco2(k,j))*cou*al*dt
                    dzoutco2(k,j) = dqdz
                    ! Radiative outflow boundary condition for h2o
                    dqdz = (h2o(k,j,nz) - h2o(k,j,nzm))*dz
                    outh2o(k,j) = ((crr*h2o(k,j,nz)) + (cll*h2o(k,j,nzm)))
                    outh2o(k,j) = outh2o(k,j) - (ga*dqdz+ro*dzouth2o(k,j))*cou*al*dt
                    dzouth2o(k,j) = dqdz
                end do
            end if
        end do

        do k=1,nxm
            ! Check for any overlap of grid cell with outlet
            ! if ((xm(k).gt.(oheight-(0.5d0*olen))).and.(xm(k).lt.(oheight+(0.5d0*olen)))) then
            if ((xc(k+1).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
                do j=xstart(2)-lvlhalo,xend(2)+lvlhalo
                    ! Radiative outflow boundary condition for vy
                    dqdz = (vy(k,j,nz) - vy(k,j,nzm))*dz
                    outvy(k,j) = ((crr*vy(k,j,nz)) + (cll*vy(k,j,nzm)))
                    outvy(k,j) = outvy(k,j) - (ga*dqdz+ro*dzoutvy(k,j))*cou*al*dt
                    dzoutvy(k,j) = dqdz
                    ! Radiative outflow boundary condition for vz
                    ! This is different because it lies on the outlet instead of half a grid cell away
                    ! This is a bit inaccurate, but this is the best we can do :'(   -[GSY]
                    ! If this is not "upwind" then the outlet is not stable and solution blows up
                    dqdz = (vz(k,j,nz) - vz(k,j,nzm))*dz 
                    outvz(k,j) = vz(k,j,nz)
                    outvz(k,j) = outvz(k,j) - (ga*dqdz+ro*dzoutvz(k,j))*cou*al*dt
                    dzoutvz(k,j) = dqdz
                end do
            end if
        end do
    end if

    return
    
end subroutine CalcOutletBC_RK3

subroutine SetInletBC

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays
    use ventilation_arrays
    use mpih

    implicit none

    integer :: k,j
    real    :: cell_ratio

    ! Wall with inlet vent
    if (xstart(3).eq.1) then
        vx(:,:,0) = -vx(:,:,1)                  ! Dirchlet condition
        vy(:,:,0) = -vy(:,:,1)                  ! Dirchlet condition
        vz(:,:,1) = 0.0d0                       ! Dirchlet condition
        temp(:,:,0) = temp(:,:,1)               ! Adiabatic
        co2(:,:,0) = co2(:,:,1)                 ! Zero flux
        h2o(:,:,0) = h2o(:,:,1)                 ! Zero flux
    
        do k=1,nxm
            !Checking for all cells overlapping with the inlet
            if ((xc(k+1).gt.(iheight-(0.5d0*ilen))).and.(xc(k).lt.(iheight+(0.5d0*ilen)))) then
                cell_ratio = (min(xc(k+1),(iheight+(0.5d0*ilen))) - max(xc(k),(iheight-(0.5d0*ilen))))*dx/dx3c(k)
                vz(k,:,1) = ivel*cell_ratio! Dirchlet condition
            end if
        end do

        ! ********** DEBUG ********** !
        do k=1,nx
            if ((xc(k).gt.(iheight-(0.5d0*ilen))).and.(xc(k).lt.(iheight+(0.5d0*ilen)))) then
                temp(k,:,0) = -2.0d0 - temp(k,:,1)! Dirchlet condition
            end if
        end do
        ! ********** DEBUG ********** !

    end if

end subroutine SetInletBC
    
subroutine SetOutletBC_RK3

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vx,vy,vz,temp,co2,h2o
    use ventilation_arrays
    use mpih

    implicit none

    integer :: i,j,k
    real    :: inlet_flux, inlet_area, outlet_flux, outlet_area
    real    :: res_dummy,c11,c22

    c11 = real(1+odesc)
    c22 = real(odesc)

    ! Wall with outlet vent
    if (xend(3).eq.nzm) then
        do k=1,nx
            ! Check if the node lies within the outlet dimensions
            if ((xc(k).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
                ! Apply outflow boundary condition for vx
                vx(k,:,nz) = c11*outvx(k,:) - c22*vx(k,:,nzm)
                ! Apply outflow boundary condition for temp
                temp(k,:,nz) = c11*outtemp(k,:) - c22*temp(k,:,nzm)
                ! Apply outflow boundary condition for co2
                co2(k,:,nz) = c11*outco2(k,:) - c22*co2(k,:,nzm)
                ! Apply outflow boundary condition for h20
                h2o(k,:,nz) = c11*outh2o(k,:) - c22*h2o(k,:,nzm)
            else
                ! Apply dirchlet boundary condition for vx on wall
                vx(k,:,nz) = -vx(k,:,nzm)
                ! Apply adiabatic boundary condition for temp on wall
                temp(k,:,nz) = temp(k,:,nzm)
                ! Apply zero flux boundary condition for co2 on wall
                co2(k,:,nz) = co2(k,:,nzm)
                ! Apply zero flux boundary condition for h2o on wall
                h2o(k,:,nz) = h2o(k,:,nzm)
            end if
        end do

        do k=1,nxm
            ! Check for any overlap of grid cell with outlet
            ! if ((xm(k).gt.(oheight-(0.5d0*olen))).and.(xm(k).lt.(oheight+(0.5d0*olen)))) then
            if ((xc(k+1).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
                ! Apply outflow boundary condition for vy
                vy(k,:,nz) = c11*outvy(k,:) - c22*vy(k,:,nzm)
                ! Apply outflow boundary condition for vz
                ! This is different because it lies directly on the outlet instead of half a grid cell away
                vz(k,:,nz) = outvz(k,:)
                ! vz(k,:,nz) = vz(k,:,nzm)
            else
                ! Apply dirchlet boundary condition for vy on wall
                vy(k,:,nz) = -vy(k,:,nzm)
                ! Apply dirchlet boundary condition for vz on wall
                vz(k,:,nz) = 0.0d0
            end if
        end do

    end if
            
end subroutine SetOutletBC_RK3

subroutine CorrectOutletFlux

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays
    use ventilation_arrays
    use mpih

    implicit none

    integer :: j,k
    real    :: inlet_flux,inlet_area,outlet_flux,outlet_area
    real    :: cell_ratio,res_dummy

    !Get the outlet flux
    inlet_area  = ilen*ylen
    inlet_flux  = ivel*inlet_area
    outlet_area = 0.0d0
    outlet_flux = 0.0d0
    if (xend(3).eq.nzm) then
        do k=1,nxm
            do j=xstart(2),xend(2)
                ! if ((xm(k).gt.(oheight-(0.5d0*olen))).and.(xm(k).lt.(oheight+(0.5d0*olen)))) then
                if ((xc(k+1).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
                    outlet_flux = outlet_flux + (vz(k,j,nz)*dx3c(k)/(dx*dy))        ! Compute vz outlet flux
                    outlet_area = outlet_area + (dx3c(k)/(dx*dy))                   ! Compute outlet area
                end if
            end do
        end do
    end if

    call MpiAllSumRealScalar(outlet_flux,res_dummy)
    outlet_flux = res_dummy

    call MpiAllSumRealScalar(outlet_area,res_dummy)
    outlet_area = res_dummy

    ! Correct the outlet flux
    if (xend(3).eq.nzm) then
        do k=1,nxm
            ! if ((xm(k).gt.(oheight-(0.5d0*olen))).and.(xm(k).lt.(oheight+(0.5d0*olen)))) then
            if ((xc(k+1).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
                vz(k,:,nz) = vz(k,:,nz) + ((inlet_flux - outlet_flux)/outlet_area)    ! Correct vz outlet flux
                ! vz(k,:,nz) = (inlet_flux)/outlet_area    ! Correct vz outlet flux
            end if
        end do
    end if

    ! DEBUG PURPOSES ==============================================================================
    if (ismaster) then
        open(100,file='Results/debug_flux.out',status='unknown',position='append',access='sequential')
        if ((ntime.eq.0).and.(.not. readflow)) then
            write(100,'(7(A14,X))') &
            'Time','Step','InFlux','InArea','OutFlux','OutArea','Correction'
        end if
        write(100,'((E14.6,X),(I14,X),5(E14.6,X))') &
        time,ntime,inlet_flux,inlet_area,outlet_flux,outlet_area,((inlet_flux - outlet_flux)/outlet_area)
        close(100)
    end if
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

subroutine CalcOutletBC_CKN

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vx,vy,vz,temp,co2,h2o
    use ventilation_arrays
    use mpih

    implicit none

    integer :: i,j,k

    ! Wall with outlet vent
    if (xend(3).eq.nzm) then
        do k=1,nx
            ! Check if the node lies within the outlet dimensions
            if ((xc(k).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
                ! Handy hack to avoid halo update for out quantities
                ! Radiative outflow boundary condition for vx
                outvx(k,:) = (vx(k,:,nz) + vx(k,:,nzm))*0.5d0
                dzoutvx(k,:) = (vx(k,:,nz) - vx(k,:,nzm))*dz
                ! Radiative outflow boundary condition for temp
                outtemp(k,:) = (temp(k,:,nz) + temp(k,:,nzm))*0.5d0
                dzouttemp(k,:) = (temp(k,:,nz) - temp(k,:,nzm))*dz
                ! Radiative outflow boundary condition for co2
                outco2(k,:) = (co2(k,:,nz) + co2(k,:,nzm))*0.5d0
                dzoutco2(k,:) = (co2(k,:,nz) - co2(k,:,nzm))*dz
                ! Radiative outflow boundary condition for h2o
                outh2o(k,:) = (h2o(k,:,nz) + h2o(k,:,nzm))*0.5d0
                dzouth2o(k,:) = (h2o(k,:,nz) - h2o(k,:,nzm))*dz
            end if
        end do

        do k=1,nxm
            ! Check for any overlap of grid cell with outlet
            if ((xc(k+1).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
                ! Radiative outflow boundary condition for vy
                outvy(k,:) = (vy(k,:,nz) + vy(k,:,nzm))*0.5d0
                dzoutvy(k,:) = (vy(k,:,nz) - vy(k,:,nzm))*dz
                ! Radiative outflow boundary condition for vz
                ! This is different because it lies on the outlet instead of half a grid cell away
                outvz(k,:) = vz(k,:,nz)
                dzoutvz(k,:) = (vz(k,:,nz) - vz(k,:,nzm))*dz ! This is a bit inaccurate, but this is the best we can do :'(   -[GSY]
            end if
        end do
    end if

    return
    
end subroutine CalcOutletBC_CKN

subroutine SetOutletBC_CKN

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vx,vy,vz,temp,co2,h2o
    use ventilation_arrays
    use mpih

    implicit none

    integer :: i,j,k
    real    :: inlet_flux, inlet_area, outlet_flux, outlet_area
    real    :: cou,res_dummy

    ! Apparently this is somehow related to the Courant number and determines the speed of 
    ! the Courant waves. Not sure why this is a fixed number instead of parameter*(delta_z/delta_t)
    ! Any disturbance at the outlet which travels in waves slower than this speed gets advected out
    ! Also called radiative/advective/non-reflecting boundary condition
    cou = 0.3d0

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
                ! Apply outflow boundary condition for vx
                vx(k,:,nz) = (((2.0d0*outvx(k,:) - vx(k,:,nzm))/dt) + (cou*(vx(k,:,nzm)*dz - dzoutvx(k,:))))/(1.0d0/dt + cou*dz)
                ! Apply outflow boundary condition for temp
                temp(k,:,nz) = (((2.0d0*outtemp(k,:) - temp(k,:,nzm))/dt) + (cou*(temp(k,:,nzm)*dz - dzouttemp(k,:))))/(1.0d0/dt + cou*dz)
                ! Apply outflow boundary condition for co2
                co2(k,:,nz) = (((2.0d0*outco2(k,:) - co2(k,:,nzm))/dt) + (cou*(co2(k,:,nzm)*dz - dzoutco2(k,:))))/(1.0d0/dt + cou*dz)
                ! Apply outflow boundary condition for h20
                h2o(k,:,nz) = (((2.0d0*outh2o(k,:) - h2o(k,:,nzm))/dt) + (cou*(h2o(k,:,nzm)*dz - dzouth2o(k,:))))/(1.0d0/dt + cou*dz)
            end if
        end do

        do k=1,nxm
            ! Check for any overlap of grid cell with outlet
            if ((xc(k+1).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
                ! Apply outflow boundary condition for vy
                vy(k,:,nz) = (((2.0d0*outvy(k,:) - vy(k,:,nzm))/dt) + (cou*(vy(k,:,nzm)*dz - dzoutvy(k,:))))/(1.0d0/dt + cou*dz)
                ! Apply outflow boundary condition for vz
                ! This is different because it lies directly on the outlet instead of half a grid cell away
                vz(k,:,nz) = ((2.0d0*outvz(k,:)/dt) + (cou*(vz(k,:,nzm)*dz - dzoutvz(k,:))))/(2.0d0/dt + cou*dz)
            end if
        end do

    end if
            
end subroutine SetOutletBC_CKN

subroutine CalcOutletBC_CKN_UPWIND

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vx,vy,vz,temp,co2,h2o
    use ventilation_arrays
    use mpih

    implicit none

    integer :: i,j,k

    ! Wall with outlet vent
    if (xend(3).eq.nzm) then
        do k=1,nx
            ! Check if the node lies within the outlet dimensions
            if ((xc(k).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
                ! Handy hack to avoid halo update for out quantities
                ! Radiative outflow boundary condition for vx
                outvx(k,:) = vx(k,:,nz)
                dzoutvx(k,:) = (vx(k,:,nz) - vx(k,:,nzm))
                ! Radiative outflow boundary condition for temp
                outtemp(k,:) = temp(k,:,nz)
                dzouttemp(k,:) = (temp(k,:,nz) - temp(k,:,nzm))
                ! Radiative outflow boundary condition for co2
                outco2(k,:) = co2(k,:,nz)
                dzoutco2(k,:) = (co2(k,:,nz) - co2(k,:,nzm))
                ! Radiative outflow boundary condition for h2o
                outh2o(k,:) = h2o(k,:,nz)
                dzouth2o(k,:) = (h2o(k,:,nz) - h2o(k,:,nzm))
            end if
        end do

        do k=1,nxm
            ! Check for any overlap of grid cell with outlet
            if ((xc(k+1).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
                ! Radiative outflow boundary condition for vy
                outvy(k,:) = vy(k,:,nz)
                dzoutvy(k,:) = (vy(k,:,nz) - vy(k,:,nzm))
                ! Radiative outflow boundary condition for vz
                ! This is different because it lies on the outlet instead of half a grid cell away
                outvz(k,:) = vz(k,:,nz)
                dzoutvz(k,:) = (vz(k,:,nz) - vz(k,:,nzm)) 
                ! This is a bit inaccurate, but this is the best we can do :'(   -[GSY]
            end if
        end do
    end if

    return
    
end subroutine CalcOutletBC_CKN_UPWIND

subroutine SetOutletBC_CKN_UPWIND

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vx,vy,vz,temp,co2,h2o
    use ventilation_arrays
    use mpih

    implicit none

    integer :: i,j,k
    real    :: inlet_flux, inlet_area, outlet_flux, outlet_area
    real    :: cou,bet,res_dummy

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
                ! Apply outflow boundary condition for vx
                vx(k,:,nz) = (outvx(k,:) + (bet*(vx(k,:,nzm) - dzoutvx(k,:))))/(1.0d0 + bet)
                ! Apply outflow boundary condition for temp
                temp(k,:,nz) = (outtemp(k,:) + (bet*(temp(k,:,nzm) - dzouttemp(k,:))))/(1.0d0 + bet)
                ! Apply outflow boundary condition for co2
                co2(k,:,nz) = (outco2(k,:) + (bet*(co2(k,:,nzm) - dzoutco2(k,:))))/(1.0d0 + bet)
                ! Apply outflow boundary condition for h20
                h2o(k,:,nz) = (outh2o(k,:) + (bet*(h2o(k,:,nzm) - dzouth2o(k,:))))/(1.0d0 + bet)
            end if
        end do

        do k=1,nxm
            ! Check for any overlap of grid cell with outlet
            if ((xc(k+1).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
                ! Apply outflow boundary condition for vy
                vy(k,:,nz) = (outvy(k,:) + (bet*(vy(k,:,nzm) - dzoutvy(k,:))))/(1.0d0 + bet)
                ! Apply outflow boundary condition for vz
                ! This is different because it lies directly on the outlet instead of half a grid cell away
                vz(k,:,nz) = (outvz(k,:) + (bet*(vz(k,:,nzm) - dzoutvz(k,:))))/(1.0d0 + bet)
            end if
        end do

    end if
            
end subroutine SetOutletBC_CKN_UPWIND
! ********** DEBUG ROUTINES **********