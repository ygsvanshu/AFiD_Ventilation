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

subroutine InitVents

    use param
    use vent_arrays
    use mpih

    implicit none

    integer :: k

    icell(:) = 0.0d0
    ocell(:) = 0.0d0

    iarea = ilen*ylen
    oarea = olen*ylen

    ixmst = nxm
    ixmen = 1
    oxmst = nxm
    oxmen = 1

    ixfst = nxm
    ixfen = 1
    oxfst = nxm
    oxfen = 1

    do k=1,nxm
        !Checking for all grid points within the inlet
        if ((xm(k).gt.(iheight-(0.5d0*ilen))).and.(xm(k).lt.(iheight+(0.5d0*ilen)))) then
            ixmst       = min(ixmst,k)
            ixmen       = max(ixmen,k)
        end if
        !Checking for all cells overlapping with the inlet
        if ((xc(k+1).gt.(iheight-(0.5d0*ilen))).and.(xc(k).lt.(iheight+(0.5d0*ilen)))) then
            ixfst       = min(ixfst,k)
            ixfen       = max(ixfen,k)
            icell(k)    = (min(xc(k+1),(iheight+(0.5d0*ilen))) - max(xc(k),(iheight-(0.5d0*ilen))))*dx/dx3c(k)
        end if
        !Checking for all grid points within the outlet
        if ((xm(k).gt.(oheight-(0.5d0*olen))).and.(xm(k).lt.(oheight+(0.5d0*olen)))) then
            oxmst       = min(oxmst,k)
            oxmen       = max(oxmen,k)
        end if
        !Checking for all cells overlapping with the outlet
        if ((xc(k+1).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
            oxfst       = min(oxfst,k)
            oxfen       = max(oxfen,k)
            ocell(k)    = (min(xc(k+1),(oheight+(0.5d0*olen))) - max(xc(k),(oheight-(0.5d0*olen))))*dx/dx3c(k)
        end if
    end do

    ixcst = nx
    ixcen = 1
    oxcst = nx
    oxcen = 1

    do k=1,nx
        !Checking for all grid points within the inlet
        if ((xc(k).gt.(iheight-(0.5d0*ilen))).and.(xc(k).lt.(iheight+(0.5d0*ilen)))) then
            ixcst       = min(ixcst,k)
            ixcen       = max(ixcen,k)
        end if
        !Checking for all grid points within the outlet
        if ((xc(k).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
            oxcst       = min(oxcst,k)
            oxcen       = max(oxcen,k)
        end if
    end do

end subroutine InitVents

subroutine SetWallBCs

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vx,vy,vz,temp,co2,h2o
    use vent_arrays
    use mpih

    implicit none

    ! Top and bottom walls (X=0,X=ALX3D)
    vx(1,:,:) = 0.d0
    vx(nx,:,:) = 0.d0

    ! Side walls without vents (Y=0,Y=YLEN)
    if (xstart(2).eq.1) then
        vx(:,0,:)   = -vx(:,1,:)                    ! Dirchlet condition
        vy(:,1,:)   = 0.d0                          ! Dirchlet condition
        vz(:,0,:)   = -vz(:,1,:)                    ! Dirchlet condition
        temp(:,0,:) = temp(:,1,:)                   ! Adiabatic
        co2(:,0,:)  = co2(:,1,:)                    ! Zero flux
        h2o(:,0,:)  = h2o(:,1,:)                    ! Zero flux
    end if

    if (xend(2).eq.nym) then
        vx(:,ny,:)   = -vx(:,nym,:)                 ! Dirchlet condition
        vy(:,ny,:)   = 0.d0                         ! Dirchlet condition
        vz(:,ny,:)   = -vz(:,nym,:)                 ! Dirchlet condition
        temp(:,ny,:) = temp(:,nym,:)                ! Adiabatic
        co2(:,ny,:)  = co2(:,nym,:)                 ! Zero flux
        h2o(:,ny,:)  = h2o(:,nym,:)                 ! Zero flux
    end if

    ! Front and back walls with vents (Z=0,Z=ZLEN)
    if (xstart(3).eq.1) then
        vx(:,:,0) = -vx(:,:,1)                      ! Dirchlet condition
        vy(:,:,0) = -vy(:,:,1)                      ! Dirchlet condition
        vz(:ixfst-1,:,1) = 0.0d0                    ! Dirchlet condition
        vz(ixfen+1:,:,1) = 0.0d0                    ! Dirchlet condition
        temp(:,:,0) = temp(:,:,1)                   ! Adiabatic
        co2(:,:,0) = co2(:,:,1)                     ! Zero flux
        h2o(:,:,0) = h2o(:,:,1)                     ! Zero flux
    end if

    if (xend(3).eq.nzm) then
        vx(:oxcst-1,:,nz)   = -vx(:oxcst-1,:,nzm)   ! Dirchlet condition
        vy(:oxfst-1,:,nz)   = -vy(:oxfst-1,:,nzm)   ! Dirchlet condition
        vz(:oxfst-1,:,nz)   = 0.0d0                 ! Dirchlet condition
        temp(:oxcst-1,:,nz) = temp(:oxcst-1,:,nzm)  ! Adiabatic
        co2(:oxcst-1,:,nz)  = co2(:oxcst-1,:,nzm)   ! Zero flux
        h2o(:oxcst-1,:,nz)  = h2o(:oxcst-1,:,nzm)   ! Zero flux

        vx(oxcst+1:,:,nz)   = -vx(oxcst+1:,:,nzm)   ! Dirchlet condition
        vy(oxfst+1:,:,nz)   = -vy(oxfst+1:,:,nzm)   ! Dirchlet condition
        vz(oxfst+1:,:,nz)   = 0.0d0                 ! Dirchlet condition
        temp(oxcst+1:,:,nz) = temp(oxcst+1:,:,nzm)  ! Adiabatic
        co2(oxcst+1:,:,nz)  = co2(oxcst+1:,:,nzm)   ! Zero flux
        h2o(oxcst+1:,:,nz)  = h2o(oxcst+1:,:,nzm)   ! Zero flux
    end if

    return
    
end subroutine SetWallBCs

subroutine CalcOutletBC

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vx,vy,vz,temp,co2,h2o
    use vent_arrays
    use mpih

    implicit none

    integer :: i,j,k

    ! This is the first part of the upwind Crank Nicolson based update of outlet nodes
    ! Here we store the derivative for the current timestep at the wall with outlet vent
    if (xend(3).eq.nzm) then
        do k=oxcst,oxcen
            outvx(k,:) = (vx(k,:,nz) - vx(k,:,nzm))
            outtemp(k,:) = (temp(k,:,nz) - temp(k,:,nzm))
            outco2(k,:) = (co2(k,:,nz) - co2(k,:,nzm))
            outh2o(k,:) = (h2o(k,:,nz) - h2o(k,:,nzm))
        end do
        do k=oxfst,oxfen
            outvy(k,:) = (vy(k,:,nz) - vy(k,:,nzm))
            outvz(k,:) = (vz(k,:,nz) - vz(k,:,nzm)) 
        end do
    end if

    return
    
end subroutine CalcOutletBC

subroutine SetOutletBC

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vx,vy,vz,temp,co2,h2o
    use vent_arrays
    use mpih

    implicit none

    integer :: i,j,k
    real    :: wavefactor

    ! This is the second part of the upwind Crank Nicolson based update of outlet nodes
    
    ! Apparently this is somehow related to the Courant number and determines the speed of 
    ! the Courant waves. Not sure why this is a fixed number instead of parameter*(delta_z/delta_t)
    ! Any disturbance at the outlet which travels in waves slower than this speed gets advected out
    ! Also called radiative/advective/non-reflecting boundary condition
    wavefactor = 0.5*wavespeed*al*dt*dz

    ! Wall with outlet vent
    if (xend(3).eq.nzm) then
        do k=oxcst,oxcen
            vx(k,:,nz) = (vx(k,:,nz) + (wavefactor*(vx(k,:,nzm) - outvx(k,:))))/(1.0d0 + wavefactor)
            temp(k,:,nz) = (temp(k,:,nz) + (wavefactor*(temp(k,:,nzm) - outtemp(k,:))))/(1.0d0 + wavefactor)
            co2(k,:,nz) = (co2(k,:,nz) + (wavefactor*(co2(k,:,nzm) - outco2(k,:))))/(1.0d0 + wavefactor)
            h2o(k,:,nz) = (h2o(k,:,nz) + (wavefactor*(h2o(k,:,nzm) - outh2o(k,:))))/(1.0d0 + wavefactor)
        end do
        do k=oxfst,oxfen
            vy(k,:,nz) = (vy(k,:,nz) + (wavefactor*(vy(k,:,nzm) - outvy(k,:))))/(1.0d0 + wavefactor)
            vz(k,:,nz) = (vz(k,:,nz) + (wavefactor*(vz(k,:,nzm) - outvz(k,:))))/(1.0d0 + wavefactor)
        end do
    end if    
            
end subroutine SetOutletBC

subroutine AddBC_Vx

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vx,rhsx
    use mpih

    implicit none

    integer :: k,j,i

    do j=xstart(2),xend(2)
        do i=xstart(3),xend(3)
            rhsx(1,j,i)  = -vx(1,j,i)
            rhsx(nx,j,i) = -vx(nx,j,i)
        end do
    end do

    return

end subroutine AddBC_Vx

subroutine AddBC_Vy

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vy,rhsx
    use mpih

    implicit none

    integer :: k,j,i

    if (xstart(2).eq.1) then
        do i=xstart(3),xend(3)
            do k=1,nxm
                rhsx(k,1,i) = -vy(k,1,i)
            end do
        end do
    end if

    return

end subroutine AddBC_Vy

subroutine AddBC_Vz

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vz,rhsx
    use vent_arrays, only: isvel,tsvel,ixfst,ixfen
    use mpih

    implicit none

    integer :: k,j,i
    real    :: cvel,ctime

    ! This ensures that the inlet velocity slowly increases from zero using a smooth third order spline.
    ctime = time + al*dt
    if (ctime.lt.(tsvel+tvel)) then
        cvel = isvel + ((ivel-isvel)*((3.0d0*(((ctime-tsvel)/tvel)**2)) - (2.0d0*(((ctime-tsvel)/tvel)**3))))
    else 
        cvel = ivel
    end if

    if (xstart(3).eq.1) then
        do j=xstart(2),xend(2)
            do k=1,nxm
                if ((k.ge.ixfst).and.(k.le.ixfen)) then
                    rhsx(k,j,1) = cvel - vz(k,j,1)
                else
                    rhsx(k,j,1) = -vz(k,j,1)
                end if
            end do
        end do
    end if

    return

end subroutine AddBC_Vz

subroutine AddOutletBC(qua,vsc,kst,ken)

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: rhsx
    use vent_arrays
    use mpih

    implicit none

    real, intent(in), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: qua
    real, intent(in) :: vsc
    integer, intent(in) :: kst,ken
    real    :: betadx,wavefactor,numerator,denominator
    integer :: k,j,i

    ! FOR THE RADIATIVE B.C. AT THE OUTLET
    betadx      = 0.5d0*al*dt/vsc
    wavespeed   = 0.3
    wavefactor  = wavespeed*al*dt*dz
    numerator   = 0.5d0*wavefactor*ap1gi*betadx
    denominator = ((1.0d0 - ac1gi*betadx)*(1.0d0 - wavefactor)) + (ap1gi*wavefactor*betadx)
    am1oi       = (am1gi*betadx*(1-wavefactor))/denominator

    if (xend(3).eq.nzm) then
        do j=xstart(2),xend(2)
            do k=kst,ken
                rhsx(k,j,nzm) = (rhsx(k,j,nzm) + (numerator*(qua(k,j,nz)-qua(k,j,nzm))))/denominator
            end do
        end do
    end if

    return

end subroutine AddOutletBC

subroutine CorrectOutletFlux

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays
    use vent_arrays
    use mpih

    implicit none

    integer :: j,k
    real    :: res_dummy

    !Get the inlet flux
    iflux = 0.0d0
    if (xstart(3).eq.1) then
        do k=1,nxm
            do j=xstart(2),xend(2)
                iflux = iflux + (vz(k,j,1)*dx3c(k)/(dx*dy))
            end do
        end do
    end if

    call MpiAllSumRealScalar(iflux,res_dummy)
    iflux = res_dummy

    !Get the outlet flux
    oflux = 0.0d0
    if (xend(3).eq.nzm) then
        do k=1,nxm
            do j=xstart(2),xend(2)
                oflux = oflux + (vz(k,j,nz)*dx3c(k)/(dx*dy))
            end do
        end do
    end if

    call MpiAllSumRealScalar(oflux,res_dummy)
    oflux = res_dummy

    ! Correct the outlet flux
    if (xend(3).eq.nzm) then
        do k=oxfst,oxfen
            do j=xstart(2),xend(2)
                vz(k,j,nz) = vz(k,j,nz) + (ocell(k)*(iflux-oflux)/oarea)
            end do
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
        time,ntime,iflux,iarea,oflux,oarea,((iflux - oflux)/oarea)
        close(100)
    end if
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