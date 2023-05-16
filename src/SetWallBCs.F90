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

    igrid(:) = 0.0d0
    ogrid(:) = 0.0d0

    iarea = ilen*ylen
    oarea = olen*ylen

    varea = 0.0d0
    harea = nym*olen/dy

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

    do k=max(ixcst,2),min(ixcen,nxm)
        igrid(k) = (min(xc(k+1),(iheight+(0.5d0*ilen))) - max(xc(k-1),(iheight-(0.5d0*ilen))))*0.5d0*dx/g3rc(k)
    end do

    do k=max(oxcst,2),min(oxcen,nxm)
        ogrid(k) = (min(xc(k+1),(oheight+(0.5d0*olen))) - max(xc(k-1),(oheight-(0.5d0*olen))))*0.5d0*dx/g3rc(k)
        varea = varea + ((min(xc(k+1),(oheight+(0.5d0*olen))) - max(xc(k-1),(oheight-(0.5d0*olen))))*0.5d0)
    end do

end subroutine InitVents

subroutine SetWallBCs

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays
    use vent_arrays
    use mpih

    implicit none

    integer :: k

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

    if (xstart(3).eq.1) then
        vx(:,:,0) = -vx(:,:,1)                  ! Dirchlet condition
        vy(:,:,0) = -vy(:,:,1)                  ! Dirchlet condition
        temp(:,:,0) = temp(:,:,1)               ! Adiabatic
        co2(:,:,0) = co2(:,:,1)                 ! Zero flux
        h2o(:,:,0) = h2o(:,:,1)                 ! Zero flux
        do k=1,nxm                              
            if ((k.lt.ixfst).or.(k.gt.ixfen)) then
                vz(k,:,1) = 0.0d0               ! Dirchlet condition
            end if
        end do
    end if

    if (xend(3).eq.nzm) then
        do k=1,nx                           
            if ((k.lt.oxcst).or.(k.gt.oxcen)) then
                vx(k,:,nz)   = -vx(k,:,nzm)     ! Dirchlet condition
                temp(k,:,nz) = temp(k,:,nzm)    ! Adiabatic 
                co2(k,:,nz)  = co2(k,:,nzm)     ! Zero flux 
                h2o(k,:,nz)  = h2o(k,:,nzm)     ! Zero flux
            end if
        end do
        do k=1,nxm                              
            if ((k.lt.oxfst).or.(k.gt.oxfen)) then
                vy(k,:,nz)   = -vy(k,:,nzm)     ! Dirchlet condition
                vz(k,:,nz)   = 0.0d0            ! Dirchlet condition
            end if
        end do
    end if

    return
    
end subroutine SetWallBCs

subroutine SetInletBC

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays
    use vent_arrays
    use mpih

    implicit none

    integer :: k,j
    real    :: cvel

    ! This ensures that the inlet velocity slowly increases from zero using a smooth third order spline.
    if (time.lt.(tsvel+tvel)) then
        cvel = isvel + ((ivel-isvel)*((3.0d0*(((time-tsvel)/tvel)**2)) - (2.0d0*(((time-tsvel)/tvel)**3))))
    else 
        cvel = ivel
    end if

    ! Wall with inlet vent
    if (xstart(3).eq.1) then
        do k=ixfst,ixfen
            vz(k,:,1)  = cvel*icell(k)  ! Dirchlet condition
        end do
    end if

end subroutine SetInletBC

subroutine CalcOutletBC

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vx,vy,vz,temp,co2,h2o
    use vent_arrays
    use mpih

    implicit none

    integer :: k
    real    :: n1,n2,n3,dd

    n1 = 1.0d0/(al*dt)
    n2 = 0.5d0*ocou*dz
    n3 = 0.0d0 !0.5d0*dzq/ren
    dd = 1.0d0/(n1 + n2 - n3)
    
    if (xend(3).eq.nzm) then
        do k=oxcst,oxcen
            outvx(k,:)   = (n1*vx(k,:,nz))   - (n2*(vx(k,:,nz)  - vx(k,:,nzm))  ) + (n3*(vx(k,:,nz)   - 2.0d0*vx(k,:,nzm)   + vx(k,:,nzm-1))  )
            outtemp(k,:) = (n1*temp(k,:,nz)) - (n2*(temp(k,:,nz)- temp(k,:,nzm))) + (n3*(temp(k,:,nz) - 2.0d0*temp(k,:,nzm) + temp(k,:,nzm-1)))
            outco2(k,:)  = (n1*co2(k,:,nz))  - (n2*(co2(k,:,nz) - co2(k,:,nzm)) ) + (n3*(co2(k,:,nz)  - 2.0d0*co2(k,:,nzm)  + co2(k,:,nzm-1)) )
            outh2o(k,:)  = (n1*h2o(k,:,nz))  - (n2*(h2o(k,:,nz) - h2o(k,:,nzm)) ) + (n3*(h2o(k,:,nz)  - 2.0d0*h2o(k,:,nzm)  + h2o(k,:,nzm-1)) )
        end do

        do k=oxfst,oxfen
            outvy(k,:)   = (n1*vy(k,:,nz))   - (n2*(vy(k,:,nz)  - vy(k,:,nzm))  ) + (n3*(vy(k,:,nz)   - 2.0d0*vy(k,:,nzm)   + vy(k,:,nzm-1))  )
            outvz(k,:)   = (n1*vz(k,:,nz))   - (n2*(vz(k,:,nz)  - vz(k,:,nzm))  ) + (n3*(vz(k,:,nz)   - 2.0d0*vz(k,:,nzm)   + vz(k,:,nzm-1))  )
        end do
    end if

end subroutine CalcOutletBC

subroutine SetOutletBC

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vx,vy,vz,temp,co2,h2o
    use vent_arrays
    use mpih

    implicit none

    integer :: k
    real    :: n1,n2,n3,dd
    
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

    n1 = 1.0d0/(al*dt)
    n2 = 0.5d0*ocou*dz
    n3 = 0.0d0 !0.5d0*dzq/ren
    dd = 1.0d0/(n1 + n2 - n3)

    ! Wall with outlet vent
    if (xend(3).eq.nzm) then
        do k=oxcst,oxcen
            vx(k,:,nz)   =  dd*(outvx(k,:)   + (n2*vx(k,:,nzm)  ) - (n3*(2.0d0*vx(k,:,nzm)   - vx(k,:,nzm-1))  ))
            temp(k,:,nz) =  dd*(outtemp(k,:) + (n2*temp(k,:,nzm)) - (n3*(2.0d0*temp(k,:,nzm) - temp(k,:,nzm-1))))
            co2(k,:,nz)  =  dd*(outco2(k,:)  + (n2*co2(k,:,nzm) ) - (n3*(2.0d0*co2(k,:,nzm)  - co2(k,:,nzm-1)) ))
            h2o(k,:,nz)  =  dd*(outh2o(k,:)  + (n2*h2o(k,:,nzm) ) - (n3*(2.0d0*h2o(k,:,nzm)  - h2o(k,:,nzm-1)) ))
        end do

        do k=oxfst,oxfen
            vy(k,:,nz)   =  dd*(outvy(k,:)   + (n2*vy(k,:,nzm)  ) - (n3*(2.0d0*vy(k,:,nzm)   - vy(k,:,nzm-1))  ))
            vz(k,:,nz)   =  dd*(outvz(k,:)   + (n2*vz(k,:,nzm)  ) - (n3*(2.0d0*vz(k,:,nzm)   - vz(k,:,nzm-1))  ))
        end do
    end if
            
end subroutine SetOutletBC

subroutine CorrectOutletFlux

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays
    use vent_arrays
    use mpih

    implicit none

    integer :: j,k
    real    :: vflux,hflux
    real    :: res_dummy

    ! =============================================================================================
    ! vflux = 0.0d0

    ! if (xend(3).eq.nzm) then
    !     do k=1,nxm
    !         do j=xstart(2),xend(2)
    !             vflux = vflux + (0.25d0*(vx(k,j,nzm)+vx(k+1,j,nzm)+vx(k,j,nz)+vx(k+1,j,nz))*dx3c(k)/(dx*dy))
    !         end do
    !     end do
    ! end if

    ! call MpiAllSumRealScalar(vflux,res_dummy)
    ! vflux = res_dummy

    ! ! Correct the vertical flux
    ! if (xend(3).eq.nzm) then
    !     do k=max(2,oxcst),min(nxm,oxcen)
    !         do j=xstart(2),xend(2)
    !             vx(k,j,nz) = vx(k,j,nz) - (2.0d0*ogrid(k)*vflux/varea)
    !         end do
    !     end do
    ! end if

    ! vflux = 0.0d0

    ! if (xend(3).eq.nzm) then
    !     do k=1,nxm
    !         do j=xstart(2),xend(2)
    !             vflux = vflux + (0.25d0*(vx(k,j,nzm)+vx(k+1,j,nzm)+vx(k,j,nz)+vx(k+1,j,nz))*dx3c(k)/(dx*dy))
    !         end do
    !     end do
    ! end if

    ! call MpiAllSumRealScalar(vflux,res_dummy)
    ! vflux = res_dummy

    ! =============================================================================================

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
            write(100,'(8(A14,X))') &
            'Time','Step','InFlux','InArea','OutFlux','OutArea','Correction','VFLUX'
        end if
        write(100,'((E14.6,X),(I14,X),6(E14.6,X))') &
        time,ntime,iflux,iarea,oflux,oarea,((iflux - oflux)/oarea),vflux
        close(100)
    end if
    ! DEBUG PURPOSES ==============================================================================

end subroutine CorrectOutletFlux

subroutine OutVisc(x,z,visc)

    use param

    implicit none

    real, intent(in) :: x,z
    real, intent(out) :: visc
    real    :: delx,delz,dist
    real    :: obot,otop

    visc = 1.0d0

    obot = (oheight - (0.5d0*olen))
    otop = (oheight + (0.5d0*olen))

    if (x.lt.obot) then
        delx = obot - x
    else if (x.gt.otop) then
        delx = x - otop
    else 
        delx = 0.0d0
    end if

    delz = z - zlen

    dist = ((delx**2.0d0) + (delz**2.0d0))**0.5d0

    if (dist.lt.0) then
        visc = ovsc
    else if (dist.gt.odst) then
        visc = 1
    else
        visc = ((ovsc-1.0d0)*((2.0d0*((dist/odst)**3)) - (3.0d0*((dist/odst)**2)) + 1)) + 1
    end if

    return

end subroutine OutVisc

subroutine InitVisc

    use param
    use decomp_2d, only: xstart,xend
    use vent_arrays, only: outvscx,outvscy,outvscz

    implicit none

    integer :: i,j,k
    real    :: visc

    outvscx(:,:) = 1.0d0
    outvscy(:,:) = 1.0d0
    outvscz(:,:) = 1.0d0

    do k=1,nxm
        do i=xstart(3),xend(3)
            
            call OutVisc(xc(k),zm(i),visc)
            outvscx(k,i) = visc

            call OutVisc(xm(k),zm(i),visc)
            outvscy(k,i) = visc

            call OutVisc(xm(k),zc(i),visc)
            outvscz(k,i) = visc

        end do
    end do

    return

end subroutine InitVisc

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
