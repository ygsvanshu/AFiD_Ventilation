!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: GlobalQuantities.F90                           !
!    CONTAINS: subroutine GlobalQuantities                !
!                                                         ! 
!    PURPOSE: Calculate maximum velocity and temperature. !
!     volume averaged Nusselt number and Reynolds number  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GlobalQuantities

    use param
    use local_arrays,only: vy,vx,vz,temp,co2,h2o
    use decomp_2d,only: xstart,xend
    use ibm_arrays,only: ibm_body
    use mpih

    implicit none

    integer :: ic,ip,jc,jp,kc,kp
    real    :: volf,volc,ibmc
    real    :: vxcen,vycen,vzcen,tempcen,co2cen,h2ocen
    real    :: vxrms,vyrms,vzrms,vlrms,temprms,co2rms,h2orms
    real    :: vxm,vym,vzm
    real    :: vent_area,cell_area
    real    :: outflux_vz,outflux_temp,outflux_co2,outflux_h2o
    real    :: rradpr
    real    :: res_dummy

    vmax(1) = -huge(0.0d0)
    vmax(2) = -huge(0.0d0)
    vmax(3) = -huge(0.0d0)
    
    tempmax = -huge(0.0d0)
    tempmin =  huge(0.0d0)

    co2max  = -huge(0.0d0)
    co2min  =  huge(0.0d0)

    h2omax  = -huge(0.0d0)
    h2omin  =  huge(0.0d0)

    vxm     = 0.0d0
    vym     = 0.0d0
    vzm     = 0.0d0
    tempm   = 0.0d0
    co2m    = 0.0d0
    h2om    = 0.0d0

    vxrms   = 0.0d0
    vyrms   = 0.0d0
    vzrms   = 0.0d0
    vlrms   = 0.0d0
    temprms = 0.0d0
    co2rms  = 0.0d0
    h2orms  = 0.0d0

    volf    = 0.0d0

    do kc=1,nxm
        kp=kc+1
        do jc=xstart(2),xend(2)
            jp=jc+1
            do ic=xstart(3),xend(3)
                ip=ic+1

                ! Check for not ibm_body
                if (ibm_body(kc,jc,ic).lt.1.0d-8) then

                    ! Compute the cell volume
                    volc = dx3c(kc)/(dx*dy*dz)
                    volf = volf + volc

                    vmax(1) = max(vmax(1),abs(vz(kc,jc,ic)))
                    vmax(2) = max(vmax(2),abs(vy(kc,jc,ic)))
                    vmax(3) = max(vmax(3),abs(vx(kc,jc,ic)))

                    tempmax = max(tempmax,temp(kc,jc,ic))
                    tempmin = min(tempmin,temp(kc,jc,ic))

                    vxcen   = (vx(kc,jc,ic)+vx(kp,jc,ic))*0.5d0
                    vycen   = (vy(kc,jc,ic)+vy(kc,jp,ic))*0.5d0
                    vzcen   = (vz(kc,jc,ic)+vz(kc,jc,ip))*0.5d0
                    tempcen = (temp(kc,jc,ic)+temp(kp,jc,ic))*0.5d0
                    co2cen  = (co2(kc,jc,ic)+co2(kp,jc,ic))*0.5d0
                    h2ocen  = (h2o(kc,jc,ic)+h2o(kp,jc,ic))*0.5d0

                    vxm     = vxm   + volc*vxcen
                    vym     = vym   + volc*vycen
                    vzm     = vzm   + volc*vzcen
                    tempm   = tempm + volc*tempcen
                    co2m    = co2m  + volc*co2cen
                    h2om    = h2om  + volc*h2ocen

                    vxcen   = (vx(kc,jc,ic)**2+vx(kp,jc,ic)**2)*0.5d0
                    vycen   = (vy(kc,jc,ic)**2+vy(kc,jp,ic)**2)*0.5d0
                    vzcen   = (vz(kc,jc,ic)**2+vz(kc,jc,ip)**2)*0.5d0
                    tempcen = (temp(kc,jc,ic)**2+temp(kp,jc,ic)**2)*0.5d0
                    co2cen  = (co2(kc,jc,ic)**2+co2(kp,jc,ic)**2)*0.5d0
                    h2ocen  = (h2o(kc,jc,ic)**2+h2o(kp,jc,ic)**2)*0.5d0

                    vxrms   = vxrms   + volc*vxcen
                    vyrms   = vyrms   + volc*vxcen
                    vzrms   = vzrms   + volc*vxcen
                    vlrms   = vlrms   + volc*(vxcen+vycen+vzcen)
                    temprms = temprms + volc*tempcen
                    co2rms  = co2rms  + volc*co2cen
                    h2orms  = h2orms  + volc*h2ocen

                end if
            enddo
        enddo
    enddo

    ! Compute the flux of Temp,CO2,and H2O from the outflow vent

    outflux_temp  = 0.0d0
    outflux_co2   = 0.0d0
    outflux_h2o   = 0.0d0
    vent_area     = 0.0d0
                
    if (xend(3).eq.nzm) then
        do kc=1,nxm
            do jc=xstart(2),xend(2)
                ! Check for any overlap of grid cell with outlet
                if ((xc(kc+1).gt.(oheight-(0.5d0*olen))).and.(xc(kc).lt.(oheight+(0.5d0*olen)))) then
                    cell_area    = (dx3c(kc)/(dx*dy))
                    vent_area    = vent_area + cell_area
                    tempcen      = (temp(kc,jc,nzm)+temp(kp,jc,nzm)+temp(kc,jc,nz)+temp(kp,jc,nz))*0.25d0
                    co2cen       = (co2(kc,jc,nzm)+co2(kp,jc,nzm)+co2(kc,jc,nz)+co2(kp,jc,nz))*0.25d0
                    h2ocen       = (h2o(kc,jc,nzm)+h2o(kp,jc,nzm)+h2o(kc,jc,nz)+h2o(kp,jc,nz))*0.25d0
                    outflux_vz   = outflux_vz + (vz(kc,jc,nz)*cell_area)
                    outflux_temp = outflux_temp + (tempcen*vz(kc,jc,nz)*cell_area)
                    outflux_co2  = outflux_co2  + (co2cen*vz(kc,jc,nz)*cell_area)
                    outflux_h2o  = outflux_h2o  + (h2ocen*vz(kc,jc,nz)*cell_area)
                end if
            end do
        end do
    end if

    call MpiMaxRealScalar(vmax(1),res_dummy)
    vmax(1) = res_dummy
    call MpiMaxRealScalar(vmax(2),res_dummy)
    vmax(2) = res_dummy
    call MpiMaxRealScalar(vmax(3),res_dummy)
    vmax(3) = res_dummy

    call MpiMinRealScalar(tempmin,res_dummy)
    tempmin = res_dummy
    call MpiMaxRealScalar(tempmax,res_dummy)
    tempmax = res_dummy

    call MpiMinRealScalar(co2min,res_dummy)
    co2min = res_dummy
    call MpiMaxRealScalar(co2max,res_dummy)
    co2max = res_dummy

    call MpiMinRealScalar(h2omin,res_dummy)
    h2omin = res_dummy
    call MpiMaxRealScalar(h2omax,res_dummy)
    h2omax = res_dummy

    call MpiSumRealScalar(vxm,res_dummy)
    vxm = res_dummy
    call MpiSumRealScalar(vym,res_dummy)
    vym = res_dummy
    call MpiSumRealScalar(vzm,res_dummy)
    vzm = res_dummy
    call MpiSumRealScalar(tempm,res_dummy)
    tempm = res_dummy
    call MpiSumRealScalar(co2m,res_dummy)
    co2m = res_dummy
    call MpiSumRealScalar(h2om,res_dummy)
    h2om = res_dummy
    
    call MpiSumRealScalar(vxrms,res_dummy)
    vxrms = res_dummy
    call MpiSumRealScalar(vyrms,res_dummy)
    vyrms = res_dummy
    call MpiSumRealScalar(vzrms,res_dummy)
    vzrms = res_dummy
    call MpiSumRealScalar(vlrms,res_dummy)
    vlrms = res_dummy
    call MpiSumRealScalar(temprms,res_dummy)
    temprms = res_dummy
    call MpiSumRealScalar(co2rms,res_dummy)
    co2rms = res_dummy
    call MpiSumRealScalar(h2orms,res_dummy)
    h2orms = res_dummy

    call MpiSumRealScalar(volf,res_dummy)
    volf = res_dummy

    call MpiSumRealScalar(outflux_vz,res_dummy)
    outflux_vz = res_dummy
    call MpiSumRealScalar(outflux_temp,res_dummy)
    outflux_temp = res_dummy
    call MpiSumRealScalar(outflux_co2,res_dummy)
    outflux_co2 = res_dummy
    call MpiSumRealScalar(outflux_h2o,res_dummy)
    outflux_h2o = res_dummy

    call MpiSumRealScalar(vent_area,res_dummy)
    vent_area = res_dummy
    
    if(ismaster) then

        vxm   = vxm/volf
        vym   = vym/volf
        vzm   = vzm/volf
        tempm = tempm/volf
        co2m  = co2m/volf
        h2om  = h2om/volf

        vxrms   = dsqrt(vxrms/volf)
        vyrms   = dsqrt(vyrms/volf)
        vzrms   = dsqrt(vzrms/volf)
        vlrms   = dsqrt(vlrms/volf)
        temprms = dsqrt(temprms/volf)
        co2rms  = dsqrt(co2rms/volf)
        h2orms  = dsqrt(h2orms/volf)

        outflux_vz   = outflux_vz/vent_area
        outflux_temp = outflux_temp/vent_area
        outflux_co2  = outflux_co2/vent_area
        outflux_h2o  = outflux_h2o/vent_area

        open(95,file='Results/vol_stats.out',status='unknown',position='append',access='sequential')
        if ((ntime.eq.0).and.(.not.readflow)) then
            write(95,'(18(A14,X))') &
            'Time','Mean_Vx','Mean_Vy','Mean_Vz','Mean_Temp','Mean_CO2','Mean_H2O',&
            'RMS_Vx','RMS_Vy','RMS_Vz','RMS_Vel','RMS_Temp','RMS_CO2','RMS_H2O',&
            'Outflux Vz','OutFlux_Temp','OutFlux_CO2','OutFlux_H2O'
        end if
        write(95,'(18(E14.6,X))') &
        time,vxm,vym,vzm,tempm,co2m,h2om,&
        vzrms,vyrms,vzrms,vlrms,temprms,co2rms,h2orms,&
        outflux_vz,outflux_temp,outflux_co2,outflux_h2o
        close(95)

    endif

    return  

    end     
