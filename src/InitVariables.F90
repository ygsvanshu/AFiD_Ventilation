!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: InitVariables.F90                              !
!    CONTAINS: subroutine InitVariables                   !
!                                                         ! 
!    PURPOSE: Initialization routine. Sets to zero all    !
!     variables used in the code                          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InitVariables

    use param
    use local_arrays
    use stat_arrays
    use ventilation_arrays
    use decomp_2d
    use AuxiliaryRoutines

    implicit none
      
    !-------------------------------------------------
    ! Arrays for grid making
    !-------------------------------------------------

    call AllocateReal1DArray(zc,1,nz)
    call AllocateReal1DArray(zm,1,nz)
    call AllocateReal1DArray(ak1,1,nz)
    call AllocateReal1DArray(ao,1,nz)

    call AllocateReal1DArray(yc,1,ny)
    call AllocateReal1DArray(ym,1,ny)
    call AllocateReal1DArray(ak2,1,ny)
    call AllocateReal1DArray(ap,1,ny)

    call AllocateReal1DArray(xc,1,nx)
    call AllocateReal1DArray(xm,1,nx)
    call AllocateReal1DArray(g3rc,1,nx)
    call AllocateReal1DArray(g3rm,1,nx)
    ! =ModR10=Robert=2020-11-02=====================================================
    call AllocateReal1DArray(dx3c,1,nx) ! Additional Grid parameters
    call AllocateReal1DArray(dx3m,1,nx) ! Additional Grid parameters
    ! =End=of=ModR10================================================================

    call AllocateReal1DArray(udx3c,1,nx)
    call AllocateReal1DArray(udx3m,1,nx)

    call AllocateReal1DArray(ap3ck,1,nx)
    call AllocateReal1DArray(ac3ck,1,nx)
    call AllocateReal1DArray(am3ck,1,nx)

    call AllocateReal1DArray(ap3sk,1,nx)
    call AllocateReal1DArray(ac3sk,1,nx)
    call AllocateReal1DArray(am3sk,1,nx)
    
    call AllocateReal1DArray(ap3ssk,1,nx)
    call AllocateReal1DArray(ac3ssk,1,nx)
    call AllocateReal1DArray(am3ssk,1,nx)

    call AllocateReal1DArray(amphk,1,nx)
    call AllocateReal1DArray(acphk,1,nx)
    call AllocateReal1DArray(apphk,1,nx)

    call AllocateInt1dArray(kmc,1,nx)
    call AllocateInt1dArray(kpc,1,nx)
    call AllocateInt1dArray(kmv,1,nx)
    call AllocateInt1dArray(kpv,1,nx)

    ! =ModV13=Vanshu=2020=11=12======================================================
    call AllocateReal1dArray(dcpdxc,1,nx)
    call AllocateReal1dArray(dccdxc,1,nx)
    call AllocateReal1dArray(dcmdxc,1,nx)

    call AllocateReal1dArray(dcpdxm,1,nx)
    call AllocateReal1dArray(dccdxm,1,nx)
    call AllocateReal1dArray(dcmdxm,1,nx)
    ! =End=of=ModV13================================================================

    !-------------------------------------------------
    ! Arrays for statistics    
    !-------------------------------------------------

    !! Modified on 10/02/2020 by Vanshu [Modification #16]
    !! Copied from earlier version on 31/12/2019 by Vanshu [Modification #2]

    if (statcalc) then
        ! X-CUT
        call AllocateReal2DArray(vx_m1_xcut,xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal2DArray(vy_m1_xcut,xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal2DArray(vz_m1_xcut,xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal2DArray(temp_m1_xcut,xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal2DArray(co2_m1_xcut,xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal2DArray(h2o_m1_xcut,xstart(2),xend(2),xstart(3),xend(3))

        call AllocateReal2DArray(vx_m2_xcut,xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal2DArray(vy_m2_xcut,xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal2DArray(vz_m2_xcut,xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal2DArray(temp_m2_xcut,xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal2DArray(co2_m2_xcut,xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal2DArray(h2o_m2_xcut,xstart(2),xend(2),xstart(3),xend(3))

        call AllocateReal2DArray(vx_m3_xcut,xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal2DArray(vy_m3_xcut,xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal2DArray(vz_m3_xcut,xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal2DArray(temp_m3_xcut,xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal2DArray(co2_m3_xcut,xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal2DArray(h2o_m3_xcut,xstart(2),xend(2),xstart(3),xend(3))

        call AllocateReal2DArray(vx_m4_xcut,xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal2DArray(vy_m4_xcut,xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal2DArray(vz_m4_xcut,xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal2DArray(temp_m4_xcut,xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal2DArray(co2_m4_xcut,xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal2DArray(h2o_m4_xcut,xstart(2),xend(2),xstart(3),xend(3))
        
        ! Y-CUT
        call AllocateReal2DArray(vx_m1_ycut,1,nx,xstart(3),xend(3))
        call AllocateReal2DArray(vy_m1_ycut,1,nx,xstart(3),xend(3))
        call AllocateReal2DArray(vz_m1_ycut,1,nx,xstart(3),xend(3))
        call AllocateReal2DArray(temp_m1_ycut,1,nx,xstart(3),xend(3))
        call AllocateReal2DArray(co2_m1_ycut,1,nx,xstart(3),xend(3))
        call AllocateReal2DArray(h2o_m1_ycut,1,nx,xstart(3),xend(3))

        call AllocateReal2DArray(vx_m2_ycut,1,nx,xstart(3),xend(3))
        call AllocateReal2DArray(vy_m2_ycut,1,nx,xstart(3),xend(3))
        call AllocateReal2DArray(vz_m2_ycut,1,nx,xstart(3),xend(3))
        call AllocateReal2DArray(temp_m2_ycut,1,nx,xstart(3),xend(3))
        call AllocateReal2DArray(co2_m2_ycut,1,nx,xstart(3),xend(3))
        call AllocateReal2DArray(h2o_m2_ycut,1,nx,xstart(3),xend(3))

        call AllocateReal2DArray(vx_m3_ycut,1,nx,xstart(3),xend(3))
        call AllocateReal2DArray(vy_m3_ycut,1,nx,xstart(3),xend(3))
        call AllocateReal2DArray(vz_m3_ycut,1,nx,xstart(3),xend(3))
        call AllocateReal2DArray(temp_m3_ycut,1,nx,xstart(3),xend(3))
        call AllocateReal2DArray(co2_m3_ycut,1,nx,xstart(3),xend(3))
        call AllocateReal2DArray(h2o_m3_ycut,1,nx,xstart(3),xend(3))

        call AllocateReal2DArray(vx_m4_ycut,1,nx,xstart(3),xend(3))
        call AllocateReal2DArray(vy_m4_ycut,1,nx,xstart(3),xend(3))
        call AllocateReal2DArray(vz_m4_ycut,1,nx,xstart(3),xend(3))
        call AllocateReal2DArray(temp_m4_ycut,1,nx,xstart(3),xend(3))
        call AllocateReal2DArray(co2_m4_ycut,1,nx,xstart(3),xend(3))
        call AllocateReal2DArray(h2o_m4_ycut,1,nx,xstart(3),xend(3))

        ! Z-CUT
        call AllocateReal2DArray(vx_m1_zcut,1,nx,xstart(2),xend(2))
        call AllocateReal2DArray(vy_m1_zcut,1,nx,xstart(2),xend(2))
        call AllocateReal2DArray(vz_m1_zcut,1,nx,xstart(2),xend(2))
        call AllocateReal2DArray(temp_m1_zcut,1,nx,xstart(2),xend(2))
        call AllocateReal2DArray(co2_m1_zcut,1,nx,xstart(2),xend(2))
        call AllocateReal2DArray(h2o_m1_zcut,1,nx,xstart(2),xend(2))

        call AllocateReal2DArray(vx_m2_zcut,1,nx,xstart(2),xend(2))
        call AllocateReal2DArray(vy_m2_zcut,1,nx,xstart(2),xend(2))
        call AllocateReal2DArray(vz_m2_zcut,1,nx,xstart(2),xend(2))
        call AllocateReal2DArray(temp_m2_zcut,1,nx,xstart(2),xend(2))
        call AllocateReal2DArray(co2_m2_zcut,1,nx,xstart(2),xend(2))
        call AllocateReal2DArray(h2o_m2_zcut,1,nx,xstart(2),xend(2))

        call AllocateReal2DArray(vx_m3_zcut,1,nx,xstart(2),xend(2))
        call AllocateReal2DArray(vy_m3_zcut,1,nx,xstart(2),xend(2))
        call AllocateReal2DArray(vz_m3_zcut,1,nx,xstart(2),xend(2))
        call AllocateReal2DArray(temp_m3_zcut,1,nx,xstart(2),xend(2))
        call AllocateReal2DArray(co2_m3_zcut,1,nx,xstart(2),xend(2))
        call AllocateReal2DArray(h2o_m3_zcut,1,nx,xstart(2),xend(2))

        call AllocateReal2DArray(vx_m4_zcut,1,nx,xstart(2),xend(2))
        call AllocateReal2DArray(vy_m4_zcut,1,nx,xstart(2),xend(2))
        call AllocateReal2DArray(vz_m4_zcut,1,nx,xstart(2),xend(2))
        call AllocateReal2DArray(temp_m4_zcut,1,nx,xstart(2),xend(2))
        call AllocateReal2DArray(co2_m4_zcut,1,nx,xstart(2),xend(2))
        call AllocateReal2DArray(h2o_m4_zcut,1,nx,xstart(2),xend(2))

    end if

    call AllocateReal2DArray(outvx,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo)
    call AllocateReal2DArray(outvy,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo)
    call AllocateReal2DArray(outvz,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo)
    call AllocateReal2DArray(outtemp,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo)
    call AllocateReal2DArray(outco2,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo)
    call AllocateReal2DArray(outh2o,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo)

    call AllocateReal2DArray(dzoutvx,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo)
    call AllocateReal2DArray(dzoutvy,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo)
    call AllocateReal2DArray(dzoutvz,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo)
    call AllocateReal2DArray(dzouttemp,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo)
    call AllocateReal2DArray(dzoutco2,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo)
    call AllocateReal2DArray(dzouth2o,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo)

    !-------------------------------------------------
    ! Arrays with ghost cells
    !-------------------------------------------------
    call AllocateReal3DArray(vy,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
    call AllocateReal3DArray(vz,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
    call AllocateReal3DArray(vx,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
    call AllocateReal3DArray(pr,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
    call AllocateReal3DArray(temp,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
    call AllocateReal3DArray(co2,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
    call AllocateReal3DArray(h2o,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
    call AllocateReal3DArray(dphhalo,1,nxm,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
    
    ! For IBM of person body and breath
    call AllocateReal3DArray(ibm_body,1,nxm,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
    call AllocateReal3DArray(ibm_breath,1,nxm,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)

    !-----------------------------------------------
    ! Arrays without ghost cells
    !-----------------------------------------------
    call AllocateReal3DArray(rhs,1,nx,xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(dph,1,nxm,xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(dq,1,nx,xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(qcap,1,nx,xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(hro,1,nx,xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(qco2,1,nx,xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(qh2o,1,nx,xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(rux,1,nx,xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(ruy,1,nx,xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(ruz,1,nx,xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(rutemp,1,nx,xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(ruco2,1,nx,xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(ruh2o,1,nx,xstart(2),xend(2),xstart(3),xend(3))

    return 

    end   

