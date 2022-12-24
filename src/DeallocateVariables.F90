!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: DeallocateVariables.F90                        !
!    CONTAINS: subroutine DeallocateVariables             !
!                                                         ! 
!    PURPOSE: Finalization routine. Deallocates all       !
!     variables used in the code                          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine DeallocateVariables

    use param
    use local_arrays
    use stat_arrays
    use ibm_arrays
    use AuxiliaryRoutines

    implicit none
    
    call DestroyReal1DArray(zc)
    call DestroyReal1DArray(zm)
    call DestroyReal1DArray(ak1)
    call DestroyReal1DArray(ao)

    call DestroyReal1DArray(yc)
    call DestroyReal1DArray(ym)
    call DestroyReal1DArray(ak2)
    call DestroyReal1DArray(ap)

    call DestroyReal1DArray(xc)
    call DestroyReal1DArray(xm)
    call DestroyReal1DArray(g3rc)
    call DestroyReal1DArray(g3rm)
    ! =ModR10=Robert=2020-11-02=====================================================
    call DestroyReal1DArray(dx3c) ! Additional Grid parameters
    call DestroyReal1DArray(dx3m) ! Additional Grid parameters
    ! =End=of=ModR10================================================================

    call DestroyReal1DArray(udx3c)
    call DestroyReal1DArray(udx3m)

    call DestroyReal1DArray(ap3ck)
    call DestroyReal1DArray(ac3ck)
    call DestroyReal1DArray(am3ck)

    call DestroyReal1DArray(ap3sk)
    call DestroyReal1DArray(ac3sk)
    call DestroyReal1DArray(am3sk)

    call DestroyReal1DArray(ap3ssk)
    call DestroyReal1DArray(ac3ssk)
    call DestroyReal1DArray(am3ssk)

    call DestroyReal1DArray(amphk)
    call DestroyReal1DArray(acphk)
    call DestroyReal1DArray(apphk)

    call DestroyInt1dArray(kmc)
    call DestroyInt1dArray(kpc)
    call DestroyInt1dArray(kmv)
    call DestroyInt1dArray(kpv)

    ! =ModV13=Vanshu=2020=11=12=====================================================
    call DestroyReal1dArray(dcpdxc)
    call DestroyReal1dArray(dccdxc)
    call DestroyReal1dArray(dcmdxc)

    call DestroyReal1dArray(dcpdxm)
    call DestroyReal1dArray(dccdxm)
    call DestroyReal1dArray(dcmdxm)

    if (statcalc) then
        ! X-CUT
        call DestroyReal2DArray(vx_m1_xcut)
        call DestroyReal2DArray(vy_m1_xcut)
        call DestroyReal2DArray(vz_m1_xcut)
        call DestroyReal2DArray(temp_m1_xcut)
        call DestroyReal2DArray(co2_m1_xcut)
        call DestroyReal2DArray(h2o_m1_xcut)

        call DestroyReal2DArray(vx_m2_xcut)
        call DestroyReal2DArray(vy_m2_xcut)
        call DestroyReal2DArray(vz_m2_xcut)
        call DestroyReal2DArray(temp_m2_xcut)
        call DestroyReal2DArray(co2_m2_xcut)
        call DestroyReal2DArray(h2o_m2_xcut)

        call DestroyReal2DArray(vx_m3_xcut)
        call DestroyReal2DArray(vy_m3_xcut)
        call DestroyReal2DArray(vz_m3_xcut)
        call DestroyReal2DArray(temp_m3_xcut)
        call DestroyReal2DArray(co2_m3_xcut)
        call DestroyReal2DArray(h2o_m3_xcut)

        call DestroyReal2DArray(vx_m4_xcut)
        call DestroyReal2DArray(vy_m4_xcut)
        call DestroyReal2DArray(vz_m4_xcut)
        call DestroyReal2DArray(temp_m4_xcut)
        call DestroyReal2DArray(co2_m4_xcut)
        call DestroyReal2DArray(h2o_m4_xcut)
        
        ! Y-CUT
        call DestroyReal2DArray(vx_m1_ycut)
        call DestroyReal2DArray(vy_m1_ycut)
        call DestroyReal2DArray(vz_m1_ycut)
        call DestroyReal2DArray(temp_m1_ycut)
        call DestroyReal2DArray(co2_m1_ycut)
        call DestroyReal2DArray(h2o_m1_ycut)

        call DestroyReal2DArray(vx_m2_ycut)
        call DestroyReal2DArray(vy_m2_ycut)
        call DestroyReal2DArray(vz_m2_ycut)
        call DestroyReal2DArray(temp_m2_ycut)
        call DestroyReal2DArray(co2_m2_ycut)
        call DestroyReal2DArray(h2o_m2_ycut)

        call DestroyReal2DArray(vx_m3_ycut)
        call DestroyReal2DArray(vy_m3_ycut)
        call DestroyReal2DArray(vz_m3_ycut)
        call DestroyReal2DArray(temp_m3_ycut)
        call DestroyReal2DArray(co2_m3_ycut)
        call DestroyReal2DArray(h2o_m3_ycut)

        call DestroyReal2DArray(vx_m4_ycut)
        call DestroyReal2DArray(vy_m4_ycut)
        call DestroyReal2DArray(vz_m4_ycut)
        call DestroyReal2DArray(temp_m4_ycut)
        call DestroyReal2DArray(co2_m4_ycut)
        call DestroyReal2DArray(h2o_m4_ycut)

        ! Z-CUT
        call DestroyReal2DArray(vx_m1_zcut)
        call DestroyReal2DArray(vy_m1_zcut)
        call DestroyReal2DArray(vz_m1_zcut)
        call DestroyReal2DArray(temp_m1_zcut)
        call DestroyReal2DArray(co2_m1_zcut)
        call DestroyReal2DArray(h2o_m1_zcut)

        call DestroyReal2DArray(vx_m2_zcut)
        call DestroyReal2DArray(vy_m2_zcut)
        call DestroyReal2DArray(vz_m2_zcut)
        call DestroyReal2DArray(temp_m2_zcut)
        call DestroyReal2DArray(co2_m2_zcut)
        call DestroyReal2DArray(h2o_m2_zcut)

        call DestroyReal2DArray(vx_m3_zcut)
        call DestroyReal2DArray(vy_m3_zcut)
        call DestroyReal2DArray(vz_m3_zcut)
        call DestroyReal2DArray(temp_m3_zcut)
        call DestroyReal2DArray(co2_m3_zcut)
        call DestroyReal2DArray(h2o_m3_zcut)

        call DestroyReal2DArray(vx_m4_zcut)
        call DestroyReal2DArray(vy_m4_zcut)
        call DestroyReal2DArray(vz_m4_zcut)
        call DestroyReal2DArray(temp_m4_zcut)
        call DestroyReal2DArray(co2_m4_zcut)
        call DestroyReal2DArray(h2o_m4_zcut)
    end if

    call DestroyReal3DArray(vx)
    call DestroyReal3DArray(vy)
    call DestroyReal3DArray(vz)
    call DestroyReal3DArray(temp)
    call DestroyReal3DArray(co2)
    call DestroyReal3DArray(h2o)

    call DestroyReal3DArray(pr)
    call DestroyReal3DArray(rhs)

    call DestroyReal3DArray(dq)
    call DestroyReal3DArray(qcap)
    call DestroyReal3DArray(dph)
    call DestroyReal3DArray(hro)
    call DestroyReal3DArray(qco2)
    call DestroyReal3DArray(qh2o)
    call DestroyReal3DArray(dphhalo)

    call DestroyReal3DArray(rux)
    call DestroyReal3DArray(ruy)
    call DestroyReal3DArray(ruz)
    call DestroyReal3DArray(rutemp)
    call DestroyReal3DArray(ruco2)
    call DestroyReal3DArray(ruh2o)
    
    call DestroyReal3DArray(ibm_body)
    call DestroyReal3DArray(ibm_breath)

    return 

    end   

