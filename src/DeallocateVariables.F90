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
    use vent_arrays
    use movie_indices
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
    call DestroyReal1DArray(dx3c)
    call DestroyReal1DArray(dx3m)

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

    call DestroyReal1dArray(icell)
    call DestroyReal1dArray(ocell)

    call DestroyReal1dArray(igrid)
    call DestroyReal1dArray(ogrid)

    call DestroyReal2DArray(outvx)
    call DestroyReal2DArray(outvy)
    call DestroyReal2DArray(outvz)
    call DestroyReal2DArray(outtemp)
    call DestroyReal2DArray(outco2)
    call DestroyReal2DArray(outh2o)

    call DestroyReal2DArray(outvscx)
    call DestroyReal2DArray(outvscy)
    call DestroyReal2DArray(outvscz)

    call DestroyReal2DArray(mov_xcut)
    call DestroyReal2DArray(mov_ycut)
    call DestroyReal2DArray(mov_zcut)
    call DestroyReal2DArray(mov_icut)
    call DestroyReal2DArray(mov_ocut)

    if (statcalc) then
        call DestroyReal3DArray(stat3d_vx_m1)
        call DestroyReal3DArray(stat3d_vy_m1)
        call DestroyReal3DArray(stat3d_vz_m1)
        call DestroyReal3DArray(stat3d_pr_m1)
        call DestroyReal3DArray(stat3d_temp_m1)
        call DestroyReal3DArray(stat3d_co2_m1)
        call DestroyReal3DArray(stat3d_h2o_m1)
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

    return 

    end   

