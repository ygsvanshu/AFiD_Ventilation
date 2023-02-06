!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateCO2.F90                       !
!    CONTAINS: subroutine ImplicitAndUpdateCO2            !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the carbon dioxide and call the implicit solver.    !
!     After this routine, the carbon dioxide has been     !
!     updated to the new timestep                         !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ImplicitAndUpdateCO2

    use param
    use local_arrays, only: co2,qco2,ruco2,rhsx
    use ibm_arrays, only: ibm_gx_px
    use vent_arrays, only: oxcst,oxcen
    use decomp_2d, only: xstart,xend

    implicit none

    integer :: im,ic,ip
    integer :: jm,jc,jp
    integer :: km,kc,kp
    real    :: alpec
    real    :: dxxt,dyyt,dzzt
    real    :: app,acc,amm

    alpec=al/pec

    !$OMP   PARALLEL DO &
    !$OMP   DEFAULT(none) &
    !$OMP   SHARED(xstart,xend,nxm,co2) &
    !$OMP   SHARED(dyq,dzq,am3ck,ac3ck,ap3ck) &
    !$OMP   SHARED(ga,ro,alpec,dt) &
    !$OMP   SHARED(rhs,ruco2,qco2) &
    !$OMP   PRIVATE(ic,jc,kc,im,jm,km,ip,jp,kp) &
    !$OMP   PRIVATE(amm,acc,app) &
    !$OMP   PRIVATE(dxxt,dyyt,dzzt)

    do ic=xstart(3),xend(3)
        im=ic-1
        ip=ic+1
        do jc=xstart(2),xend(2)
            jm=jc-1
            jp=jc+1
            do kc=2,nxm
                km=kc-1
                kp=kc+1

                if (kc.eq.1)  km=kc+1
                if (kc.eq.nx) kp=kc-1

                ! Second derivative in x-direction of co2
                dxxt = co2(kp,jc,ic)*ap3ssk(kc) + co2(kc,jc,ic)*ac3ssk(kc) + co2(km,jc,ic)*am3ssk(kc)
                ! Second derivative in y-direction of co2
                dzzt = (co2(kc,jc,ip) -2.0*co2(kc,jc,ic) + co2(kc,jc,im))*dzq
                ! Second derivative in z-direction of co2
                dyyt = (co2(kc,jp,ic) -2.0*co2(kc,jc,ic) + co2(kc,jm,ic))*dyq

                !    Calculate right hand side of Eq. 5 (VO96)
                rhsx(kc,jc,ic) = (ga*qco2(kc,jc,ic) + ro*ruco2(kc,jc,ic) + alpec*(dxxt+dyyt+dzzt))*dt

                !    Store the non-linear terms for the calculation of the next timestep
                ruco2(kc,jc,ic) = qco2(kc,jc,ic)

            enddo
        enddo
    enddo

    !$OMP END PARALLEL DO

    call AddBodyIBM(co2,ibm_gx_px,0.0d0)
    call AddOutletBC(co2,pec,oxcst,oxcen)
    call SolveImpEqnUpdate_CO2

    return

    end
