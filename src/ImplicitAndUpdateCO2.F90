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
    use local_arrays, only: co2,qco2,ruco2,rhs
    use decomp_2d, only: xstart,xend

    implicit none

    integer :: jc,kc,ic
    integer :: km,kp
    real    :: alpec,dxxt
    real    :: app,acc,amm

    alpec=al/pec

    !$OMP   PARALLEL DO &
    !$OMP   DEFAULT(none) &
    !$OMP   SHARED(xstart,xend,nxm,co2) &
    !$OMP   SHARED(kmv,kpv,am3ck,ac3ck,ap3ck) &
    !$OMP   SHARED(ga,ro,alpec,dt) &
    !$OMP   SHARED(rhs,ruco2,qco2) &
    !$OMP   PRIVATE(ic,jc,kc,km,kp) &
    !$OMP   PRIVATE(amm,acc,app) &
    !$OMP   PRIVATE(dxxt)

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nx

                km = kc - 1
                kp = kc + 1

                if (kc.eq.1)  km = kc + 1
                if (kc.eq.nx) kp = kc - 1

                ! Calculate second derivative of co2erature in the x-direction.
                ! This is the only term calculated implicitly for co2erature.

                dxxt = co2(kp,jc,ic)*ap3ck(kc) + co2(kc,jc,ic)*ac3ck(kc) + co2(km,jc,ic)*am3ck(kc)

                ! Calculate right hand side of Eq. 5 (VO96)

                rhs(kc,jc,ic) = (ga*qco2(kc,jc,ic) + ro*ruco2(kc,jc,ic) + alpec*dxxt)*dt

                ! Store the non-linear terms for the calculation of 
                ! the next timestep

                ruco2(kc,jc,ic) = qco2(kc,jc,ic)

            enddo
        enddo
    enddo

    !$OMP END PARALLEL DO

    !  Solve equation and update co2
    call SolveImpEqnUpdate_CO2

    return

    end
