!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateH2O.F90                       !
!    CONTAINS: subroutine ImplicitAndUpdateH2O            !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the water vapour and call the implicit solver.      !
!     After this routine, the water vapour has been       !
!     updated to the new timestep                         !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ImplicitAndUpdateH2O

    use param
    use local_arrays, only: h2o,qh2o,ruh2o,rhs
    use decomp_2d, only: xstart,xend

    implicit none

    integer :: jc,kc,ic
    integer :: km,kp
    real    :: alpec,dxxt
    real    :: app,acc,amm

    alpec=al/pec

    !$OMP   PARALLEL DO &
    !$OMP   DEFAULT(none) &
    !$OMP   SHARED(xstart,xend,nxm,h2o) &
    !$OMP   SHARED(kmv,kpv,am3ck,ac3ck,ap3ck) &
    !$OMP   SHARED(ga,ro,alpec,dt) &
    !$OMP   SHARED(rhs,ruh2o,qh2o) &
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

                ! Calculate second derivative of h2oerature in the x-direction.
                ! This is the only term calculated implicitly for h2oerature.

                dxxt = h2o(kp,jc,ic)*ap3ck(kc) + h2o(kc,jc,ic)*ac3ck(kc) + h2o(km,jc,ic)*am3ck(kc)

                ! Calculate right hand side of Eq. 5 (VO96)

                rhs(kc,jc,ic) = (ga*qh2o(kc,jc,ic) + ro*ruh2o(kc,jc,ic) + alpec*dxxt)*dt

                ! Store the non-linear terms for the calculation of 
                ! the next timestep

                ruh2o(kc,jc,ic) = qh2o(kc,jc,ic)

            enddo
        enddo
    enddo

    !$OMP END PARALLEL DO

    !  Solve equation and update water vapour
    call SolveImpEqnUpdate_H2O

    return

    end
