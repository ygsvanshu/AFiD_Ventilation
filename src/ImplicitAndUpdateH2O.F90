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
    use local_arrays, only: h2o,qh2o,ruh2o,rhsx
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
    !$OMP   SHARED(xstart,xend,nxm,h2o) &
    !$OMP   SHARED(dyq,dzq,am3ck,ac3ck,ap3ck) &
    !$OMP   SHARED(ga,ro,alpec,dt) &
    !$OMP   SHARED(rhs,ruh2o,qh2o) &
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

                ! Second derivative in x-direction of h2o
                dxxt = h2o(kp,jc,ic)*ap3ssk(kc) + h2o(kc,jc,ic)*ac3ssk(kc) + h2o(km,jc,ic)*am3ssk(kc)
                ! Second derivative in y-direction of h2o
                dzzt = (h2o(kc,jc,ip) -2.0*h2o(kc,jc,ic) + h2o(kc,jc,im))*dzq
                ! Second derivative in z-direction of h2o
                dyyt = (h2o(kc,jp,ic) -2.0*h2o(kc,jc,ic) + h2o(kc,jm,ic))*dyq

                !    Calculate right hand side of Eq. 5 (VO96)
                rhsx(kc,jc,ic) = (ga*qh2o(kc,jc,ic) + ro*ruh2o(kc,jc,ic) + alpec*(dxxt+dyyt+dzzt))*dt

                !    Store the non-linear terms for the calculation of the next timestep
                ruh2o(kc,jc,ic) = qh2o(kc,jc,ic)

            enddo
        enddo
    enddo

    !$OMP END PARALLEL DO

    call AddBodyIBM(h2o,ibm_gx_px,0.0d0)
    call AddOutletBC(h2o,pec,oxcst,oxcen)
    call SolveImpEqnUpdate_H2O

    return

    end
