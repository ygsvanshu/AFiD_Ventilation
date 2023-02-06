!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateTemp.F90                      !
!    CONTAINS: subroutine ImplicitAndUpdateTemp           !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the temperature and call the implicit solver.       !
!     After this routine, the temperature has been        !
!     updated to the new timestep                         !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ImplicitAndUpdateTemp

    use param
    use local_arrays, only: temp,hro,rutemp,rhsx
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
    !$OMP   SHARED(xstart,xend,nxm,temp) &
    !$OMP   SHARED(dyq,dzq,am3ck,ac3ck,ap3ck) &
    !$OMP   SHARED(ga,ro,alpec,dt) &
    !$OMP   SHARED(rhs,rutemp,hro) &
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

                ! Second derivative in x-direction of temp
                dxxt = temp(kp,jc,ic)*ap3ssk(kc) + temp(kc,jc,ic)*ac3ssk(kc) + temp(km,jc,ic)*am3ssk(kc)
                ! Second derivative in y-direction of temp
                dzzt = (temp(kc,jc,ip) -2.0*temp(kc,jc,ic) + temp(kc,jc,im))*dzq
                ! Second derivative in z-direction of temp
                dyyt = (temp(kc,jp,ic) -2.0*temp(kc,jc,ic) + temp(kc,jm,ic))*dyq

                !    Calculate right hand side of Eq. 5 (VO96)
                rhsx(kc,jc,ic) = (ga*hro(kc,jc,ic) + ro*rutemp(kc,jc,ic) + alpec*(dxxt+dyyt+dzzt))*dt

                !    Store the non-linear terms for the calculation of the next timestep
                rutemp(kc,jc,ic) = hro(kc,jc,ic)

            enddo
        enddo
    enddo

    !$OMP END PARALLEL DO

    call AddBodyIBM(temp,ibm_gx_px,1.0d0)
    call AddOutletBC(temp,pec,oxcst,oxcen)
    call SolveImpEqnUpdate_Temp

    return

    end
