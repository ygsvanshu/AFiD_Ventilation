!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateVY.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateVY             !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the y (horizontal) dimension        !
!     and call the implicit solver                        !
!     After this routine, the velocity field in y has     !
!     been updated to the new timestep                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ImplicitAndUpdateVY

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vy,ruy,pr,rhs,dph
    use vent_arrays, only: outvscy

    implicit none

    integer :: kc,jmm,jc,ic
    integer :: kpp,kmm
    real    :: alre,udy
    real    :: amm,acc,app
    real    :: dyp,dxxvy

    alre=al/ren
    udy=dy*al

    !$OMP   PARALLEL DO &
    !$OMP   DEFAULT(none) &
    !$OMP   SHARED(xstart,xend,nxm,vy,pr) &
    !$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk) &
    !$OMP   SHARED(dy,al,ga,ro,alre,dt,dph) &
    !$OMP   SHARED(udy,udx3m,rhs,ruy) &
    !$OMP   PRIVATE(ic,jc,kc,kmm,kpp,jmm) &
    !$OMP   PRIVATE(amm,acc,app) &
    !$OMP   PRIVATE(dyp,dxxvy)

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            jmm=jc-1
            do kc=1,nxm

                kmm=kmv(kc)
                kpp=kpv(kc)
                amm=am3sk(kc)
                acc=ac3sk(kc)
                app=ap3sk(kc)

                !   Second derivative in x-direction of vy
                
                if (kc.eq.1) then
                    dxxvy = vy(kpp,jc,ic)*app + vy(kc,jc,ic)*acc + 0.d0*amm
                else if (kc.eq.nxm) then
                    dxxvy = 0.0*app + vy(kc,jc,ic)*acc + vy(kmm,jc,ic)*amm
                else
                    dxxvy = vy(kpp,jc,ic)*app + vy(kc,jc,ic)*acc + vy(kmm,jc,ic)*amm
                end if

                !   component of grad(pr) along y direction

                dyp = (pr(kc,jc,ic)-pr(kc,jmm,ic))*udy

                !    Calculate right hand side of Eq. 5 (VO96)

                rhs(kc,jc,ic) = (ga*dph(kc,jc,ic) + ro*ruy(kc,jc,ic) + outvscy(kc,ic)*alre*dxxvy - dyp)*dt

                !    Store the non-linear terms for the calculation of 
                !    the next timestep

                ruy(kc,jc,ic)=dph(kc,jc,ic)
            enddo
        enddo
    enddo

    !$OMP END PARALLEL DO

    call SolveImpEqnUpdate_Y

    return
end
