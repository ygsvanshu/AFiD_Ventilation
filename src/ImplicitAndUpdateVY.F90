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
    use local_arrays, only: vy,ruy,pr,rhsx,dph
    use ibm_arrays, only: ibm_gy_px
    use vent_arrays, only: oxfst,oxfen
    use decomp_2d, only: xstart,xend

    implicit none

    integer :: im,ic,ip
    integer :: jm,jc,jp
    integer :: km,kc,kp
    real    :: alre,udy
    real    :: amm,acc,app
    real    :: dyp,dxxvy,dyyvy,dzzvy

    alre=al/ren
    udy=dy*al

    !$OMP   PARALLEL DO &
    !$OMP   DEFAULT(none) &
    !$OMP   SHARED(xstart,xend,nxm,vy,pr) &
    !$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk) &
    !$OMP   SHARED(dy,al,ga,ro,alre,dt,dph) &
    !$OMP   SHARED(udy,dyq,dzq,rhsx,ruy) &
    !$OMP   PRIVATE(ic,jc,kc,im,jm,km,ip,jp,kp) &
    !$OMP   PRIVATE(amm,acc,app) &
    !$OMP   PRIVATE(dyp,dxxvy,dyyvy,dzzvy)

    do ic=xstart(3),xend(3)
        im=ic-1
        ip=ic+1
        do jc=xstart(2),xend(2)
            jm=jc-1
            jp=jc+1
            do kc=1,nxm
                km=kmv(kc)
                kp=kpv(kc)

                amm=am3sk(kc)
                acc=ac3sk(kc)
                app=ap3sk(kc)

                !   Second derivative in x-direction of vy
                if (kc.eq.1) then
                    dxxvy = vy(kp,jc,ic)*app + vy(kc,jc,ic)*acc + 0.d0*amm
                else if (kc.eq.nxm) then
                    dxxvy = 0.0*app + vy(kc,jc,ic)*acc + vy(km,jc,ic)*amm
                else
                    dxxvy = vy(kp,jc,ic)*app + vy(kc,jc,ic)*acc + vy(km,jc,ic)*amm
                end if

                !   Second derivative in y-direction of vy
                dyyvy=(vy(kc,jp,ic) - 2.0*vy(kc,jc,ic) + vy(kc,jm,ic))*dyq

                !   Second derivative in z-direction of vy
                dzzvy=(vy(kc,jc,ip) - 2.0*vy(kc,jc,ic) + vy(kc,jc,im))*dzq

                !   component of grad(pr) along y direction
                dyp = (pr(kc,jc,ic)-pr(kc,jm,ic))*udy

                !    Calculate right hand side of Eq. 5 (VO96)
                rhsx(kc,jc,ic) = (ga*dph(kc,jc,ic) + ro*ruy(kc,jc,ic) + alre*(dxxvy+dyyvy+dzzvy) - dyp)*dt

                !    Store the non-linear terms for the calculation of the next timestep
                ruy(kc,jc,ic)=dph(kc,jc,ic)
            enddo
        enddo
    enddo

    !$OMP END PARALLEL DO

    call AddBodyIBM(vy,ibm_gy_px,0.0d0)
    call AddBC_Vy
    call AddOutletBC(vy,ren,oxfst,oxfen)
    call SolveImpEqnUpdate_Y

    return
end
