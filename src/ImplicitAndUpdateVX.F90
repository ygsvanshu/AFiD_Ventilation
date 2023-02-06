!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateVX.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateVX             !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the X (vertical) direction and call !
!     the implicit solver. After this routine, the        !
!     vertical velocity has been updated to the new       !
!     timestep                                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ImplicitAndUpdateVX

    use param
    use local_arrays, only: vx,rhsx,rux,qcap,pr
    use ibm_arrays, only: ibm_gx_px
    use vent_arrays, only: oxcst,oxcen
    use decomp_2d, only: xstart,xend

    implicit none

    integer :: im,ic,ip
    integer :: jm,jc,jp
    integer :: km,kc,kp
    real    :: alre,udx3
    real    :: amm,acc,app,dxp
    real    :: dxxvx,dyyvx,dzzvx

    alre=al/ren

    !$OMP   PARALLEL DO &
    !$OMP   DEFAULT(none) &
    !$OMP   SHARED(xstart,xend,nxm,vx,pr) &
    !$OMP   SHARED(am3ck,ac3ck,ap3ck) &
    !$OMP   SHARED(al,ga,ro,alre,dt,qcap) &
    !$OMP   SHARED(udx3c,dyq,dzq,rhsx,rux) &
    !$OMP   PRIVATE(ic,jc,kc,im,jm,km,ip,jp,kp) &
    !$OMP   PRIVATE(amm,acc,app,udx3) &
    !$OMP   PRIVATE(dxxvx,dyyvx,dzzvx,dxp)
    
    do ic=xstart(3),xend(3)
        im=ic-1
        ip=ic+1
        do jc=xstart(2),xend(2)
            jm=jc-1
            jp=jc+1
            do kc=2,nxm
                km=kc-1
                kp=kc+1

                udx3 = al*udx3c(kc)

                amm=am3ck(kc)
                acc=ac3ck(kc)
                app=ap3ck(kc)

                ! Second derivative in x-direction of vx
                dxxvx = vx(kp,jc,ic)*app + vx(kc,jc,ic)*acc + vx(km,jc,ic)*amm
                ! Second derivative in y-direction of vx
                dyyvx = (vx(kc,jm,ic) -2.0*vx(kc,jc,ic) + vx(kc,jp,ic))*dyq
                ! Second derivative in z-direction of vx
                dzzvx = (vx(kc,jc,im) -2.0*vx(kc,jc,ic) + vx(kc,jc,ip))*dzq
                
                ! Component of grad(pr) along x direction
                dxp = (pr(kc,jc,ic) - pr(km,jc,ic))*udx3

                ! Calculate right hand side of Eq. 5 (VO96)
                rhsx(kc,jc,ic) = (ga*qcap(kc,jc,ic) + ro*rux(kc,jc,ic) + alre*(dxxvx+dyyvx+dzzvx) - dxp)*dt 

                ! Store the non-linear terms for the calculation of the next timestep
                rux(kc,jc,ic)=qcap(kc,jc,ic)
            enddo
        enddo
    enddo

    !$OMP END PARALLEL DO

    ! Solve equation and update velocity
    call AddBodyIBM(vx,ibm_gx_px,0.0d0)
    call AddBC_Vx
    call AddOutletBc(vx,ren,oxcst,oxcen)
    call SolveImpEqnUpdate_X

    return

    end
