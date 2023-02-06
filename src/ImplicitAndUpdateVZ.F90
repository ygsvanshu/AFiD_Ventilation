!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateVZ.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateVZ             !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the z (horizontal) dimension        !
!     and call the implicit solver.                       !
!     After this routine, the velocity field in z has     !
!     been updated to the new timestep                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ImplicitAndUpdateVZ

    use param
    use local_arrays, only: vz,dq,ruz,rhsx,pr
    use ibm_arrays, only: ibm_gz_px
    use vent_arrays, only: oxfst,oxfen
    use decomp_2d, only: xstart,xend

    implicit none

    integer :: im,ic,ip
    integer :: jm,jc,jp
    integer :: km,kc,kp
    real    :: alre,amm,acc,app,udz
    real    :: dxxvz,dyyvz,dzzvz,dzp

    alre=al/ren
    udz=dz*al

    !$OMP   PARALLEL DO &
    !$OMP   DEFAULT(none) &
    !$OMP   SHARED(xstart,xend,nxm,vz,pr) &
    !$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk) &
    !$OMP   SHARED(dz,al,ga,ro,alre,dt,dq) &
    !$OMP   SHARED(udz,dyq,dzq,rhsx,ruz) &
    !$OMP   PRIVATE(ic,jc,kc,im,jm,km,ip,jp,kp) &
    !$OMP   PRIVATE(amm,acc,app) &
    !$OMP   PRIVATE(dxxvz,dzzvz,dyyvz,dzp)

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

                !   Second derivative in x-direction of vz
                if (kc.eq.1) then
                    dxxvz = vz(kp,jc,ic)*app + vz(kc,jc,ic)*acc + 0.0d0*amm
                else if (kc.eq.nxm) then           
                    dxxvz = 0.0d0*app + vz(kc,jc,ic)*acc + vz(km,jc,ic)*amm
                else
                    dxxvz=vz(kp,jc,ic)*app + vz(kc,jc,ic)*acc + vz(km,jc,ic)*amm
                end if

                !   Second derivative in y-direction of vz
                dyyvz=(vz(kc,jp,ic) - 2.0*vz(kc,jc,ic) + vz(kc,jm,ic))*dyq
                !   Second derivative in z-direction of vz
                dzzvz=(vz(kc,jc,ip) - 2.0*vz(kc,jc,ic) + vz(kc,jc,im))*dzq
      
                !   component of grad(pr) along z direction
                dzp=(pr(kc,jc,ic)-pr(kc,jc,im))*dz*al

                !    Calculate right hand side of Eq. 5 (VO96)
                rhsx(kc,jc,ic) = (ga*dq(kc,jc,ic) + ro*ruz(kc,jc,ic) + alre*(dxxvz+dyyvz+dzzvz) - dzp)*dt

                !    Store the non-linear terms for the calculation of the next timestep
                ruz(kc,jc,ic)=dq(kc,jc,ic)

            enddo
        enddo
    enddo

    !$OMP END PARALLEL DO

    call AddBodyIBM(vz,ibm_gz_px,0.0d0)
    call AddBC_Vz
    call AddOutletBC(vz,ren,oxfst,oxfen)
    call SolveImpEqnUpdate_Z

    return

    end
