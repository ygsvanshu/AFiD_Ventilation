!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsCO2.F90                           !
!    CONTAINS: subroutine ExplicitTermsCO2                !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the carbon dioxide.                                 !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExplicitTermsCO2

    use param
    use local_arrays, only: vy,vx,vz,co2,qco2
    use decomp_2d, only: xstart,xend

    implicit none

    integer :: jc,kc,ic
    integer :: km,kp,jm,jp,im,ip
    real    :: htx,hty,htz,udy,udz
    real    :: udzq,udyq
    real    :: dyyt,dzzt
    real    :: fcc,fmm

    udz=dz*0.25
    udy=dy*0.25
    udzq=dzq/pec
    udyq=dyq/pec

    !$OMP   PARALLEL DO &
    !$OMP   DEFAULT(none) &
    !$OMP   SHARED(xstart,xend,vz,vy,vx,nxm) &
    !$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk,udz) &
    !$OMP   SHARED(udy,udzq,udyq,udx3c,co2,qco2) &
    !$OMP   PRIVATE(ic,jc,kc,im,ip,km,kp,jm,jp) &
    !$OMP   PRIVATE(htx,hty,htz,dyyt,dzzt)

    do ic=xstart(3),xend(3)
        im=ic-1
        ip=ic+1
        do jc=xstart(2),xend(2)
            jm=jc-1
            jp=jc+1
            do kc=2,nxm
                km=kc-1
                kp=kc+1

                fcc = dx3c(kc)/g3rc(kc)    
                fmm = dx3c(km)/g3rc(kc)

                htz = (((fmm*vz(km,jc,ip))+(fcc*vz(kc,jc,ip)))*(co2(kc,jc,ip)+co2(kc,jc,ic))- &
                    ((fmm*vz(km,jc,ic))+(fcc*vz(kc,jc,ic)))*(co2(kc,jc,ic)+co2(kc,jc,im)) &
                    )*udz

                hty = (((fcc*vy(kc,jp,ic))+(fmm*vy(km,jp,ic)))*(co2(kc,jp,ic)+co2(kc,jc,ic))- &
                    ((fcc*vy(kc,jc,ic))+(fmm*vy(km,jc,ic)))*(co2(kc,jc,ic)+co2(kc,jm,ic)) &
                    )*udy

                htx = ((vx(kp,jc,ic)+vx(kc,jc,ic))*(co2(kp,jc,ic)+co2(kc,jc,ic))- &
                    (vx(kc,jc,ic)+vx(km,jc,ic))*(co2(kc,jc,ic)+co2(km,jc,ic)) &
                    )*udx3c(kc)*0.25d0

                dzzt = (co2(kc,jc,ip) -2.0*co2(kc,jc,ic) + co2(kc,jc,im))*udzq
                dyyt = (co2(kc,jp,ic) -2.0*co2(kc,jc,ic) + co2(kc,jm,ic))*udyq

                qco2(kc,jc,ic) = -(htx+hty+htz)+dyyt+dzzt

            enddo
            
            dzzt = (co2(1,jc,ip) -2.0*co2(1,jc,ic) + co2(1,jc,im))*udzq
            dyyt = (co2(1,jp,ic) -2.0*co2(1,jc,ic) + co2(1,jm,ic))*udyq

            qco2(1,jc,ic) = dyyt+dzzt

            dzzt = (co2(nx,jc,ip) -2.0*co2(nx,jc,ic) + co2(nx,jc,im))*udzq
            dyyt = (co2(nx,jp,ic) -2.0*co2(nx,jc,ic) + co2(nx,jm,ic))*udyq

            qco2(nx,jc,ic) = dyyt+dzzt

        enddo
    enddo

    !$OMP  END PARALLEL DO

    return

end