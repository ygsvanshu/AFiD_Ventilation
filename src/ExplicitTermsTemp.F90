!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsTemp.F90                          !
!    CONTAINS: subroutine ExplicitTermsTemp               !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the temperature.                                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExplicitTermsTemp

    use param
    use local_arrays, only: vx,vy,vz,temp,hro
    use decomp_2d, only: xstart,xend

    implicit none

    integer :: km,kc,kp,jm,jc,jp,im,ic,ip
    real    :: htx,hty,htz,udy,udz
    real    :: fcc,fmm

    udz=dz*0.25
    udy=dy*0.25

    !$OMP   PARALLEL DO &
    !$OMP   DEFAULT(none) &
    !$OMP   SHARED(xstart,xend,vz,vy,vx,nxm) &
    !$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk,udz) &
    !$OMP   SHARED(udy,udzq,udyq,udx3c,temp,hro) &
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
                !
                !    rho vz term
                !
                !                d  rho q_z
                !             -----------
                !                d   z      
                !
                htz = (((fmm*vz(km,jc,ip))+(fcc*vz(kc,jc,ip)))*(temp(kc,jc,ip)+temp(kc,jc,ic))- &
                        ((fmm*vz(km,jc,ic))+(fcc*vz(kc,jc,ic)))*(temp(kc,jc,ic)+temp(kc,jc,im)) &
                        )*udz
                !
                !    rho vy term
                !
                !                d  rho q_y 
                !             -----------
                !                d   y      
                !
                hty = (((fcc*vy(kc,jp,ic))+(fmm*vy(km,jp,ic)))*(temp(kc,jp,ic)+temp(kc,jc,ic))- &
                        ((fcc*vy(kc,jc,ic))+(fmm*vy(km,jc,ic)))*(temp(kc,jc,ic)+temp(kc,jm,ic)) &
                        )*udy
                !
                !    rho vx term
                !
                !                 d  rho q_x 
                !                -----------
                !                 d   x      
                !
                htx = ((vx(kp,jc,ic)+vx(kc,jc,ic))*(temp(kp,jc,ic)+temp(kc,jc,ic))- &
                        (vx(kc,jc,ic)+vx(km,jc,ic))*(temp(kc,jc,ic)+temp(km,jc,ic)) &
                        )*udx3c(kc)*0.25d0

                hro(kc,jc,ic) = -(htx+hty+htz)
                
            enddo
        enddo
    enddo

    !$OMP  END PARALLEL DO

    return
    
    end