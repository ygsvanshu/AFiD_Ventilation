!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: ExplicitTermsVZ.F90                            !
!    CONTAINS: subroutine ExplicitTermsVZ                 !
!                                                         !
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the z (horizontal) dimension        !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExplicitTermsVZ

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vx,vy,vz,dq
    use vent_arrays, only: outvscz

    implicit none

    integer :: kc,kp,jpp,jmm,jc,ic,imm,ipp
    integer :: kmm,kpp
    real    :: hzx,hzy,hzz,udy,udz
    real    :: udyq,udzq
    real    :: dzzvz,dyyvz
    real    :: fpp,fpc,fmc,fmm

    udyq=dyq/ren
    udzq=dzq/ren

    udy=dy*0.25
    udz=dz*0.25

    !$OMP  PARALLEL DO &
    !$OMP  DEFAULT(none) &
    !$OMP  SHARED(xstart,xend,vz,vy,vx,dz,dy,udx3m) &
    !$OMP  SHARED(kmv,kpv,am3sk,ac3sk,ap3sk,udz) &
    !$OMP  SHARED(udy,udzq,udyq,dq,nxm) &
    !$OMP  PRIVATE(ic,jc,kc,imm,ipp,kmm,kp,kpp) &
    !$OMP  PRIVATE(jmm,jpp) &
    !$OMP  PRIVATE(hzz,hzy,hzx,dzzvz,dyyvz)

    do ic=xstart(3),xend(3)
        imm=ic-1
        ipp=ic+1
        do jc=xstart(2),xend(2)
            jmm=jc-1
            jpp=jc+1
            do kc=1,nxm
                kmm=kmv(kc)
                kpp=kpv(kc)
                kp=kc+1

                fpp = dx3c(kpp)/g3rc(kpp)  !! GSY
                fpc = dx3c(kc)/g3rc(kpp)   !! GSY
                fmc = dx3c(kc)/g3rc(kc)    !! GSY
                fmm = dx3c(kmm)/g3rc(kc)   !! GSY

                hzx=((vx(kp,jc,ic)+vx(kp,jc,imm))*((fpp*vz(kpp,jc,ic))+(fpc*vz(kc,jc,ic))) &  !! GSY
                    -(vx(kc,jc,ic)+vx(kc,jc,imm))*((fmc*vz(kc,jc,ic))+(fmm*vz(kmm,jc,ic))) &  !! GSY
                    )*udx3m(kc)*0.25d0

                hzy=( (vy(kc,jpp,ic)+vy(kc,jpp,imm)) &
                    *(vz(kc,jpp,ic)+vz(kc,jc,ic)) &
                    -(vy(kc,jc,ic)+vy(kc,jc,imm)) &
                    *(vz(kc,jc,ic)+vz(kc,jmm,ic)) &
                    )*udy

                hzz=( (vz(kc,jc,ipp)+vz(kc,jc,ic)) &
                    *(vz(kc,jc,ipp)+vz(kc,jc,ic)) &
                    -(vz(kc,jc,imm)+vz(kc,jc,ic)) &
                    *(vz(kc,jc,imm)+vz(kc,jc,ic)) &
                    )*udz
                
                dzzvz=(vz(kc,jc,ipp) -2.0*vz(kc,jc,ic) + vz(kc,jc,imm))*udzq*outvscz(kc,ic)
                dyyvz=(vz(kc,jpp,ic) -2.0*vz(kc,jc,ic) + vz(kc,jmm,ic))*udyq*outvscz(kc,ic)

                dq(kc,jc,ic)=-(hzx+hzy+hzz)+dyyvz+dzzvz

            enddo
        enddo
    enddo

    !$OMP END PARALLEL DO

    return

    end