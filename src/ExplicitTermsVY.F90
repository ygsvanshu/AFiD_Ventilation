!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: ExplicitTermsVY.F90                            !
!    CONTAINS: subroutine ExplicitTermsVY                 !
!                                                         !
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the y (horizontal) dimension        !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExplicitTermsVY

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vx,vy,vz,dph
    use vent_arrays, only: outvscy


    implicit none

    integer :: kc,kp,jpp,jmm,jc,ic,imm,ipp
    integer :: kpp,kmm
    real    :: udzq,udyq
    real    :: udy,udz,hyx,hyy,hyz
    real    :: dyyvy, dzzvy
    real    :: fpp,fpc,fmc,fmm

    udyq=dyq/ren
    udzq=dzq/ren

    udy=dy*0.25
    udz=dz*0.25

    !$OMP  PARALLEL DO &
    !$OMP  DEFAULT(none) &
    !$OMP  SHARED(xstart,xend,vz,vy,vx,dz,dy) &
    !$OMP  SHARED(kmv,kpv,am3sk,ac3sk,ap3sk,udz) &
    !$OMP  SHARED(udy,udzq,udyq,udx3m,dph,nxm) &
    !$OMP  PRIVATE(ic,jc,kc,imm,ipp,kmm,kp,kpp) &
    !$OMP  PRIVATE(jmm,jpp) &
    !$OMP  PRIVATE(hyx,hyy,hyz,dyyvy,dzzvy)

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

                hyx=((vx(kp,jc,ic)+vx(kp,jmm,ic))*((fpp*vy(kpp,jc,ic))+(fpc*vy(kc,jc,ic))) &
                    -(vx(kc,jc,ic)+vx(kc,jmm,ic))*((fmc*vy(kc,jc,ic))+(fmm*vy(kmm,jc,ic))) &
                    )*udx3m(kc)*0.25d0
            
                hyy=( (vy(kc,jpp,ic)+vy(kc,jc,ic)) &
                    *(vy(kc,jpp,ic)+vy(kc,jc,ic)) &
                    -(vy(kc,jmm,ic)+vy(kc,jc,ic)) &
                    *(vy(kc,jmm,ic)+vy(kc,jc,ic)) &
                    )*udy

                hyz=( (vy(kc,jc,ipp)+vy(kc,jc,ic)) &
                    *(vz(kc,jc,ipp)+vz(kc,jmm,ipp)) &
                    -(vy(kc,jc,ic)+vy(kc,jc,imm)) &
                    *(vz(kc,jc,ic)+vz(kc,jmm,ic)) &
                    )*udz

                dyyvy=(vy(kc,jpp,ic) -2.0*vy(kc,jc,ic) + vy(kc,jmm,ic))*udyq*outvscy(kc,ic)
                dzzvy=(vy(kc,jc,ipp) -2.0*vy(kc,jc,ic) + vy(kc,jc,imm))*udzq*outvscy(kc,ic)

                dph(kc,jc,ic)=-(hyx+hyy+hyz)+dyyvy+dzzvy

            enddo
        enddo
    enddo

    !$OMP END PARALLEL DO

    return
    
    end