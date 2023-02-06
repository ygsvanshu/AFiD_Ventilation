!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsVX.F90                            !
!    CONTAINS: subroutine ExplicitTermsVX                 !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the x (vertical) dimension          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExplicitTermsVX

    use param
    use local_arrays, only: vz,vy,vx,temp,co2,h2o,qcap
    use decomp_2d, only: xstart,xend

    implicit none

    integer :: jc,kc
    integer :: km,kp,jmm,jpp,ic,imm,ipp
    real    :: hxx,hxy,hxz
    real    :: udz,udy,tempit
    real    :: fcc,fmm

    udy=dy*0.25
    udz=dz*0.25

    !$OMP   PARALLEL DO &
    !$OMP   DEFAULT(none) &
    !$OMP   SHARED(xstart,xend,nxm,vz,vy,vx,dz,dy) &
    !$OMP   SHARED(am3sk,ac3sk,ap3sk) &
    !$OMP   SHARED(udz,udy,udx3c,qcap,temp) &
    !$OMP   PRIVATE(ic,jc,kc,imm,ipp,km,kp) &
    !$OMP   PRIVATE(jmm,jpp,tempit) &
    !$OMP   PRIVATE(hxx,hxy,hxz)

    do ic=xstart(3),xend(3)
        imm=ic-1
        ipp=ic+1
        do jc=xstart(2),xend(2)
            jmm=jc-1
            jpp=jc+1
            do kc=2,nxm
                km=kc-1
                kp=kc+1

                fcc = dx3c(kc)/g3rc(kc)
                fmm = dx3c(km)/g3rc(kc)
                
                ! vx vz term
                !
                !  d  q_x q_t 
                ! -----------
                !  d   t      
                
                hxz=((((fcc*vz(kc,jc,ipp))+(fmm*vz(km,jc,ipp))) &
                    *(vx(kc,jc,ipp)+vx(kc,jc,ic))) &
                    -(((fcc*vz(kc,jc,ic))+(fmm*vz(km,jc,ic))) &
                    *(vx(kc,jc,ic)+vx(kc,jc,imm))))*udz
                
                ! vx vy term
                !
                !  d  q_x q_r 
                ! -----------
                !  d   r      
                
                hxy=((((fcc*vy(kc,jpp,ic))+(fmm*vy(km,jpp,ic))) &
                    *(vx(kc,jpp,ic)+vx(kc,jc,ic))) &
                    -(((fcc*vy(kc,jc,ic))+(fmm*vy(km,jc,ic))) &
                    *(vx(kc,jc,ic)+vx(kc,jmm,ic))))*udy
                
                ! vx vx term
                !
                !  d  q_x q_x 
                ! -----------
                !  d   x      
                
                hxx=((vx(kp,jc,ic)+vx(kc,jc,ic))*(vx(kp,jc,ic)+vx(kc,jc,ic)) &
                    -(vx(kc,jc,ic)+vx(km,jc,ic))*(vx(kc,jc,ic)+vx(km,jc,ic)) &
                    )*udx3c(kc)*0.25d0

                qcap(kc,jc,ic) = -(hxx+hxy+hxz) + temp(kc,jc,ic) + (lambda_h2o*h2o(kc,jc,ic)) + (lambda_co2*co2(kc,jc,ic))

            enddo
        enddo
    enddo
    
    !$OMP END PARALLEL DO

    return
    
    end
