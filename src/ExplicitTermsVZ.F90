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
      use local_arrays, only: vx,vy,vz,dq
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: kc,kp,jpp,jmm,jc,ic,imm,ipp
      integer :: kmm,kpp
      real    :: hzx,hzy,hzz,udy,udz
      real    :: udyq,udzq
      real    :: dzzvz,dyyvz
! =ModV21=Vanshu=2022=02=21======================================================
      real    :: fpp,fpc,fmc,fmm
! =End=of=ModV21=================================================================

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

! =ModV21=Vanshu=2022=02=21======================================================
         fpp = dx3c(kpp)/g3rc(kpp)  !! GSY
         fpc = dx3c(kc)/g3rc(kpp)   !! GSY
         fmc = dx3c(kc)/g3rc(kc)    !! GSY
         fmm = dx3c(kmm)/g3rc(kc)   !! GSY
!
!     vz vx term
!
!                 d  q_t q_x 
!                -----------
!                 d   x      
!
      hzx=((vx(kp,jc,ic)+vx(kp,jc,imm))*((fpp*vz(kpp,jc,ic))+(fpc*vz(kc,jc,ic))) &  !! GSY
          -(vx(kc,jc,ic)+vx(kc,jc,imm))*((fmc*vz(kc,jc,ic))+(fmm*vz(kmm,jc,ic))) &  !! GSY
          )*udx3m(kc)*0.25d0
!
!     vz vy term
!
!                 d  q_t q_r 
!                ------------
!                 d   r      
!
      hzy=( (vy(kc,jpp,ic)+vy(kc,jpp,imm)) &
           *(vz(kc,jpp,ic)+vz(kc,jc,ic)) &
           -(vy(kc,jc,ic)+vy(kc,jc,imm)) &
           *(vz(kc,jc,ic)+vz(kc,jmm,ic)) &
          )*udy
!      
!     vz vz term
!
!                 d  q_t q_t 
!                ------------
!                 d   t      
!
      hzz=( (vz(kc,jc,ipp)+vz(kc,jc,ic)) &
           *(vz(kc,jc,ipp)+vz(kc,jc,ic)) &
           -(vz(kc,jc,imm)+vz(kc,jc,ic)) &
           *(vz(kc,jc,imm)+vz(kc,jc,ic)) &
          )*udz

! =End=of=ModV21=================================================================

!
!   11 second derivative of vz
!
            dzzvz=(vz(kc,jc,ipp) &
                  -2.0*vz(kc,jc,ic) &
                  +vz(kc,jc,imm))*udzq
!
!   22 second derivative of vz
!
            dyyvz=(vz(kc,jpp,ic) &
                  -2.0*vz(kc,jc,ic) &
                  +vz(kc,jmm,ic))*udyq

!   Sum up ExplicitTermsVZ 
!
            dq(kc,jc,ic)=-(hzx+hzy+hzz)+dyyvz+dzzvz
! ==============================================================================
!  End of ModR02
! ==============================================================================
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      return
      end
!
