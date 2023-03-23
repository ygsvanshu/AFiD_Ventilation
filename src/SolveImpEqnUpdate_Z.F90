!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_Z.F90                        !
!    CONTAINS: subroutine SolveImpEqnUpdate_Z             !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for velocity  !
!     in any the vertical direction, and updates it to    !
!     time t+dt                                           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SolveImpEqnUpdate_Z

    use param
    use local_arrays, only : vz,rhs,ibm_body
    use decomp_2d, only: xstart,xend

    implicit none

    real    :: amkl(nxm) ,apkl(nxm),ackl(nxm)
    real    :: amkT(nxm-1),ackT(nxm),apkT(nxm-1)
    real    :: appk(nxm-2)
    real    :: rx1d(nxm)
    integer :: ipkv(nxm)
    integer :: jc,kc,info,ic
    real    :: betadx,ackl_b,ibmz
    real    :: visc

    betadx=0.5d0*al*dt/ren

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)

            do kc=1,nxm
                call OutVisc(xm(kc),zc(ic),visc)
                ibmz     = 0.5d0*(ibm_body(kc,jc,ic)+ibm_body(kc,jc,ic-1))
                ackl_b   = 1.0d0/(1.0d0-ac3sk(kc)*betadx*visc)
                amkl(kc) = -am3sk(kc)*betadx*visc*ackl_b*(1.0d0-ibmz)
                ackl(kc) = 1.0d0
                apkl(kc) = -ap3sk(kc)*betadx*visc*ackl_b*(1.0d0-ibmz)
                rx1d(kc) = rhs(kc,jc,ic)*ackl_b*(1.0d0-ibmz) - vz(kc,jc,ic)*ibmz
            enddo

            amkT = amkl(2:nxm)
            ackT = ackl(1:nxm)
            apkT = apkl(1:(nxm-1))
            
            call dgttrf(nxm,amkT,ackT,apkT,appk,ipkv,info)
            call dgttrs('N',nxm,1,amkT,ackT,apkT,appk,ipkv,rx1d,nxm,info)
            
            do kc=1,nxm
                vz(kc,jc,ic) = vz(kc,jc,ic) + rx1d(kc)
            end do
        end do
    end do

    return
    end
