!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_Y.F90                        !
!    CONTAINS: subroutine SolveImpEqnUpdate_Y             !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for velocity  !
!     in any the vertical direction, and updates it to    !
!     time t+dt                                           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SolveImpEqnUpdate_Y

    use param
    use local_arrays, only : vy,rhs,ibm_body
    use decomp_2d, only: xstart,xend

    implicit none

    real    :: amkl(nxm) ,apkl(nxm),ackl(nxm)
    real    :: amkT(nxm-1),ackT(nxm),apkT(nxm-1)
    real    :: appk(nxm-2)
    real    :: rx1d(nxm)
    integer :: ipkv(nxm)
    integer :: jc,kc,info,ic
    real    :: betadx,ackl_b,ibmy
    real    :: visc

    betadx=0.5d0*al*dt/ren

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)

            do kc=1,nxm
                call OutVisc(xm(kc),zm(ic),visc)
                ibmy     = 0.5d0*(ibm_body(kc,jc,ic)+ibm_body(kc,jc-1,ic))
                ackl_b   = 1.0d0/(1.0d0-ac3sk(kc)*betadx*visc)
                amkl(kc) = -am3sk(kc)*betadx*visc*ackl_b*(1.0d0-ibmy)
                ackl(kc) = 1.0d0
                apkl(kc) = -ap3sk(kc)*betadx*visc*ackl_b*(1.0d0-ibmy)
                rx1d(kc) = rhs(kc,jc,ic)*ackl_b*(1.0d0-ibmy) - vy(kc,jc,ic)*ibmy
            enddo

            amkT = amkl(2:nxm)
            ackT = ackl(1:nxm)
            apkT = apkl(1:(nxm-1))
            
            call dgttrf(nxm,amkT,ackT,apkT,appk,ipkv,info)
            call dgttrs('N',nxm,1,amkT,ackT,apkT,appk,ipkv,rx1d,nxm,info)

            do kc=1,nxm
                vy(kc,jc,ic) = vy(kc,jc,ic) + rx1d(kc)
            end do
        end do
    end do

    return
    end
