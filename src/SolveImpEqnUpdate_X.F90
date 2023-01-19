!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_X.F90                        !
!    CONTAINS: subroutine SolveImpEqnUpdate_X             !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for velocity  !
!     in any the vertical direction, and updates it to    !
!     time t+dt                                           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SolveImpEqnUpdate_X

    use param
    use local_arrays, only : vx,rhs
    use ibm_arrays, only: ibm_body
    use decomp_2d, only: xstart,xend

    implicit none

    real    :: amkl(nx) ,apkl(nx),ackl(nx)
    real    :: amkT(nxm),ackT(nx),apkT(nxm)
    real    :: appk(nxm-1)
    real    :: rx1d(nx)
    integer :: ipkv(nx)
    integer :: jc,kc,info,ic
    real    :: betadx,ackl_b,ibmx

    betadx=beta*al

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)

            amkl(1) = 0.0d0
            apkl(1) = 0.0d0
            ackl(1) = 1.0d0
            rx1d(1) = 0.0d0 - vx(1,jc,ic)

            do kc=2,nxm
                ibmx     = 0.5d0*(ibm_body(kc,jc,ic)+ibm_body(kc-1,jc,ic))
                ackl_b   = 1.0d0/(1.0d0-ac3ck(kc)*betadx)
                amkl(kc) = -am3ck(kc)*betadx*ackl_b*(1.0d0-ibmx)
                ackl(kc) = 1.0d0
                apkl(kc) = -ap3ck(kc)*betadx*ackl_b*(1.0d0-ibmx)
                rx1d(kc) = rhs(kc,jc,ic)*ackl_b*(1.0d0-ibmx) - vx(kc,jc,ic)*ibmx
            enddo
        
            amkl(nx) = 0.0d0
            apkl(nx) = 0.0d0
            ackl(nx) = 1.0d0
            rx1d(nx) = 0.0d0 - vx(nx,jc,ic)

            amkT = amkl(2:nx)
            apkT = apkl(1:(nx-1))
            ackT = ackl(1:nx)

            call dgttrf(nx,amkT,ackT,apkT,appk,ipkv,info)
            call dgttrs('N',nx,1,amkT,ackT,apkT,appk,ipkv,rx1d,nx,info)
            
            do kc=2,nxm
                rhs(kc,jc,ic) = rx1d(kc)
            end do

            do kc=2,nxm
                vx(kc,jc,ic) = vx(kc,jc,ic) + rhs(kc,jc,ic)
            end do
        end do
    end do

    return
    end
