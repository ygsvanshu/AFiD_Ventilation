!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_CO2.F90                      !
!    CONTAINS: subroutine SolveImpEqnUpdate_CO2           !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for           !
!    carbon dioxide, and updates it to time t+dt          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SolveImpEqnUpdate_CO2

    use param
    use local_arrays, only: co2,rhs
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

    betadx=0.5d0*al*dt/pec

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)

            do kc=1,nx
                ibmx     = 0.5d0*(ibm_body(kc,jc,ic)+ibm_body(kc-1,jc,ic))
                ackl_b   = 1.0d0/(1.0d0-ac3ssk(kc)*betadx)
                amkl(kc) = -am3ssk(kc)*betadx*ackl_b*(1.0d0-ibmx)
                ackl(kc) = 1.0d0
                apkl(kc) = -ap3ssk(kc)*betadx*ackl_b*(1.0d0-ibmx)
                rx1d(kc) = rhs(kc,jc,ic)*ackl_b*(1.0d0-ibmx) - co2(kc,jc,ic)*ibmx
            enddo

            apkl(1)  = amkl(1)  + apkl(1)
            amkl(nx) = amkl(nx) + apkl(nx)

            amkT=amkl(2:nx)
            apkT=apkl(1:nxm)
            ackT=ackl(1:nx)

            call dgttrf(nx,amkT,ackT,apkT,appk,ipkv,info)
            call dgttrs('N',nx,1,amkT,ackT,apkT,appk,ipkv,rx1d,nx,info)
            
            do kc=1,nx
                rhs(kc,jc,ic) = rx1d(kc)
            end do

            do kc=1,nx
                co2(kc,jc,ic) = co2(kc,jc,ic) + rhs(kc,jc,ic)
            end do
        end do
    end do

    return
    end