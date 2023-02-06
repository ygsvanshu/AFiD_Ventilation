!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_Temp.F90                     !
!    CONTAINS: subroutine SolveImpEqnUpdate_Temp          !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for           !
!     temperature, and updates it to time t+dt            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SolveImpEqnUpdate_Temp

    use param
    use local_arrays, only : temp,rhsx,rhsy,rhsz
    use ibm_arrays, only: ibm_gx_px,ibm_gx_py,ibm_gx_pz
    use vent_arrays, only: oxcst,oxcen
    use implicit_decomp
    use decomp_2d

    implicit none

    real    :: amkl(1:nx ),apkl(1:nx) ,ackl(1:nx )
    real    :: amjl(1:nym),apjl(1:nym),acjl(1:nym)
    real    :: amil(1:nzm),apil(1:nzm),acil(1:nzm)

    real    :: amkT(nxm)  ,ackT(nx) ,apkT(nxm)
    real    :: amjT(nym-1),acjT(nym),apjT(nym-1)
    real    :: amiT(nzm-1),aciT(nzm),apiT(nzm-1)

    real    :: appk(nxm-1)
    real    :: appj(nym-2)
    real    :: appi(nzm-2)

    real    :: rx1d(1:nx)
    real    :: ry1d(1:nym)
    real    :: rz1d(1:nzm)

    integer :: ipkv(nx)
    integer :: ipjv(nym)
    integer :: ipiv(nzm)

    integer :: ic,jc,kc,info
    real    :: ackl_b,acjl_b,acil_b
    real    :: betadx

    betadx = 0.5d0*al*dt/pec

    call transpose_x_to_y(rhsx,rhsy,decomp_diff)
    call transpose_y_to_z(rhsy,rhsz,decomp_diff)

    do jc=zstart(2),zend(2)
        do kc=zstart(1),zend(1)
            do ic=1,nzm
                acil_b   = 1.0d0/(1.0d0-ac1gi*betadx)
                if (ibm_gx_pz(kc,jc,ic)) then
                    amil(ic) = 0.0d0
                    acil(ic) = 1.0d0
                    apil(ic) = 0.0d0
                    rz1d(ic) = rhsz(kc,jc,ic)
                else
                    amil(ic) = -am1gi*betadx*acil_b
                    acil(ic) = 1.0d0
                    apil(ic) = -ap1gi*betadx*acil_b
                    rz1d(ic) = rhsz(kc,jc,ic)*acil_b
                end if
            enddo

            ! NEUMANN BOUNDARY CONDITION @ Z = 0
            acil(1) = acil(1) + amil(1)

            if ((kc.ge.oxcst).and.(kc.le.oxcen)) then
                ! RADIATIVE BOUNDARY CONDITION @ Z = ZLEN
                amil(nzm) = -am1oi
            else
                ! NEUMANN BOUNDARY CONDITION @ Z = 0
                acil(nzm) = acil(nzm) + apil(nzm)
            end if

            amiT = amil(2:nzm)
            aciT = acil(1:nzm)
            apiT = apil(1:nzm-1)

            call dgttrf(nzm,amiT,aciT,apiT,appi,ipiv,info)
            call dgttrs('N',nzm,1,amiT,aciT,apiT,appi,ipiv,rz1d,nzm,info)

            do ic=1,nzm                
                rhsz(kc,jc,ic) = rz1d(ic)
            enddo

        end do
    end do

    call transpose_z_to_y(rhsz,rhsy,decomp_diff)

    do ic=ystart(3),yend(3)
        do kc=ystart(1),yend(1)
            do jc=1,nym
                acjl_b   = 1.0d0/(1.0d0-ac2gj*betadx)
                if (ibm_gx_py(kc,jc,ic)) then
                    amjl(jc) = 0.0d0
                    acjl(jc) = 1.0d0
                    apjl(jc) = 0.0d0
                    ry1d(jc) = rhsy(kc,jc,ic)
                else
                    amjl(jc) = -am2gj*betadx*acjl_b
                    acjl(jc) = 1.0d0
                    apjl(jc) = -ap2gj*betadx*acjl_b
                    ry1d(jc) = rhsy(kc,jc,ic)*acjl_b
                end if
            enddo

            ! DIRICHLET BOUNDARY CONDITION @ Y = 0
            acjl(1) = acjl(1) + amjl(1)

            ! DIRICHLET BOUNDARY CONDITION @ Y = YLEN
            acjl(nym) = acjl(nym) + apjl(nym)

            amjT     = amjl(2:nym)
            acjT     = acjl(1:nym)
            apjT     = apjl(1:nym-1)

            call dgttrf(nym,amjT,acjT,apjT,appj,ipjv,info)    
            call dgttrs('N',nym,1,amjT,acjT,apjT,appj,ipjv,ry1d,nym,info)

            do jc=1,nym
                rhsy(kc,jc,ic) = ry1d(jc)
            end do

        end do
    end do

    call transpose_y_to_x(rhsy,rhsx,decomp_diff)

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nx
                ackl_b   = 1.0d0/(1.0d0-ac3ssk(kc)*betadx)
                if (ibm_gx_px(kc,jc,ic)) then
                    amkl(kc) = 0.0d0
                    ackl(kc) = 1.0d0
                    apkl(kc) = 0.0d0
                    rx1d(kc) = rhsx(kc,jc,ic)
                else
                    amkl(kc) = -am3ssk(kc)*betadx*ackl_b
                    ackl(kc) = 1.0d0
                    apkl(kc) = -ap3ssk(kc)*betadx*ackl_b
                    rx1d(kc) = rhsx(kc,jc,ic)*ackl_b
                end if
            enddo

            ! NEUMANN BOUNDARY CONDITION @ X = 0
            apkl(1)  = amkl(1)  + apkl(1)
            
            ! NEUMANN BOUNDARY CONDITION @ X = ALX3
            amkl(nx) = amkl(nx) + apkl(nx)

            amkT = amkl(2:nx)
            ackT = ackl(1:nx)
            apkT = apkl(1:nxm)

            call dgttrf(nx,amkT,ackT,apkT,appk,ipkv,info)
            call dgttrs('N',nx,1,amkT,ackT,apkT,appk,ipkv,rx1d,nx,info)
            
            do kc=1,nx
                rhsx(kc,jc,ic) = rx1d(kc)
            end do
        end do
    end do

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nx
                temp(kc,jc,ic) = temp(kc,jc,ic) + rhsx(kc,jc,ic)
            end do
        end do
    end do

    return
    end
