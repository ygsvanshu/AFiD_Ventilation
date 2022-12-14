!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CreateInitialConditions.F90                    !
!    CONTAINS: subroutine CreateInitialConditions         !
!                                                         !
!    PURPOSE: Initialization routine. Sets initial        !
!     conditions for velocity and temperature             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CreateInitialConditions

    use param
    use local_arrays, only: vy,vx,vz,temp,co2,h2o
    use decomp_2d, only: xstart,xend
    use ventilation_arrays
    use mpih

    implicit none

    integer :: i,j,k
    real    :: xxx,yyy

    ! x -> wall normal
    ! y -> spanwise
    ! z -> streamwise

    do k=1,nxm
        do j=xstart(2),xend(2)
            do i=xstart(3),xend(3)
                vx(k,j,i) = 0.0d0
                vy(k,j,i) = 0.0d0
                vz(k,j,i) = 0.0d0
                temp(k,j,i) = 0.0d0
                co2(k,j,i) = 0.0d0
                h2o(k,j,i) = 0.0d0
            enddo
        enddo
    enddo

    ! if (xend(3).eq.nzm) then
    !     do k=1,nxm
    !         if ((xc(k+1).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
    !             do j=max(1,xstart(2)-lvlhalo),min(nym,xend(2)+lvlhalo)
    !                 xxx = (xm(k) - (oheight-(0.5d0*olen)))/olen
    !                 yyy = ym(j)/ylen
    !                 vz(k,j,nz) = (4.0d0*ivel*dsin(pi*xxx)*dsin(pi*xxx)*dsin(pi*yyy)*dsin(pi*yyy))
    !             enddo
    !         end if
    !     end do
    ! end if

    outvx(:,:)   = 0.0d0
    outvy(:,:)   = 0.0d0
    outvz(:,:)   = 0.0d0
    outtemp(:,:) = 0.0d0
    outco2(:,:)  = 0.0d0
    outh2o(:,:)  = 0.0d0

    ! do k=1,nxm
    !     if ((xc(k+1).gt.(oheight-(0.5d0*olen))).and.(xc(k).lt.(oheight+(0.5d0*olen)))) then
    !         do j=max(1,xstart(2)-lvlhalo),min(nym,xend(2)+lvlhalo)
    !             xxx = (xm(k) - (oheight-(0.5d0*olen)))/olen
    !             yyy = ym(j)/ylen
    !             outvz(k,j) = (4.0d0*ivel*dsin(pi*xxx)*dsin(pi*xxx)*dsin(pi*yyy)*dsin(pi*yyy))
    !         enddo
    !     end if
    ! end do

    ! =ModV17=Vanshu=2020=11=17======================================================
    ! No need to write it twice (waste of disk space)
    ! call HdfWriteRealHalo3D('Results/Initial/temp.h5',temp)
    ! call HdfWriteRealHalo3D('Results/Initial/vx.h5',vx)
    ! call HdfWriteRealHalo3D('Results/Initial/vy.h5',vy)
    ! call HdfWriteRealHalo3D('Results/Initial/vz.h5',vz)
    ! =End=of=ModV17=================================================================

    ! filname = trim('Results/Temp/continua_temp00000000.h5')
    ! call HdfWriteRealHalo3D(filname,temp)
    ! filname = trim('Results/Vx/continua_vx00000000.h5')
    ! call HdfWriteRealHalo3D(filname,vx)
    ! filname = trim('Results/Vy/continua_vy00000000.h5')
    ! call HdfWriteRealHalo3D(filname,vy)
    ! filname = trim('Results/Vz/continua_vz00000000.h5')
    ! call HdfWriteRealHalo3D(filname,vz)

    return

end subroutine CreateInitialConditions