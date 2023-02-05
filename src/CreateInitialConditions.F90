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

    return

end subroutine CreateInitialConditions