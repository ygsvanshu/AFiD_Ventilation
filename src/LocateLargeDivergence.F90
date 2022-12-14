!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: LocateLargeDivergence.F90                      !
!    CONTAINS: subroutine LocateLargeDivergence           !
!                                                         ! 
!    PURPOSE: Debugging routine. Output the location(s)   !
!     of excessive divergence.                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine LocateLargeDivergence

    use param
    use local_arrays, only: vy,vx,vz
    use mpih
    use decomp_2d, only: xstart,xend,nrank

    implicit none

    integer :: jc,kc,kp,jp,ic,ip
    real    :: dvxdx,dvydy,dvzdz,dqcap
        
    if(nrank.eq.0) write(*,*) "!! LARGE DIVERGENCE !!"
    do ic=xstart(3),xend(3)
        ip=ic+1
        do jc=xstart(2),xend(2)
            jp=jc+1
            do kc=1,nxm
                kp=kc+1

                dvxdx = (vx(kp,jc,ic)-vx(kc,jc,ic))*udx3m(kc)
                dvydy = (vy(kc,jp,ic)-vy(kc,jc,ic))*dy
                dvzdz = (vz(kc,jc,ip)-vz(kc,jc,ic))*dz
                dqcap = dvxdx + dvydy + dvzdz

                if (abs(dqcap).gt.resid) then
                    write(6,'(A6,I6,3(A4,I6),4(A10,F10.6))') &
                    " Rank = ",nrank," i = ",ic," j = ",jc," k = ",kc, &
                    " dvx/dx = ",dvxdx," dvy/dy = ",dvydy," dvz/dz = ",dvzdz," Div = ",dqcap
                endif

            enddo
        enddo
    enddo
    
    return     
    end         
