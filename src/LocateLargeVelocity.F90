!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: LocateLargeDivergence.F90                      !
!    CONTAINS: subroutine LocateLargeDivergence           !
!                                                         ! 
!    PURPOSE: Debugging routine. Output the location(s)   !
!     of excessive divergence.                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine LocateLargeVelocity

    use param
    use local_arrays, only: vy,vx,vz
    use mpih
    use decomp_2d, only: xstart,xend,nrank

    implicit none

    integer :: ic,jc,kc

    if(nrank.eq.0) write(*,*) "!! LARGE VELOCITY !!"
    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                
                if (vx(kc,jc,ic).gt.limitVel.or.vy(kc,jc,ic).gt.limitVel.or.vz(kc,jc,ic).gt.limitVel) then
                    write(6,'(A6,I6,3(A4,I6),3(A6,F10.6))') &
                    " Rank = ",nrank," i = ",ic," j = ",jc," k = ",kc, &
                    " vx = ",vx(kc,jc,ic)," vy = ",vy(kc,jc,ic)," vz = ",vz(kc,jc,ic)
                endif

            enddo
        enddo
    enddo
    
    return     
    end         
