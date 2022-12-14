!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CreateResultDirectories.F90                    !
!    CONTAINS: subroutine CreateResultDirectories         !
!                                                         !
!    PURPOSE: Create directories for saving results       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Created on 05/02/2020 by Vanshu [Modification #14]

subroutine CreateResultDirectories

    implicit none

    call system('mkdir -p Results')
! =ModV17=Vanshu=2020=11=17======================================================
    ! call system('mkdir -p Results/Initial')
! =End=of=ModV17=================================================================
    call system('mkdir -p Results/Vx')
    call system('mkdir -p Results/Vy')
    call system('mkdir -p Results/Vz')
    call system('mkdir -p Results/Temp')
    call system('mkdir -p Results/CO2')
    call system('mkdir -p Results/H2O')
! ==============================================================================
!  ModR04 Robert 2020-08-05
!     Write statistics snapshots
!  ModR05 Robert 2020-09-09
!     Write Gridinfo snapshots
! ==============================================================================
    call system('mkdir -p Results/Stats')
    call system('mkdir -p Results/Grid')
! ==============================================================================
!  End of ModR04 & ModR05
! ==============================================================================

end subroutine CreateResultDirectories

!! End [Modification #14]
