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
    call system('mkdir -p Results/Vx')
    call system('mkdir -p Results/Vy')
    call system('mkdir -p Results/Vz')
    call system('mkdir -p Results/Pr')
    call system('mkdir -p Results/Temp')
    call system('mkdir -p Results/CO2')
    call system('mkdir -p Results/H2O')
    call system('mkdir -p Results/Grid')

end subroutine CreateResultDirectories

!! End [Modification #14]
