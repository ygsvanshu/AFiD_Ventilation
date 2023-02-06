!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: QuitRoutine.F90                                !
!    CONTAINS: subroutine QuitRoutine, NotifyError        !
!                                                         ! 
!    PURPOSE: Routines to exit the program and write the  !
!     data if necessary                                   !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine QuitRoutine(tin,normalexit,errorcode)

    use hdf5
    use mpih
    use param
    use decomp_2d, only: nrank, decomp_2d_finalize
    use decomp_2d_fft
    use implicit_decomp

    implicit none
    logical, intent(in) :: normalexit
    integer :: errorcode
    real :: tin(3)

    if (ismaster) write(6,*) ''

    if(errorcode.ne.100) then !EP skip if already finalized

        tin(3) = MPI_WTIME()
        if(ismaster) then
            call NotifyError(errorcode) 
        endif

        if(normalexit) then
            if(nrank.eq.0) write(6,'(a,f10.2,a)') 'Total Iteration Time = ',tin(3) -tin(2),' sec.'
            ! ==============================================================================
            !  ModR04 Robert 2020-08-05
            !  Write statistics snapshots: WriteStats -> WriteStatsEnd
            ! ==============================================================================
            ! ==============================================================================
            !  ModV19 Vanshu 2022-01-20
            !  Write statistics snapshots only if time greater than tsta
            !  No need to write flow field snapshot as it is already written by 
            !  WriteFlowField subroutine.
            ! ==============================================================================
            if (statcalc.and.(time.gt.tsta)) then
                ! call WriteStatsSnap
                call WriteStatsEnd
            endif
            ! call WriteFlowFieldSnapshot
            ! ==============================================================================
            !  End of ModV19
            ! ==============================================================================
            ! ==============================================================================
            !  End of ModR04
            ! ==============================================================================
            call WriteFlowField
            call WriteOutlet
        else
            call Movie_xcut(errorcode)
            call Movie_ycut(errorcode)
            call Movie_zcut(errorcode)
            call Movie_outlet(errorcode)
            call MPI_Abort(MPI_COMM_WORLD,1)
        endif
        
        call DeallocateVariables
        call HdfClose
        call decomp_info_finalize(decomp_diff)
        call decomp_2d_fft_finalize
        call decomp_2d_finalize

    endif

end subroutine QuitRoutine

subroutine NotifyError(errorcode)

    use param
    
    implicit none
    integer, intent(in) :: errorcode

    if(errorcode.eq.166) then 
        write(6,*) 'dt too small, DT= ', dt 
        168 format(10x,e14.7)
    else if(errorcode.eq.165) then
        write(6,164) 
        164 format(10x,'cfl too large  ')
    else if(errorcode.eq.266) then
        write(6,268)
        268 format(10x,'velocities diverged')
        call LocateLargeVelocity
    else if(errorcode.eq.169) then
        write(6,178) 
        write(6,179) 
        write(6,180)               
        178 format(10x,'too large local residue for mass conservation at:')
        179 format(10x,'Probably the matrix in SolvePressureCorrection becomes singular')
        180 format(10x,'Try changing nxm or str3')
        call LocateLargeDivergence
    else if(errorcode.eq.333) then
        write(*,*) "time greater than tmax"
        write(*,*) "statistics and continuation updated"
    else if(errorcode.eq.334) then
        write(*,*) "walltime greater than walltimemax"
        write(*,*) "statistics and continuation updated"
    else if(errorcode.eq.444) then
        write(*,*) "FFT size in ny or nz is not efficient"
    else if(errorcode.eq.555) then
        write(*,*) "Maximum number of timesteps reached"
        write(*,*) "statistics and continuation updated"
    else if(errorcode.eq.666) then
        write(*,*) "Person geomtery is not fully inside domain"
    else if(errorcode.eq.667) then
        write(*,*) "Breath position is not fully inside domain"
    else 
        write(*,*) "!Unknown Error!"
    end if

    return

    end

     
