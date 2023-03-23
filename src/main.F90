program AFiD

    use mpih
    use param
    use local_arrays, only: vx,vy,vz,temp,co2,h2o,pr,ibm_body
    use vent_arrays
    use decomp_2d
    use decomp_2d_fft
    use stat_arrays, only: nstatsamples,tstat,tinterval
        
    implicit none
        
    integer :: errorcode, nthreads, ierror
    real    :: instCFL,dmax
    real    :: ti(2), tin(3), minwtdt
    real   :: ts
    integer :: prow=0,pcol=0
    integer :: lfactor,lfactor2
    character(100) :: arg
    
    real  :: progress_t(2)
    real   :: progress_runtime
    real   :: progress_eta

    integer :: progress_h     !! Hours
    integer :: progress_m     !! Minutes
    integer :: progress_s     !! Seconds

    integer :: progress_eta_h   !! Hours
    integer :: progress_eta_m   !! Minutes
    integer :: progress_eta_s   !! Seconds

    !*******************************************************
    !******* Read input file bou.in by all processes********
    !*******************************************************

    call ReadInputFile

    if (command_argument_count().eq.2) then
        call get_command_argument(1,arg)
        read(arg,'(i10)') prow
        call get_command_argument(2,arg)
        read(arg,'(i10)') pcol
    endif

    call decomp_2d_init(nxm,nym,nzm,prow,pcol,(/ .false.,.false.,.false. /))

    ts=MPI_WTIME()
    tin(1) = MPI_WTIME()

    call MpiBarrier

    call HdfStart

    if (nrank.eq.master) ismaster = .true.

    if (ismaster) write(6,*) 'MPI tasks=', nproc

    if (ismaster) then

    ! "Nothing wrong with a little bit of artistic liberty :3 " - Vanshu

    write(6,'(A77)') '============================================================================='
    write(6,'(A77)') '                           __________ _     _  _____                         '
    write(6,'(A77)') '                          /    _____/´ ) _ / \|  __ \                        '
    write(6,'(A77)') '                         / /| |      /__(_)\_/| |  \ \                       '
    write(6,'(A77)') '========================/ /=| |===============| |===\ \======================'
    write(6,'(A77)') '=======================/ ___   __/======| |===| |===/ /======================'
    write(6,'(A77)') '                      / /   | |         | |   | |__/ /                       '
    write(6,'(A77)') '                     /_/    |_|         |_|   |_____/                        '
    write(6,'(A77)') ' AFiD 2.0                                                                    '
    write(6,'(A77)') '============================================================================='
    write(6,'(A77)') '          _          _     _  _       _   _      _           _   _           '
    write(6,'(A77)') '         |_) |_| \/ (_` | /  (_`     / \ |_     |_ |  | | | | \ (_`          '
    write(6,'(A77)') '         |   | | /  ._) | \_ ._)     \_/ |      |  |_ |_| | |_/ ._)          '
    write(6,'(A77)') ''
    write(6,'(A77)') '============================================================================='
    write(6,'(A77)') ''
    write(6,'(A77)') '                            V E N T I L A T I O N                            '
    write(6,'(A77)') ''
    write(6,'(A77)') '**** 3D CELL ASPECT RATIOS ============================================= ****'
    write(6,'(A77)') ''
    write(6,'(A40,F8.4)') '     Domain size in Y / Height        = ', ylen/alx3
    write(6,'(A40,F8.4)') '     Domain size in Z / Height        = ', zlen/alx3
    write(6,'(A77)') ''
    write(6,'(A77)') '**** PHYSICAL PARAMETERS =============================================== ****'
    write(6,'(A77)') ''
    write(6,'(A77)') '     Case: Periodic lateral wall boundary condition                          '
    write(6,'(A40,E10.3)') '     Rayleigh Number                  = ', ray
    write(6,'(A40,E10.3)') '     Prandtl Number                   = ', pra
    write(6,'(A77)') ''
    write(6,'(A77)') '**** VENTILATION VERSION =============================================== ****'
    write(6,'(A77)') ''
    write(6,'(A40,E10.3)') '     Expansion ratio for CO2          = ', lambda_co2
    write(6,'(A40,F8.4)')  '     Expansion ratio for H2O          = ', lambda_h2o
    write(6,'(A40,F8.4)')  '     Velocity imposed at inlet vent   = ', ivel
    write(6,'(A40,F8.4)')  '     Inlet vent dimension             = ', ilen
    write(6,'(A40,F8.4)')  '     Inlet vent position              = ', iheight
    write(6,'(A40,F8.4)')  '     Outlet vent dimension            = ', olen
    write(6,'(A40,F8.4)')  '     Outlet vent position             = ', oheight
    write(6,'(A77)') ''
    write(6,'(A77)') '**** TIME-STEPPING ===================================================== ****'
    write(6,'(A77)') ''
    if(variabletstep) then
        write(6,'(A40,F8.4)')  '     Variable dt, Fixed CFL           = ', limitCFL
        write(6,'(A40,E10.3)') '     Variable dt, maximum dt          = ', dtmax
        write(6,'(A40,E10.3)') '     Variable dt, minimum dt          = ', dtmin   
    else 
        write(6,'(A40,F8.4)')  '     Fixed dt, dt                     = ', dtmax
        write(6,'(A40,F8.4)')  '     Fixed dt, maximum CFL            = ', limitCFL          
    endif
    write(6,'(A40,I15)')       '     Maximum number of timesteps      = ', ntst

    call CreateResultDirectories

    endif

    call InitTimeMarchScheme
    call InitVariables
    call CreateGrid
    call WriteGridInfo

    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.true.,.false./), comm_zcut, ierror)
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.false.,.true./), comm_ycut, ierror)
    comm_xcut = MPI_COMM_WORLD

    if (ismaster) then
        write(6,'(A77)') ''
        write(6,'(A77)') '**** GRID RESOLUTION =================================================== ****'
        write(6,'(A77)') ''
        write(6,'(A40,I8)') '     Number of grid points in X       = ', nxm
        write(6,'(A40,I8)') '     Number of grid points in Y       = ', nym
        write(6,'(A40,I8)') '     Number of grid points in Z       = ', nzm
        write(6,'(A77)') ''
        write(6,'(A40,E10.3)') '     Average cell spacing in X        = ', 1.d0/dx
        write(6,'(A40,E10.3)') '     Uniform cell spacing in Y        = ', 1.d0/dy
        write(6,'(A40,E10.3)') '     Average cell spacing in Z        = ', 1.d0/dz
        write(6,'(A77)') ''
    endif

    time=0.d0

    call InitPressureSolver
    call InitVents  ! Initialize vent indices and cell areas

    if(readflow) then
        if(ismaster) write(6,*) 'Reading initial condition from file'
        call ReadFlowField
        call ReadOutlet 
    else
        if(ismaster) write(6,*) 'Creating initial condition'
        ntime=0
        time=0.d0
        instCFL=0.d0
        if (person_on) call CreateBodyIBM
        call CreateInitialConditions
        call SetOutletBC
    endif

    call update_halo(vx,lvlhalo)
    call update_halo(vy,lvlhalo)
    call update_halo(vz,lvlhalo)
    call update_halo(temp,lvlhalo)
    call update_halo(co2,lvlhalo)
    call update_halo(h2o,lvlhalo)
    call update_halo(pr,lvlhalo)
    call update_halo(ibm_body,lvlhalo)
  
    call CorrectOutletFlux
    ! Use the computed flux and start time to slowly ramp up 
    ! the inlet velocity. This is taken care in SetInletBC
    isvel = iflux/iarea
    tsvel = time
    call SetInletBC
    call SetWallBCs

    ! Init tstat
    if(statcalc) then
        tsta = max(tsta,time)
        call InitStats
        if (statread) then
            call ReadStats
            write(6,*) 'Reading previous stats'
        end if
        if (ismaster) write(6,*) 'Averaging starts at t =',tsta
    end if

    if (savemovie) call InitMovie

    call CheckDivergence(dmax)
    call GlobalQuantities

    if(ismaster) then   
        write(6,*)'Initial maximum divergence: ',dmax               
        tin(2) = MPI_WTIME()             
        write(6,'(a,f6.2,a)') ' Initialization Time = ', tin(2) -tin(1), ' sec.'
        write(6,*) ''
        write(6,'(A)') '-----------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+-----------+------------+------------'
        write(6,'(A)') '  Timestep |  FlowTime  |     DT     |  Max Div   |  vmax(1)   |  vmax(2)   |  vmax(3)   |  tempmean  |  tempmax   |  tempmin   |  CO2 mean  |  H2O mean  | CPU Time  |   CPU DT   |     ETA    '
        write(6,'(I10,A,11(E10.3,A),I3.2,2(A,I2.2),A,E10.3,A,I3.2,2(A,I2.2))') -1,' | ',time,' | ',dt,' | ',dmax,' | ',vmax(1),' | ',vmax(2),' | ',vmax(3),' | ',tempm,' | ',tempmax,' | ',tempmin,' | ',co2m,' | ',h2om,' | ',0,':',0,':',0,' | ',0.d0 ,' | ',0,':',0,':',0
    end if
                                                        
    !  ********* starts the time dependent calculation ***
    errorcode = 0 !EP set errocode to 0 (OK)
    minwtdt = huge(0.0d0) !EP initialize minimum time step walltime

    ! Check input for efficient FFT
    ! factorize input FFT directions. The largest factor should
    ! be relatively small to have efficient FFT's
    lfactor=2 ! initialize
    call Factorize(nym,lfactor2) ! check nym
    lfactor=max(lfactor,lfactor2)
    call Factorize(nzm,lfactor2)
    lfactor=max(lfactor,lfactor2)
    ! if largest factor larger than 7 quit the simulation     
    if (lfactor>7) errorcode=444

    !! Modified on 11/02/2020 by Vanshu [Modification #17]
    progress_t(1) = MPI_WTIME()    
    !! End [Modification #17]

    do ntime=0,ntst      

        ti(1) = MPI_WTIME()

        !EP   Determine timestep size
        call CalcMaxCFL(instCFL)

        if (variabletstep) then
            if (ntime.gt.0) then
                if (instCFL.lt.1.0d-8) then !EP prevent fp-overflow
                    dt=dtmax
                else
                    dt=limitCFL/instCFL
                endif
            else
                if (readflow) then
                    dt=limitCFL/instCFL
                else
                    dt=dtmin
                end if
            endif
            if(dt.gt.dtmax) dt=dtmax
            if(dt.lt.dtmin) errorcode = 166
        else  
            !RO! fixed time-step
            instCFL=instCFL*dt
            if(instCFL.gt.limitCFL) errorcode = 165
        endif

        call TimeMarcher

        if((ntime.eq.1).or.(mod(ntime,nout).eq.0)) then ! TOUT -> NOUT
            
            call GlobalQuantities
            if(vmax(1).gt.limitVel.and.vmax(2).gt.limitVel) errorcode = 266

            call CalcMaxCFL(instCFL)
            call CheckDivergence(dmax)

            if(.not.variabletstep) instCFL=instCFL*dt
            if(abs(dmax).gt.resid) errorcode = 169

            if((statcalc).and.(time.gt.tsta)) call CalcStats

        endif

        if(time.gt.tmax) errorcode = 333

        ti(2) = MPI_WTIME()
        minwtdt = min(minwtdt,ti(2) - ti(1))

        if(mod(ntime,nout).eq.0) then ! TOUT -> NOUT
            if(ismaster) then

                !! Modified on 11/02/2020 by Vanshu [Modification #17]

                progress_t(2) = MPI_WTIME()
                progress_runtime = progress_t(2)-progress_t(1)
                progress_h = int(progress_runtime/3600)
                progress_m = int((progress_runtime-(3600*progress_h))/60)
                progress_s = int(progress_runtime-(3600*progress_h)-(60*progress_m))
                if (time.eq.0) then
                    progress_eta = walltimemax
                else
                    progress_eta = min( (((tmax-time)/time)*progress_runtime) , ((real(ntst-ntime)/real(ntime))*progress_runtime) )
                end if
                progress_eta_h = int(progress_eta/3600)
                progress_eta_m = int((progress_eta-(3600*progress_eta_h))/60)
                progress_eta_s = int(progress_eta-(3600*progress_eta_h)-(60*progress_eta_m))

                write(6,'(A)') '-----------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+-----------+------------+------------'
                write(6,'(A)') '  Timestep |  FlowTime  |     DT     |  Max Div   |  vmax(1)   |  vmax(2)   |  vmax(3)   |  tempmean  |  tempmax   |  tempmin   |  CO2 mean  |  H2O mean  | CPU Time  |   CPU DT   |     ETA    '
                write(6,'(I10,A,11(E10.3,A),I3.2,2(A,I2.2),A,E10.3,A,I3.2,2(A,I2.2))') ntime,' | ',time,' | ',dt,' | ',dmax,' | ',vmax(1),' | ',vmax(2),' | ',vmax(3),' | ',tempm,' | ',tempmax,' | ',tempmin,' | ',co2m,' | ',h2om,' | ',progress_h,':',progress_m,':',progress_s,' | ',minwtdt ,' | ',progress_eta_h,':',progress_eta_m,':',progress_eta_s
            
            end if
            minwtdt = huge(0.0d0)
        end if

        if (savesnap) then
            if (time.ge.startsnap) then
                if (mod(time,tsnap).lt.dt) then
                    call WriteFlowFieldSnapshot
                    if(statcalc .and. (time.gt.tsta)) call WriteStatsSnap
                endif
            endif
        endif

        if (savemovie) then
            if (mod(time,tframe) < dt) then
                call Movie_xcut
                call Movie_ycut
                call Movie_zcut
                call Movie_outlet
            end if
        end if

        if( (ti(2) - tin(1)) .gt. walltimemax) errorcode = 334
        if( ntime .eq. ntst ) errorcode = 555 
        call MpiBcastInt(errorcode)

        !EP! Conditional exits
        if(errorcode.ne.0) then

            !EP! dt too small
            if(errorcode.eq.166) call QuitRoutine(tin,.false.,errorcode)

            !EP! cfl too high    
            if(errorcode.eq.165) call QuitRoutine(tin,.false.,errorcode)

            !EP! velocities diverged
            ! if(errorcode.eq.266) call QuitRoutine(tin,.false.,errorcode)
            if(errorcode.eq.266) call QuitRoutine(tin,.true.,errorcode)
            
            !EP! mass not conserved
            if(errorcode.eq.169) call QuitRoutine(tin,.false.,errorcode)

            !EP! Physical time exceeded tmax, no error; normal quit
            if(errorcode.eq.333) call QuitRoutine(tin,.true.,errorcode)

            !EP! walltime exceeded walltimemax, no error; normal quit
            if(errorcode.eq.334) call QuitRoutine(tin,.true.,errorcode)

            !RS! FFT input not correct
            if(errorcode.eq.444) call QuitRoutine(tin,.false.,errorcode)

            !RS! maximum number of timesteps reached, no error; normal quit
            if(errorcode.eq.555) call QuitRoutine(tin,.true.,errorcode)

            errorcode = 100 !EP already finalized

            exit

        endif

    enddo !EP main loop

    call QuitRoutine(tin,.true.,errorcode)

    end