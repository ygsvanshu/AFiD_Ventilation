!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ReadInputFile.F90                              !
!    CONTAINS: subroutine ReadInputFile                   !
!                                                         ! 
!    PURPOSE: Read parameters from bou.in file            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ReadInputFile

    use param

    implicit none
    
    character(len=4) :: dummy
    integer flagstat,flagbal,stst3flag
    logical fexist

    ! ==============================================================================
    !             !NEW! New BOU.IN STYLE !NEW!
    ! ==============================================================================
    call read_from_bouin
    ! ==============================================================================
    ! ==============================================================================

    nx=nxm+1
    ny=nym+1                                                          
    nz=nzm+1                                                          

    ! DEFINITIONS FOR THE NATURAL CONVECTION
    ren = dsqrt(ray/pra)
    pec = dsqrt(pra*ray)
    pi  = 2.d0*dasin(1.d0)

    if ((statread).and.(.not. readflow)) write(6,*) 'Warning: Restarting flowfield with statistics read'
    savemovie = (savemovie_x.or.savemovie_y.or.savemovie_z.or.savemovie_i.or.savemovie_o)

    return 

end subroutine ReadInputFile

! =ModR01=Robert=2020-07-31=====================================================
!  - !NEW! New BOU.IN STYLE !NEW!
! Some subroutines required for new reading style of bou in parameters
! - Subroutine read_from_bouin    Reads the parameters from Input file bou.in
! - Subroutine Scan_string        Scans for valid strings
! - Subroutine STOP_CONFIG        Stop configuration if errors are detected
!
! =ModR09=Robert=2020-10-28=====================================================
!  - Change of input parameter:              TOUT -> NOUT
!  - New nomenclature for internal variable: NOUT -> NARG

SUBROUTINE read_from_bouin

    USE param
    
    IMPLICIT NONE  
    INTEGER i, narg
    INTEGER, PARAMETER :: LARGE_INTEGER=200, nimp=100
    CHARACTER*200 line
    CHARACTER*100 ss(nimp)
    CHARACTER*1   stringdummy1

    open(unit=15,file='bou_newStyle.in',status='old')
    DO 222 i=1,LARGE_INTEGER
        read(15,'(a200)',end=444) line

        IF(line(1:3)=='101') THEN 
        ! ###### READ DATA ######  
            call scan_string (line, 1, ss, narg)
            stringdummy1=ss(1)
            if('y'==stringdummy1) then
                readflow=.true.
            elseif('n'==stringdummy1) then
                readflow=.false.
            else
                Write(*,*)"ERROR: Input value of parameter NREAD not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'            "
                call stop_config
            endif

        ELSEIF(line(1:3)=='102') THEN 
        ! ###### READ STATS ######  
            call scan_string (line, 1, ss, narg)
            stringdummy1=ss(1)
            if('y'==stringdummy1) then
                statread = .true.
            elseif('n'==stringdummy1) then
                statread = .false.
            else
                Write(*,*)"ERROR: Input value of parameter STAREAD not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'             "
                call stop_config
            endif

        ELSEIF(line(1:3)=='201') THEN 
        ! ###### GRID DIMENSION ######  
            call scan_string (line, 3, ss, narg)
            read(ss(1),*) nxm
            read(ss(2),*) nym
            read(ss(3),*) nzm

        ELSEIF(line(1:3)=='202') THEN 
        ! ###### AXIS EXTENDS ######  
            call scan_string (line, 3, ss, narg)
            read(ss(1),*) alx3
            read(ss(2),*) ylen
            read(ss(3),*) zlen

        ELSEIF(line(1:3)=='203') THEN 
        ! ###### AXIS PARAMETER 1 ######  
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) istr3

        ELSEIF(line(1:3)=='204') THEN 
        ! ###### AXIS PARAMETER 2 ######  
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) str3

        ELSEIF(line(1:3)=='301') THEN 
        ! ###### MAX STEPS ######  
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) ntst

        ELSEIF(line(1:3)=='302') THEN 
        ! ###### MAX TIME ######  
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) tmax

        ELSEIF(line(1:3)=='303') THEN 
        ! ###### MAX WALLCLOCK TIME ######  
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) walltimemax

        ELSEIF(line(1:3)=='401') THEN 
        ! ###### NUMERICAL SCHEME ######  
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) nsst

        ELSEIF(line(1:3)=='402') THEN 
        ! ###### VARY TIME STEP ######  
            call scan_string (line, 1, ss, narg)
            stringdummy1=ss(1)
            if('y'==stringdummy1) then
                variabletstep = .true.
            elseif('n'==stringdummy1) then
                variabletstep = .false.
            else
                Write(*,*)"ERROR: Input value of parameter IDTV not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'             "
                call stop_config
            endif

        ELSEIF(line(1:3)=='403') THEN 
        ! ###### FIRST TIME STEP ######  
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) dt

        ELSEIF(line(1:3)=='404') THEN 
        ! ###### MIN MAX STEP SIZE ######  
            call scan_string (line, 2, ss, narg)
            read(ss(1),*) dtmin
            read(ss(2),*) dtmax

        ELSEIF(line(1:3)=='405') THEN 
        ! ###### TIMESTEP STABILITY ######  
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) limitCFL

        ELSEIF(line(1:3)=='406') THEN 
        ! ###### RESID ######  
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) resid

        ELSEIF(line(1:3)=='407') THEN 
        ! ###### MAX VELOCITY ######  
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) limitVel

        ELSEIF(line(1:3)=='501') THEN 
        ! ###### PRANDTL NUMBER ######  
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) pra

        ELSEIF(line(1:3)=='502') THEN 
        ! ###### RAYLEIGH NUMBER ######  
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) ray

        ELSEIF(line(1:3)=='503') THEN 
        ! ###### THERMAL EXPANSION RATIO FOR CO2 ######  
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) lambda_co2

        ELSEIF(line(1:3)=='504') THEN 
        ! ###### THERMAL EXPANSION RATIO FOR H2O VAPOUR ######  
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) lambda_h2o

        ELSEIF(line(1:3)=='505') THEN 
        ! ###### INLET VENT DIMENSION AND LOCATION ######  
            call scan_string (line, 2, ss, narg)
            read(ss(1),*) ilen
            read(ss(2),*) iheight

        ELSEIF(line(1:3)=='506') THEN 
        ! ###### OUTLET VENT DIMENSION AND LOCATION ######  
            call scan_string (line, 2, ss, narg)
            read(ss(1),*) olen
            read(ss(2),*) oheight

        ELSEIF(line(1:3)=='507') THEN 
        ! ###### INLET FLOW VELOCITY ######  
            call scan_string (line, 2, ss, narg)
            read(ss(1),*) ivel
            read(ss(2),*) tvel

        ELSEIF(line(1:3)=='508') THEN 
            ! ###### OUTLET PARAMETERS ######  
            call scan_string (line, 3, ss, narg)
            read(ss(1),*) ocou
            read(ss(2),*) ovsc
            read(ss(3),*) odst

        ELSEIF(line(1:3)=='509') THEN 
        ! ###### PERSON AND BREATH IBM ENABLED ######  
            call scan_string (line, 2, ss, narg)
            stringdummy1=ss(1)
            if('y'==stringdummy1) then
                person_on = .true.
            elseif('n'==stringdummy1) then
                person_on = .false.
            else
                Write(*,*)"ERROR: Input value of parameter PERSON not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'             "
                call stop_config
            endif
            stringdummy1=ss(2)
            if('y'==stringdummy1) then
                breath_on = .true.
            elseif('n'==stringdummy1) then
                breath_on = .false.
            else
                Write(*,*)"ERROR: Input value of parameter BREATH not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'             "
                call stop_config
            endif

        ELSEIF(line(1:3)=='510') THEN 
        ! ###### FILENAME OF PERSON OBJ GEOMETRY ######  
            call scan_string (line, 1, ss, narg)
            person_objfile = ss(1)

        ELSEIF(line(1:3)=='511') THEN 
        ! ###### PERSON GEOMETRY LOCATION AND SCALE ######  
            call scan_string (line, 4, ss, narg)
            read(ss(1),*) personx
            read(ss(2),*) persony
            read(ss(3),*) personz
            read(ss(4),*) sclf

        ELSEIF(line(1:3)=='512') THEN 
        ! ###### BREATH GEOMETRY LOCATION ######  
            call scan_string (line, 3, ss, narg)
            read(ss(1),*) breathx
            read(ss(2),*) breathy
            read(ss(3),*) breathz

        ELSEIF(line(1:3)=='513') THEN 
        ! ###### BREATHING KERNEL OPTIONS ######  
            call scan_string (line, 2, ss, narg)
            read(ss(1),*) kernel_width_space
            read(ss(2),*) kernel_width_time

        ELSEIF(line(1:3)=='514') THEN 
        ! ###### BOUNDARY CONDITIONS ######  
            call scan_string (line, 2, ss, narg)
            read(ss(1),*) inslws
            read(ss(2),*) inslwn

        ELSEIF(line(1:3)=='601') THEN 
        ! ###### STATISTICS ######  
            call scan_string (line, 3, ss, narg)
            stringdummy1=ss(1)
            if('y'==stringdummy1) then
                statcalc = .true.
            elseif('n'==stringdummy1) then
                statcalc = .false.
            else
                Write(*,*)"ERROR: Input value of parameter STATON not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'             "
                call stop_config
            endif
            read(ss(2),*) tsta
            read(ss(3),*) nout
            
        ELSEIF(line(1:3)=='602') THEN 
        ! ###### 2D STATISTICS SLICES ######
            call scan_string (line, 3, ss, narg)
            read(ss(1),*) stats2Dx
            read(ss(2),*) stats2Dy
            read(ss(3),*) stats2Dz

        ELSEIF(line(1:3)=='603') THEN 
        ! ###### SNAPSHOTS ######
            call scan_string (line, 3, ss, narg)
            stringdummy1=ss(1)
            if('y'==stringdummy1) then
                savesnap = .true.
            elseif('n'==stringdummy1) then
                savesnap = .false.
            else
                Write(*,*)"ERROR: Input value of parameter SNAPON not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'             "
                call stop_config
            endif
            read(ss(2),*) startsnap
            read(ss(3),*) tsnap

        ELSEIF(line(1:3)=='604') THEN 
        ! ###### MOVIE OUTPUT ######  
            call scan_string (line, 6, ss, narg)

            stringdummy1=ss(1)
            if('y'==stringdummy1) then
                savemovie_x=.true.
            elseif('n'==stringdummy1) then
                savemovie_x=.false.
            else
                Write(*,*)"ERROR: Input value of parameter MX not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'             "
                call stop_config
            endif

            stringdummy1=ss(2)
            if('y'==stringdummy1) then
                savemovie_y=.true.
            elseif('n'==stringdummy1) then
                savemovie_y=.false.
            else
                Write(*,*)"ERROR: Input value of parameter MY not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'             "
                call stop_config
            endif

            stringdummy1=ss(3)
            if('y'==stringdummy1) then
                savemovie_z=.true.
            elseif('n'==stringdummy1) then
                savemovie_z=.false.
            else
                Write(*,*)"ERROR: Input value of parameter MZ not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'             "
                call stop_config
            endif

            stringdummy1=ss(4)
            if('y'==stringdummy1) then
                savemovie_i=.true.
            elseif('n'==stringdummy1) then
                savemovie_i=.false.
            else
                Write(*,*)"ERROR: Input value of parameter MI not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'             "
                call stop_config
            endif

            stringdummy1=ss(5)
            if('y'==stringdummy1) then
                savemovie_o=.true.
            elseif('n'==stringdummy1) then
                savemovie_o=.false.
            else
                Write(*,*)"ERROR: Input value of parameter MO not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'             "
                call stop_config
            endif

            read(ss(6),*) tframe

        ELSEIF(line(1:3)=='605') THEN 
            ! ###### 2D MOVIE SLICES ######
                call scan_string (line, 3, ss, narg)
                read(ss(1),*) movie2Dx
                read(ss(2),*) movie2Dy
                read(ss(3),*) movie2Dz

        ELSEIF(line(1:3)=='606') THEN 
            ! ###### SAVE CONTINUA CHECKPOINT ######
                call scan_string (line, 1, ss, narg)
                read(ss(1),*) tcontinua

        ENDIF 

222 CONTINUE

444 close(15)

END SUBROUTINE read_from_bouin

SUBROUTINE SCAN_STRING (S, N, SS, M)

    IMPLICIT NONE 

    INTEGER I, J, N, M, APO(2*N) 
    CHARACTER*200 S ; CHARACTER*100 SS(N)
    CHARACTER*1, PARAMETER :: PRIME="'"

    !
    ! Given the input string s, this routine determines the sub-strings ss of s 
    ! which are delimited by consecutive pairs of primes ('). n (input) is the 
    ! expected number of substrings, and m (output) is the effective number 
    ! found. Apo is a pointer to the primes. The execution is stopped whenever 
    ! i) n=0, ii) m=0, iii) m/=n, or iv) if m is odd. The maximum lenght of the 
    ! substrings is 20, that of the input string is 200.  == GS Nov 17 2007 == 
    !
    ! ----- Exits for n=0
    if(n==0) then 
        Write(*, *) "SCAN_STRING has nothing to do" 
        call stop_config 
    endif
    !
    ! ----- Looking for primes (', not ") within the input string  
    i=0
    do j=1, 200
        if(s(j:j)==prime) then 
            i=i+1 ; apo(i)=j
        endif
    enddo 
    m=i/2 
    !
    if(i==0) then 
    ! ----- Exits if no primes are found 
        Write(*, *) "SCAN_STRING has found NO primes" 
        Write(*,'(a100)') s 
        call stop_config 
    !
    elseif(mod(i,2)/=0) then 
    ! ----- Exits if an odd number of primes is found
        Write(*, *) "SCAN_STRING has found an even number of primes:", i
        call stop_config 
    !
    elseif(m/=n) then 
    ! ----- Exits if m/=n, otherwise determines the substrings
        Write(*, *) "SCAN_STRING has found ", m, " substrings"
        Write(*, *) "SCAN_STRING  expected ", n, " substrings"
        call stop_config 
    else
        do i=1, n 
            ss(i)=s(apo(2*i-1)+1:apo(2*i)-1)
        enddo
    endif	

END SUBROUTINE scan_string 

SUBROUTINE STOP_CONFIG

    use mpih

    IMPLICIT NONE
    write(* ,*) 'STOP_CONFIG: The program will STOP! -----------------------------' 
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    STOP 1

END SUBROUTINE STOP_CONFIG 	
