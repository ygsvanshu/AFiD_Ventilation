!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                           ! 
!    FILE: ReadOutlet.F90                                   !
!    CONTAINS: subroutine ReadOutlet                        !
!                                                           ! 
!    PURPOSE: Initialization routine. Reads in flow at the  !
!    outlet to continue the calculation                     !
!                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ReadOutlet

    use mpih
    use decomp_2d
    use local_arrays, only: vx,vy,vz,pr,temp,co2,h2o
    use param

    implicit none

    character*50                            :: filename,dsetname
    logical                                 :: fexist
    real,dimension(1:nx,xstart(2):xend(2))  :: readwrite_outlet
      
    filename = trim('continua_outlet.h5')
    inquire(file=filename,exist=fexist)
    if (fexist) then 
        if (ismaster) write (6,*) 'Reading outlet'
        if (xend(3).eq.nzm) then 

            dsetname = trim('vx')
            call HdfReadReal2D_Z(filename,dsetname,readwrite_outlet,xstart(2),xend(2),nx,ny)
            vx(1:nx,xstart(2):xend(2),nz) = &
            2.0d0*readwrite_outlet(1:nx,xstart(2):xend(2)) - vx(1:nx,xstart(2):xend(2),nzm)

            dsetname = trim('vy')
            call HdfReadReal2D_Z(filename,dsetname,readwrite_outlet,xstart(2),xend(2),nx,ny)
            vy(1:nx,xstart(2):xend(2),nz) = &
            2.0d0*readwrite_outlet(1:nx,xstart(2):xend(2)) - vy(1:nx,xstart(2):xend(2),nzm)

            dsetname = trim('vz')
            call HdfReadReal2D_Z(filename,dsetname,readwrite_outlet,xstart(2),xend(2),nx,ny)
            vz(1:nx,xstart(2):xend(2),nz) = readwrite_outlet(1:nx,xstart(2):xend(2))

            dsetname = trim('pr')
            call HdfReadReal2D_Z(filename,dsetname,readwrite_outlet,xstart(2),xend(2),nx,ny)
            pr(1:nx,xstart(2):xend(2),nz) = &
            2.0d0*readwrite_outlet(1:nx,xstart(2):xend(2)) - pr(1:nx,xstart(2):xend(2),nzm)

            dsetname = trim('temp')
            call HdfReadReal2D_Z(filename,dsetname,readwrite_outlet,xstart(2),xend(2),nx,ny)
            temp(1:nx,xstart(2):xend(2),nz) = &
            2.0d0*readwrite_outlet(1:nx,xstart(2):xend(2)) - temp(1:nx,xstart(2):xend(2),nzm)

            dsetname = trim('co2')
            call HdfReadReal2D_Z(filename,dsetname,readwrite_outlet,xstart(2),xend(2),nx,ny)
            co2(1:nx,xstart(2):xend(2),nz) = &
            2.0d0*readwrite_outlet(1:nx,xstart(2):xend(2)) - co2(1:nx,xstart(2):xend(2),nzm)
            
            dsetname = trim('h2o')
            call HdfReadReal2D_Z(filename,dsetname,readwrite_outlet,xstart(2),xend(2),nx,ny)
            h2o(1:nx,xstart(2):xend(2),nz) = &
            2.0d0*readwrite_outlet(1:nx,xstart(2):xend(2)) - h2o(1:nx,xstart(2):xend(2),nzm)
        
        end if
    end if

    return

end subroutine ReadOutlet