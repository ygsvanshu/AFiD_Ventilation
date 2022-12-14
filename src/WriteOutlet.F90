!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                           ! 
!    FILE: WriteOutlet.F90                  !
!    CONTAINS: subroutine WriteOutlet       !
!                                           ! 
!    PURPOSE: Save the flow at the outlet   !
!    to continue the calculation            !
!                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine WriteOutlet
	
    use param
    use decomp_2d, only: xstart,xend,nrank
    use local_arrays, only: vz,vy,vx,temp,co2,h2o
 	
    implicit none
	
	character*50    :: filename,dsetname

    if (xend(3).eq.nzm) then

        filename = trim('continua_outlet.h5')
    
        dsetname = trim("vx")
        call HdfWriteReal2D_Z(dsetname,filename,0.5d0*(vx(1:nx,xstart(2):xend(2),nzm) + vx(1:nx,xstart(2):xend(2),nz)))

        dsetname = trim("vy")
        call HdfWriteReal2D_Z(dsetname,filename,0.5d0*(vy(1:nx,xstart(2):xend(2),nzm) + vy(1:nx,xstart(2):xend(2),nz)))
    
        dsetname = trim("vz")
        call HdfWriteReal2D_Z(dsetname,filename,vz(1:nx,xstart(2):xend(2),nz))

        dsetname = trim("temp")
        call HdfWriteReal2D_Z(dsetname,filename,0.5d0*(temp(1:nx,xstart(2):xend(2),nzm) + temp(1:nx,xstart(2):xend(2),nz)))

        dsetname = trim("co2")
        call HdfWriteReal2D_Z(dsetname,filename,0.5d0*(co2(1:nx,xstart(2):xend(2),nzm) + co2(1:nx,xstart(2):xend(2),nz)))

        dsetname = trim("h2o")
        call HdfWriteReal2D_Z(dsetname,filename,0.5d0*(h2o(1:nx,xstart(2):xend(2),nzm) + h2o(1:nx,xstart(2):xend(2),nz)))
    
    end if
    
    return

end subroutine WriteOutlet
