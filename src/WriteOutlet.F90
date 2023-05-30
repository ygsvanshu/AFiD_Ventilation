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
    use movie_indices, only: mov_zcut
 	
    implicit none
	
	character*50    :: filename,dsetname

    if (xend(3).eq.nzm) then

        filename = trim('continua_outlet.h5')
        dsetname = trim('')
        
        call HdfCreatePath(dsetname,filename,comm_zcut)
    
        dsetname = trim("vx")
        mov_zcut(1:nx,xstart(2):xend(2)) = 0.5d0*(vx(1:nx,xstart(2):xend(2),nzm) + vx(1:nx,xstart(2):xend(2),nz))
        call HdfWriteReal2D_Z(dsetname,filename,mov_zcut)

        dsetname = trim("vy")
        mov_zcut(1:nx,xstart(2):xend(2)) = 0.5d0*(vy(1:nx,xstart(2):xend(2),nzm) + vy(1:nx,xstart(2):xend(2),nz))
        call HdfWriteReal2D_Z(dsetname,filename,mov_zcut)
    
        dsetname = trim("vz")
        mov_zcut(1:nx,xstart(2):xend(2)) = vz(1:nx,xstart(2):xend(2),nz)
        call HdfWriteReal2D_Z(dsetname,filename,mov_zcut)

        dsetname = trim("temp")
        mov_zcut(1:nx,xstart(2):xend(2)) = 0.5d0*(temp(1:nx,xstart(2):xend(2),nzm) + temp(1:nx,xstart(2):xend(2),nz))
        call HdfWriteReal2D_Z(dsetname,filename,mov_zcut)

        dsetname = trim("co2")
        mov_zcut(1:nx,xstart(2):xend(2)) = 0.5d0*(co2(1:nx,xstart(2):xend(2),nzm) + co2(1:nx,xstart(2):xend(2),nz))
        call HdfWriteReal2D_Z(dsetname,filename,mov_zcut)

        dsetname = trim("h2o")
        mov_zcut(1:nx,xstart(2):xend(2)) = 0.5d0*(h2o(1:nx,xstart(2):xend(2),nzm) + h2o(1:nx,xstart(2):xend(2),nz))
        call HdfWriteReal2D_Z(dsetname,filename,mov_zcut)
    
    end if
    
    return

end subroutine WriteOutlet
