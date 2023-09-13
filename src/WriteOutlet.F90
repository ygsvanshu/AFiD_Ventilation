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
    use local_arrays, only: vz,vy,vx,temp,pr,co2,h2o
    use movie_indices, only: mov_ocut
 	
    implicit none
	
    integer         :: k,j,i
	character*50    :: filename,dsetname

    if (xend(3).eq.nzm) then

        filename = trim('continua_outlet.h5')
    
        call HdfParallelCreateBlankFile(filename,comm_zcut)

        dsetname = trim("/vx")
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_ocut(k,j) = 0.5d0*(vx(k,j,nzm) + vx(k,j,nz))
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_ocut)

        dsetname = trim("/vy")
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_ocut(k,j) = 0.5d0*(vy(k,j,nzm) + vy(k,j,nz))
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_ocut)
    
        dsetname = trim("/vz")
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_ocut(k,j) = vz(k,j,nz)
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_ocut)

        dsetname = trim("/pr")
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_ocut(k,j) = 0.5d0*(pr(k,j,nzm) + pr(k,j,nz))
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_ocut)

        dsetname = trim("/temp")
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_ocut(k,j) = 0.5d0*(temp(k,j,nzm) + temp(k,j,nz))
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_ocut)

        dsetname = trim("/co2")
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_ocut(k,j) = 0.5d0*(co2(k,j,nzm) + co2(k,j,nz))
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_ocut)

        dsetname = trim("/h2o")
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_ocut(k,j) = 0.5d0*(h2o(k,j,nzm) + h2o(k,j,nz))
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_ocut)

    end if
    
    return

end subroutine WriteOutlet
