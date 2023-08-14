!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                           ! 
!    FILE: MovieRoutines.F90                                !
!    CONTAINS: subroutines Movie_xcut, Movie_ycut           !
!              and Movie_zcut                               !
!                                                           ! 
!    PURPOSE: Write down xcut, ycut and zcut snapshots      ! 
!             for flow variables                            !
!                                                           !
!    CREATED BY: G.S.Yerragolam (g.s.yerragolam@utwente.nl) !
!                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InitMovie

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: ibm_body
    use movie_indices
    use mpih

    implicit none

    integer         :: i,j,k
    real            :: cpc,cmc
    character*50    :: dsetname,filename,frame

    ! For the X grid
    do k=1,nx-1
        if ((xc(k).le.movie2Dx).and.(xc(k+1).gt.movie2Dx)) mov_xk = k
    end do
    do j=1,nym-1
        if ((ym(j).le.movie2Dy).and.(ym(j+1).gt.movie2Dy)) mov_xj = j
    end do
    do i=1,nzm-1
        if ((zm(i).le.movie2Dz).and.(zm(i+1).gt.movie2Dz)) mov_xi = i
    end do

    ! For the Y grid
    do k=1,nxm-1
        if ((xm(k).le.movie2Dx).and.(xm(k+1).gt.movie2Dx)) mov_yk = k
    end do
    do j=1,ny-1
        if ((yc(j).le.movie2Dy).and.(yc(j+1).gt.movie2Dy)) mov_yj = j
    end do
    do i=1,nzm-1
        if ((zm(i).le.movie2Dz).and.(zm(i+1).gt.movie2Dz)) mov_yi = i
    end do

    ! For the Z grid
    do k=1,nxm-1
        if ((xm(k).le.movie2Dx).and.(xm(k+1).gt.movie2Dx)) mov_zk = k
    end do
    do j=1,nym-1
        if ((ym(j).le.movie2Dy).and.(ym(j+1).gt.movie2Dy)) mov_zj = j
    end do
    do i=1,nz-1
        if ((zc(i).le.movie2Dz).and.(zc(i+1).gt.movie2Dz)) mov_zi = i
    end do

    ! For the C grid
    do k=1,nxm-1
        if ((xm(k).le.movie2Dx).and.(xm(k+1).gt.movie2Dx)) mov_ck = k
    end do
    do j=1,nym-1
        if ((ym(j).le.movie2Dy).and.(ym(j+1).gt.movie2Dy)) mov_cj = j
    end do
    do i=1,nzm-1
        if ((zm(i).le.movie2Dz).and.(zm(i+1).gt.movie2Dz)) mov_ci = i
    end do

    filename = trim('Results/movie_xcut.h5')

    dsetname = trim("vx/")
    call HdfCreatePath(dsetname,filename,comm_xcut)
    dsetname = trim("vy/")
    call HdfCreatePath(dsetname,filename,comm_xcut)
    dsetname = trim("vz/")
    call HdfCreatePath(dsetname,filename,comm_xcut)
    dsetname = trim("pr/")
    call HdfCreatePath(dsetname,filename,comm_xcut)
    dsetname = trim("temp/")
    call HdfCreatePath(dsetname,filename,comm_xcut)
    dsetname = trim("co2/")
    call HdfCreatePath(dsetname,filename,comm_xcut)
    dsetname = trim("h2o/")
    call HdfCreatePath(dsetname,filename,comm_xcut)
    
    cpc = (movie2Dx     - xm(mov_ck))/(xm(mov_ck+1) - xm(mov_ck))
    cmc = (xm(mov_ck+1) - movie2Dx)  /(xm(mov_ck+1) - xm(mov_ck))

    dsetname = trim("ibm_body")
    do j=xstart(2),xend(2)
        do i=xstart(3),xend(3)
            mov_xcut(j,i) = cmc*ibm_body(mov_ck,j,i) + cpc*ibm_body(mov_ck+1,j,i)
        end do
    end do
    call HdfWriteReal2D_X(dsetname,filename,mov_xcut)

    if ((yc(xstart(2)).le.movie2Dy).and.(yc(xend(2)+lvlhalo).gt.movie2Dy)) then
        filename = trim('Results/movie_ycut.h5')

        dsetname = trim("vx/")
        call HdfCreatePath(dsetname,filename,comm_ycut)
        dsetname = trim("vy/")
        call HdfCreatePath(dsetname,filename,comm_ycut)
        dsetname = trim("vz/")
        call HdfCreatePath(dsetname,filename,comm_ycut)
        dsetname = trim("pr/")
        call HdfCreatePath(dsetname,filename,comm_ycut)
        dsetname = trim("temp/")
        call HdfCreatePath(dsetname,filename,comm_ycut)
        dsetname = trim("co2/")
        call HdfCreatePath(dsetname,filename,comm_ycut)
        dsetname = trim("h2o/")
        call HdfCreatePath(dsetname,filename,comm_ycut)

        cpc = (movie2Dy     - ym(mov_cj))/(ym(mov_cj+1) - ym(mov_cj))
        cmc = (ym(mov_cj+1) - movie2Dy)  /(ym(mov_cj+1) - ym(mov_cj))

        dsetname = trim("ibm_body")
        do k=1,nx
            do i=xstart(3),xend(3)
                mov_ycut(k,i) = cmc*ibm_body(k,mov_cj,i) + cpc*ibm_body(k,mov_cj+1,i)
            end do
        end do
        call HdfWriteReal2D_Y(dsetname,filename,mov_ycut)
    end if

    if ((zc(xstart(3)).le.movie2Dz).and.(zc(xend(3)+lvlhalo).gt.movie2Dz)) then
        filename = trim('Results/movie_zcut.h5')

        dsetname = trim("vx/")
        call HdfCreatePath(dsetname,filename,comm_zcut)
        dsetname = trim("vy/")
        call HdfCreatePath(dsetname,filename,comm_zcut)
        dsetname = trim("vz/")
        call HdfCreatePath(dsetname,filename,comm_zcut)
        dsetname = trim("pr/")
        call HdfCreatePath(dsetname,filename,comm_zcut)
        dsetname = trim("temp/")
        call HdfCreatePath(dsetname,filename,comm_zcut)
        dsetname = trim("co2/")
        call HdfCreatePath(dsetname,filename,comm_zcut)
        dsetname = trim("h2o/")
        call HdfCreatePath(dsetname,filename,comm_zcut)

        cpc = (movie2Dz     - zm(mov_ci))/(zm(mov_ci+1) - zm(mov_ci))
        cmc = (zm(mov_ci+1) - movie2Dz)  /(zm(mov_ci+1) - zm(mov_ci))

        dsetname = trim("ibm_body")
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_zcut(k,j) = cmc*ibm_body(k,j,mov_ci) + cpc*ibm_body(k,j,mov_ci+1)
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_zcut)
    end if

    if (xend(3).eq.nzm) then
        filename = trim('Results/movie_outlet.h5')

        dsetname = trim("vx/")
        call HdfCreatePath(dsetname,filename,comm_zcut)
        dsetname = trim("vy/")
        call HdfCreatePath(dsetname,filename,comm_zcut)
        dsetname = trim("vz/")
        call HdfCreatePath(dsetname,filename,comm_zcut)
        dsetname = trim("pr/")
        call HdfCreatePath(dsetname,filename,comm_zcut)
        dsetname = trim("temp/")
        call HdfCreatePath(dsetname,filename,comm_zcut)
        dsetname = trim("co2/")
        call HdfCreatePath(dsetname,filename,comm_zcut)
        dsetname = trim("h2o/")
        call HdfCreatePath(dsetname,filename,comm_zcut)
    end if

end subroutine InitMovie

subroutine Movie_xcut

    use param
    use local_arrays, only: vx,vy,vz,pr,temp,co2,h2o
    use decomp_2d, only: xstart,xend
    use movie_indices
    use mpih

    implicit none

    integer             :: k,j,i
	real                :: cpx,cmx,cpy,cmy,cpz,cmz,cpc,cmc
    character*50        :: dsetname,filename,frame

    filename = trim('Results/movie_xcut.h5')
    write (frame,"(i5.5)") nint(time/tframe)

    cpx = (movie2Dx     - xc(mov_xk))/(xc(mov_xk+1) - xc(mov_xk))
    cmx = (xc(mov_xk+1) - movie2Dx)  /(xc(mov_xk+1) - xc(mov_xk))
    cpy = (movie2Dx     - xm(mov_yk))/(xm(mov_yk+1) - xm(mov_yk))
    cmy = (xm(mov_yk+1) - movie2Dx)  /(xm(mov_yk+1) - xm(mov_yk))
    cpz = (movie2Dx     - xm(mov_zk))/(xm(mov_zk+1) - xm(mov_zk))
    cmz = (xm(mov_zk+1) - movie2Dx)  /(xm(mov_zk+1) - xm(mov_zk))
    cpc = (movie2Dx     - xm(mov_ck))/(xm(mov_ck+1) - xm(mov_ck))
    cmc = (xm(mov_ck+1) - movie2Dx)  /(xm(mov_ck+1) - xm(mov_ck))

    dsetname = trim("vx")//'/'//trim(frame)
    do j=xstart(2),xend(2)
        do i=xstart(3),xend(3)
            mov_xcut(j,i) = cmx*vx(mov_xk,j,i) + cpx*vx(mov_xk+1,j,i)
        end do
    end do
    call HdfWriteReal2D_X(dsetname,filename,mov_xcut)

    dsetname = trim("vy")//'/'//trim(frame)
    do j=xstart(2),xend(2)
        do i=xstart(3),xend(3)
            mov_xcut(j,i) = cmy*vy(mov_yk,j,i) + cpy*vy(mov_yk+1,j,i)
        end do
    end do
    call HdfWriteReal2D_X(dsetname,filename,mov_xcut)

    dsetname = trim("vz")//'/'//trim(frame)
    do j=xstart(2),xend(2)
        do i=xstart(3),xend(3)
            mov_xcut(j,i) = cmz*vz(mov_zk,j,i) + cpz*vz(mov_zk+1,j,i)
        end do
    end do
    call HdfWriteReal2D_X(dsetname,filename,mov_xcut)

    dsetname = trim("pr")//'/'//trim(frame)
    do j=xstart(2),xend(2)
        do i=xstart(3),xend(3)
            mov_xcut(j,i) = cmc*pr(mov_ck,j,i) + cpc*pr(mov_ck+1,j,i)
        end do
    end do
    call HdfWriteReal2D_X(dsetname,filename,mov_xcut)

    dsetname = trim("temp")//'/'//trim(frame)
    do j=xstart(2),xend(2)
        do i=xstart(3),xend(3)
            mov_xcut(j,i) = cmx*temp(mov_xk,j,i) + cpx*temp(mov_xk+1,j,i)
        end do
    end do
    call HdfWriteReal2D_X(dsetname,filename,mov_xcut)

    dsetname = trim("co2")//'/'//trim(frame)
    do j=xstart(2),xend(2)
        do i=xstart(3),xend(3)
            mov_xcut(j,i) = cmx*co2(mov_xk,j,i) + cpx*co2(mov_xk+1,j,i)
        end do
    end do
    call HdfWriteReal2D_X(dsetname,filename,mov_xcut)

    dsetname = trim("h2o")//'/'//trim(frame)
    do j=xstart(2),xend(2)
        do i=xstart(3),xend(3)
            mov_xcut(j,i) = cmx*h2o(mov_xk,j,i) + cpx*h2o(mov_xk+1,j,i)
        end do
    end do
    call HdfWriteReal2D_X(dsetname,filename,mov_xcut)

end subroutine Movie_xcut

subroutine Movie_ycut

    use param
    use local_arrays, only: vx,vy,vz,pr,temp,co2,h2o
    use decomp_2d, only: xstart,xend
    use movie_indices
    use mpih

    implicit none

    integer             :: k,j,i
	real                :: cpx,cmx,cpy,cmy,cpz,cmz,cmc,cpc
    character*50        :: dsetname,filename,frame

    if ((yc(xstart(2)).le.movie2Dy).and.(yc(xend(2)+lvlhalo).gt.movie2Dy)) then

        filename = trim('Results/movie_ycut.h5')
        write (frame,"(i5.5)") nint(time/tframe)

        cpx = (movie2Dy     - ym(mov_xj))/(ym(mov_xj+1) - ym(mov_xj))
        cmx = (ym(mov_xj+1) - movie2Dy)  /(ym(mov_xj+1) - ym(mov_xj))
        cpy = (movie2Dy     - yc(mov_yj))/(yc(mov_yj+1) - yc(mov_yj))
        cmy = (yc(mov_yj+1) - movie2Dy)  /(yc(mov_yj+1) - yc(mov_yj))
        cpz = (movie2Dy     - ym(mov_zj))/(ym(mov_zj+1) - ym(mov_zj))
        cmz = (ym(mov_zj+1) - movie2Dy)  /(ym(mov_zj+1) - ym(mov_zj))
        cpc = (movie2Dy     - ym(mov_cj))/(ym(mov_cj+1) - ym(mov_cj))
        cmc = (ym(mov_cj+1) - movie2Dy)  /(ym(mov_cj+1) - ym(mov_cj))

        dsetname = trim("vx")//'/'//trim(frame)
        do k=1,nx
            do i=xstart(3),xend(3)
                mov_ycut(k,i) = cmx*vx(k,mov_xj,i) + cpx*vx(k,mov_xj+1,i)
            end do
        end do
        call HdfWriteReal2D_Y(dsetname,filename,mov_ycut)

        dsetname = trim("vy")//'/'//trim(frame)
        do k=1,nx
            do i=xstart(3),xend(3)
                mov_ycut(k,i) = cmy*vy(k,mov_yj,i) + cpy*vy(k,mov_yj+1,i)
            end do
        end do
        call HdfWriteReal2D_Y(dsetname,filename,mov_ycut)

        dsetname = trim("vz")//'/'//trim(frame)
        do k=1,nx
            do i=xstart(3),xend(3)
                mov_ycut(k,i) = cmz*vz(k,mov_zj,i) + cpz*vz(k,mov_zj+1,i)
            end do
        end do
        call HdfWriteReal2D_Y(dsetname,filename,mov_ycut)

        dsetname = trim("pr")//'/'//trim(frame)
        do k=1,nx
            do i=xstart(3),xend(3)
                mov_ycut(k,i) = cmc*pr(k,mov_cj,i) + cpc*pr(k,mov_cj+1,i)
            end do
        end do
        call HdfWriteReal2D_Y(dsetname,filename,mov_ycut)

        dsetname = trim("temp")//'/'//trim(frame)
        do k=1,nx
            do i=xstart(3),xend(3)
                mov_ycut(k,i) = cmx*temp(k,mov_xj,i) + cpx*temp(k,mov_xj+1,i)
            end do
        end do
        call HdfWriteReal2D_Y(dsetname,filename,mov_ycut)

        dsetname = trim("co2")//'/'//trim(frame)
        do k=1,nx
            do i=xstart(3),xend(3)
                mov_ycut(k,i) = cmx*co2(k,mov_xj,i) + cpx*co2(k,mov_xj+1,i)
            end do
        end do
        call HdfWriteReal2D_Y(dsetname,filename,mov_ycut)

        dsetname = trim("h2o")//'/'//trim(frame)
        do k=1,nx
            do i=xstart(3),xend(3)
                mov_ycut(k,i) = cmx*h2o(k,mov_xj,i) + cpx*h2o(k,mov_xj+1,i)
            end do
        end do
        call HdfWriteReal2D_Y(dsetname,filename,mov_ycut)

    end if

end subroutine Movie_ycut

subroutine Movie_zcut

    use param
    use local_arrays, only: vx,vy,vz,pr,temp,co2,h2o
    use decomp_2d, only: xstart,xend
    use movie_indices
    use mpih

    implicit none

    integer             :: k,j,i
	real                :: cpx,cmx,cpy,cmy,cpz,cmz,cmc,cpc
    character*50        :: filename,dsetname,frame

    if ((zc(xstart(3)).le.movie2Dz).and.(zc(xend(3)+lvlhalo).gt.movie2Dz)) then

        filename = trim('Results/movie_zcut.h5')
        write (frame,"(i5.5)") nint(time/tframe)

        cpx = (movie2Dz     - zm(mov_xi))/(zm(mov_xi+1) - zm(mov_xi))
        cmx = (zm(mov_xi+1) - movie2Dz)  /(zm(mov_xi+1) - zm(mov_xi))
        cpy = (movie2Dz     - zm(mov_yi))/(zm(mov_yi+1) - zm(mov_yi))
        cmy = (zm(mov_yi+1) - movie2Dz)  /(zm(mov_yi+1) - zm(mov_yi))
        cpz = (movie2Dz     - zc(mov_zi))/(zc(mov_zi+1) - zc(mov_zi))
        cmz = (zc(mov_zi+1) - movie2Dz)  /(zc(mov_zi+1) - zc(mov_zi))
        cpc = (movie2Dz     - zm(mov_ci))/(zm(mov_ci+1) - zm(mov_ci))
        cmc = (zm(mov_ci+1) - movie2Dz)  /(zm(mov_ci+1) - zm(mov_ci))

        dsetname = trim("vx")//'/'//trim(frame)
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_zcut(k,j) = cmx*vx(k,j,mov_xi) + cpx*vx(k,j,mov_xi+1)
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_zcut)

        dsetname = trim("vy")//'/'//trim(frame)
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_zcut(k,j) = cmy*vy(k,j,mov_yi) + cpy*vy(k,j,mov_yi+1)
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_zcut)
    
        dsetname = trim("vz")//'/'//trim(frame)
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_zcut(k,j) = cmz*vz(k,j,mov_zi) + cpz*vz(k,j,mov_zi+1)
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_zcut)

        dsetname = trim("pr")//'/'//trim(frame)
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_zcut(k,j) = cmc*pr(k,j,mov_ci) + cpc*pr(k,j,mov_ci+1)
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_zcut)

        dsetname = trim("temp")//'/'//trim(frame)
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_zcut(k,j) = cmx*temp(k,j,mov_xi) + cpx*temp(k,j,mov_xi+1)
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_zcut)

        dsetname = trim("co2")//'/'//trim(frame)
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_zcut(k,j) = cmx*co2(k,j,mov_xi) + cpx*co2(k,j,mov_xi+1)
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_zcut)

        dsetname = trim("h2o")//'/'//trim(frame)
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_zcut(k,j) = cmx*h2o(k,j,mov_xi) + cpx*h2o(k,j,mov_xi+1)
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_zcut)

    end if

end subroutine Movie_zcut

subroutine Movie_outlet

    use param
    use local_arrays, only: vx,vy,vz,pr,temp,co2,h2o
    use decomp_2d, only: xstart,xend
    use movie_indices
    use mpih

    implicit none

    integer             :: k,j,i
    character*50        :: dsetname,filename,frame

    if (xend(3).eq.nzm) then

        filename = trim('Results/movie_outlet.h5')
        write (frame,"(i5.5)") nint(time/tframe)
    
        dsetname = trim("vx")//'/'//trim(frame)
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_ocut(k,j) = 0.5d0*(vx(k,j,nzm) + vx(k,j,nz))
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_ocut)

        dsetname = trim("vy")//'/'//trim(frame)
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_ocut(k,j) = 0.5d0*(vy(k,j,nzm) + vy(k,j,nz))
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_ocut)
    
        dsetname = trim("vz")//'/'//trim(frame)
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_ocut(k,j) = vz(k,j,nz)
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_ocut)

        dsetname = trim("pr")//'/'//trim(frame)
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_ocut(k,j) = 0.5d0*(pr(k,j,nzm) + pr(k,j,nz))
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_ocut)

        dsetname = trim("temp")//'/'//trim(frame)
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_ocut(k,j) = 0.5d0*(temp(k,j,nzm) + temp(k,j,nz))
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_ocut)

        dsetname = trim("co2")//'/'//trim(frame)
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_ocut(k,j) = 0.5d0*(co2(k,j,nzm) + co2(k,j,nz))
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_ocut)

        dsetname = trim("h2o")//'/'//trim(frame)
        do k=1,nx
            do j=xstart(2),xend(2)
                mov_ocut(k,j) = 0.5d0*(h2o(k,j,nzm) + h2o(k,j,nz))
            end do
        end do
        call HdfWriteReal2D_Z(dsetname,filename,mov_ocut)

    end if

end subroutine Movie_outlet
