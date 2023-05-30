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
    mov_xcut(xstart(2):xend(2),xstart(3):xend(3)) = cmc*ibm_body(mov_ck,xstart(2):xend(2),xstart(3):xend(3)) + cpc*ibm_body(mov_ck+1,xstart(2):xend(2),xstart(3):xend(3))
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
        mov_ycut(1:nx,xstart(3):xend(3)) = cmc*ibm_body(1:nx,mov_cj,xstart(3):xend(3)) + cpc*ibm_body(1:nx,mov_cj+1,xstart(3):xend(3))
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
        mov_zcut(1:nx,xstart(2):xend(2)) = cmc*ibm_body(1:nx,xstart(2):xend(2),mov_ci) + cpc*ibm_body(1:nx,xstart(2):xend(2),mov_ci+1)
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
    mov_xcut(xstart(2):xend(2),xstart(3):xend(3)) = cmx*vx(mov_xk,xstart(2):xend(2),xstart(3):xend(3)) + cpx*vx(mov_xk+1,xstart(2):xend(2),xstart(3):xend(3))
    call HdfWriteReal2D_X(dsetname,filename,mov_xcut)

    dsetname = trim("vy")//'/'//trim(frame)
    mov_xcut(xstart(2):xend(2),xstart(3):xend(3)) = cmy*vy(mov_yk,xstart(2):xend(2),xstart(3):xend(3)) + cpy*vy(mov_yk+1,xstart(2):xend(2),xstart(3):xend(3))
    call HdfWriteReal2D_X(dsetname,filename,mov_xcut)

    dsetname = trim("vz")//'/'//trim(frame)
    mov_xcut(xstart(2):xend(2),xstart(3):xend(3)) = cmz*vz(mov_zk,xstart(2):xend(2),xstart(3):xend(3)) + cpz*vz(mov_zk+1,xstart(2):xend(2),xstart(3):xend(3))
    call HdfWriteReal2D_X(dsetname,filename,mov_xcut)

    dsetname = trim("pr")//'/'//trim(frame)
    mov_xcut(xstart(2):xend(2),xstart(3):xend(3)) = cmc*pr(mov_ck,xstart(2):xend(2),xstart(3):xend(3)) + cpc*pr(mov_ck+1,xstart(2):xend(2),xstart(3):xend(3))
    call HdfWriteReal2D_X(dsetname,filename,mov_xcut)

    dsetname = trim("temp")//'/'//trim(frame)
    mov_xcut(xstart(2):xend(2),xstart(3):xend(3)) = cmx*temp(mov_xk,xstart(2):xend(2),xstart(3):xend(3)) + cpx*temp(mov_xk+1,xstart(2):xend(2),xstart(3):xend(3))
    call HdfWriteReal2D_X(dsetname,filename,mov_xcut)

    dsetname = trim("co2")//'/'//trim(frame)
    mov_xcut(xstart(2):xend(2),xstart(3):xend(3)) = cmx*co2(mov_xk,xstart(2):xend(2),xstart(3):xend(3)) + cpx*co2(mov_xk+1,xstart(2):xend(2),xstart(3):xend(3))
    call HdfWriteReal2D_X(dsetname,filename,mov_xcut)

    dsetname = trim("h2o")//'/'//trim(frame)
    mov_xcut(xstart(2):xend(2),xstart(3):xend(3)) = cmx*h2o(mov_xk,xstart(2):xend(2),xstart(3):xend(3)) + cpx*h2o(mov_xk+1,xstart(2):xend(2),xstart(3):xend(3))
    call HdfWriteReal2D_X(dsetname,filename,mov_xcut)

end subroutine Movie_xcut

subroutine Movie_ycut

    use param
    use local_arrays, only: vx,vy,vz,pr,temp,co2,h2o
    use decomp_2d, only: xstart,xend
    use movie_indices
    use mpih

    implicit none

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
        mov_ycut(1:nx,xstart(3):xend(3)) = cmx*vx(1:nx,mov_xj,xstart(3):xend(3)) + cpx*vx(1:nx,mov_xj+1,xstart(3):xend(3))
        call HdfWriteReal2D_Y(dsetname,filename,mov_ycut)

        dsetname = trim("vy")//'/'//trim(frame)
        mov_ycut(1:nx,xstart(3):xend(3)) = cmy*vy(1:nx,mov_yj,xstart(3):xend(3)) + cpy*vy(1:nx,mov_yj+1,xstart(3):xend(3))
        call HdfWriteReal2D_Y(dsetname,filename,mov_ycut)

        dsetname = trim("vz")//'/'//trim(frame)
        mov_ycut(1:nx,xstart(3):xend(3)) = cmz*vz(1:nx,mov_zj,xstart(3):xend(3)) + cpz*vz(1:nx,mov_zj+1,xstart(3):xend(3))
        call HdfWriteReal2D_Y(dsetname,filename,mov_ycut)

        dsetname = trim("pr")//'/'//trim(frame)
        mov_ycut(1:nx,xstart(3):xend(3)) = cmc*pr(1:nx,mov_cj,xstart(3):xend(3)) + cpc*pr(1:nx,mov_cj+1,xstart(3):xend(3))
        call HdfWriteReal2D_Y(dsetname,filename,mov_ycut)

        dsetname = trim("temp")//'/'//trim(frame)
        mov_ycut(1:nx,xstart(3):xend(3)) = cmx*temp(1:nx,mov_xj,xstart(3):xend(3)) + cpx*temp(1:nx,mov_xj+1,xstart(3):xend(3))
        call HdfWriteReal2D_Y(dsetname,filename,mov_ycut)

        dsetname = trim("co2")//'/'//trim(frame)
        mov_ycut(1:nx,xstart(3):xend(3)) = cmx*co2(1:nx,mov_xj,xstart(3):xend(3)) + cpx*co2(1:nx,mov_xj+1,xstart(3):xend(3))
        call HdfWriteReal2D_Y(dsetname,filename,mov_ycut)

        dsetname = trim("h2o")//'/'//trim(frame)
        mov_ycut(1:nx,xstart(3):xend(3)) = cmx*h2o(1:nx,mov_xj,xstart(3):xend(3)) + cpx*h2o(1:nx,mov_xj+1,xstart(3):xend(3))
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
        mov_zcut(1:nx,xstart(2):xend(2)) = cmx*vx(1:nx,xstart(2):xend(2),mov_xi) + cpx*vx(1:nx,xstart(2):xend(2),mov_xi+1)
        call HdfWriteReal2D_Z(dsetname,filename,mov_zcut)

        dsetname = trim("vy")//'/'//trim(frame)
        mov_zcut(1:nx,xstart(2):xend(2)) = cmy*vy(1:nx,xstart(2):xend(2),mov_yi) + cpy*vy(1:nx,xstart(2):xend(2),mov_yi+1)
        call HdfWriteReal2D_Z(dsetname,filename,mov_zcut)
    
        dsetname = trim("vz")//'/'//trim(frame)
        mov_zcut(1:nx,xstart(2):xend(2)) = cmz*vz(1:nx,xstart(2):xend(2),mov_zi) + cpz*vz(1:nx,xstart(2):xend(2),mov_zi+1)
        call HdfWriteReal2D_Z(dsetname,filename,mov_zcut)

        dsetname = trim("pr")//'/'//trim(frame)
        mov_zcut(1:nx,xstart(2):xend(2)) = cmc*pr(1:nx,xstart(2):xend(2),mov_ci) + cpc*pr(1:nx,xstart(2):xend(2),mov_ci+1)
        call HdfWriteReal2D_Z(dsetname,filename,mov_zcut)

        dsetname = trim("temp")//'/'//trim(frame)
        mov_zcut(1:nx,xstart(2):xend(2)) = cmx*temp(1:nx,xstart(2):xend(2),mov_xi) + cpx*temp(1:nx,xstart(2):xend(2),mov_xi+1)
        call HdfWriteReal2D_Z(dsetname,filename,mov_zcut)

        dsetname = trim("co2")//'/'//trim(frame)
        mov_zcut(1:nx,xstart(2):xend(2)) = cmx*co2(1:nx,xstart(2):xend(2),mov_xi) + cpx*co2(1:nx,xstart(2):xend(2),mov_xi+1)
        call HdfWriteReal2D_Z(dsetname,filename,mov_zcut)

        dsetname = trim("h2o")//'/'//trim(frame)
        mov_zcut(1:nx,xstart(2):xend(2)) = cmx*h2o(1:nx,xstart(2):xend(2),mov_xi) + cpx*h2o(1:nx,xstart(2):xend(2),mov_xi+1)
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

    character*50        :: dsetname,filename,frame

    if (xend(3).eq.nzm) then

        filename = trim('Results/movie_outlet.h5')
        write (frame,"(i5.5)") nint(time/tframe)
    
        dsetname = trim("vx")//'/'//trim(frame)
        mov_ocut(1:nx,xstart(2):xend(2)) = 0.5d0*(vx(1:nx,xstart(2):xend(2),nzm) + vx(1:nx,xstart(2):xend(2),nz))
        call HdfWriteReal2D_Z(dsetname,filename,mov_ocut)

        dsetname = trim("vy")//'/'//trim(frame)
        mov_ocut(1:nx,xstart(2):xend(2)) = 0.5d0*(vy(1:nx,xstart(2):xend(2),nzm) + vy(1:nx,xstart(2):xend(2),nz))
        call HdfWriteReal2D_Z(dsetname,filename,mov_ocut)
    
        dsetname = trim("vz")//'/'//trim(frame)
        mov_ocut(1:nx,xstart(2):xend(2)) = vz(1:nx,xstart(2):xend(2),nz)
        call HdfWriteReal2D_Z(dsetname,filename,mov_ocut)

        dsetname = trim("pr")//'/'//trim(frame)
        mov_ocut(1:nx,xstart(2):xend(2)) = 0.5d0*(pr(1:nx,xstart(2):xend(2),nzm) + pr(1:nx,xstart(2):xend(2),nz))
        call HdfWriteReal2D_Z(dsetname,filename,mov_ocut)

        dsetname = trim("temp")//'/'//trim(frame)
        mov_ocut(1:nx,xstart(2):xend(2)) = 0.5d0*(temp(1:nx,xstart(2):xend(2),nzm) + temp(1:nx,xstart(2):xend(2),nz))
        call HdfWriteReal2D_Z(dsetname,filename,mov_ocut)

        dsetname = trim("co2")//'/'//trim(frame)
        mov_ocut(1:nx,xstart(2):xend(2)) = 0.5d0*(co2(1:nx,xstart(2):xend(2),nzm) + co2(1:nx,xstart(2):xend(2),nz))
        call HdfWriteReal2D_Z(dsetname,filename,mov_ocut)

        dsetname = trim("h2o")//'/'//trim(frame)
        mov_ocut(1:nx,xstart(2):xend(2)) = 0.5d0*(h2o(1:nx,xstart(2):xend(2),nzm) + h2o(1:nx,xstart(2):xend(2),nz))
        call HdfWriteReal2D_Z(dsetname,filename,mov_ocut)

    end if

end subroutine Movie_outlet
