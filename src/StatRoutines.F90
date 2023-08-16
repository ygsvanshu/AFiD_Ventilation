!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: StatRoutines.F90                               !
!    CONTAINS: subroutine CalcStats,WriteStats            !
!                                                         ! 
!    PURPOSE: Calculates and writes out statistics for    !
!     the flow field. All quantities are averaged in   	  !
!	  the two horizontal (homogeneous) directions.   	  !
!														  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InitStats

	use param
	use decomp_2d, only: xstart,xend
	use stat_arrays
	use mpih

    implicit none

	integer :: i,j,k

	! For the X grid
	do k=1,nx-1
		if ((xc(k).le.stats2Dx).and.(xc(k+1).gt.stats2Dx)) stat_xk = k
	end do
	do j=1,nym-1
		if ((ym(j).le.stats2Dy).and.(ym(j+1).gt.stats2Dy)) stat_xj = j
	end do
	do i=1,nzm-1
		if ((zm(i).le.stats2Dz).and.(zm(i+1).gt.stats2Dz)) stat_xi = i
	end do

	! For the Y grid
	do k=1,nxm-1
		if ((xm(k).le.stats2Dx).and.(xm(k+1).gt.stats2Dx)) stat_yk = k
	end do
	do j=1,ny-1
		if ((yc(j).le.stats2Dy).and.(yc(j+1).gt.stats2Dy)) stat_yj = j
	end do
	do i=1,nzm-1
		if ((zm(i).le.stats2Dz).and.(zm(i+1).gt.stats2Dz)) stat_yi = i
	end do

	! For the Z grid
	do k=1,nxm-1
		if ((xm(k).le.stats2Dx).and.(xm(k+1).gt.stats2Dx)) stat_zk = k
	end do
	do j=1,nym-1
		if ((ym(j).le.stats2Dy).and.(ym(j+1).gt.stats2Dy)) stat_zj = j
	end do
	do i=1,nz-1
		if ((zc(i).le.stats2Dz).and.(zc(i+1).gt.stats2Dz)) stat_zi = i
	end do

	nstatsamples 	= 0
	tstat 			= 0.0d0
	tinterval 		= 0.0d0
	
	vx_m1_xcut(:,:) = 0.0d0
	vx_m2_xcut(:,:) = 0.0d0
	vx_m3_xcut(:,:) = 0.0d0
	vx_m4_xcut(:,:) = 0.0d0

	vy_m1_xcut(:,:) = 0.0d0
	vy_m2_xcut(:,:) = 0.0d0
	vy_m3_xcut(:,:) = 0.0d0
	vy_m4_xcut(:,:) = 0.0d0

	vz_m1_xcut(:,:) = 0.0d0
	vz_m2_xcut(:,:) = 0.0d0
	vz_m3_xcut(:,:) = 0.0d0
	vz_m4_xcut(:,:) = 0.0d0

	temp_m1_xcut(:,:) = 0.0d0
	temp_m2_xcut(:,:) = 0.0d0
	temp_m3_xcut(:,:) = 0.0d0
	temp_m4_xcut(:,:) = 0.0d0

	co2_m1_xcut(:,:) = 0.0d0
	co2_m2_xcut(:,:) = 0.0d0
	co2_m3_xcut(:,:) = 0.0d0
	co2_m4_xcut(:,:) = 0.0d0

	h2o_m1_xcut(:,:) = 0.0d0
	h2o_m2_xcut(:,:) = 0.0d0
	h2o_m3_xcut(:,:) = 0.0d0
	h2o_m4_xcut(:,:) = 0.0d0

end subroutine InitStats

subroutine ReadStats

	use mpih
	use param
	use decomp_2d, only: xstart,xend
	use stat_arrays
	use hdf5

	implicit none

	character*50	:: filename,dsetname
	logical			:: fexist

	filename = trim('Results/Stats/stafield_master.h5')
	inquire(file=filename,exist=fexist)
	if (fexist) then 

		if (ismaster) then
			dsetname = trim('averaging_samples')
			call HdfSerialReadIntScalar(dsetname,filename,nstatsamples)
			dsetname = trim('averaging_time')
			call HdfSerialReadRealScalar(dsetname,filename,tstat)
			dsetname = trim('averaging_interval')
			call HdfSerialReadRealScalar(dsetname,filename,tinterval)
		end if

		call MpiBarrier
		call MpiBcastInt(nstatsamples)
		call MpiBcastReal(tstat)
		call MpiBcastReal(tinterval)

		! X-CUT

		dsetname = trim('/xcut/vx_m1')
		call HdfReadReal2D_X(filename,dsetname,vx_m1_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)
		dsetname = trim('/xcut/vy_m1')
		call HdfReadReal2D_X(filename,dsetname,vy_m1_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)
		dsetname = trim('/xcut/vz_m1')
		call HdfReadReal2D_X(filename,dsetname,vz_m1_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)
		dsetname = trim('/xcut/temp_m1')
		call HdfReadReal2D_X(filename,dsetname,temp_m1_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)
		dsetname = trim('/xcut/co2_m1')
		call HdfReadReal2D_X(filename,dsetname,co2_m1_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)
		dsetname = trim('/xcut/h2o_m1')
		call HdfReadReal2D_X(filename,dsetname,h2o_m1_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)

		dsetname = trim('/xcut/vx_m2')
		call HdfReadReal2D_X(filename,dsetname,vx_m2_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)
		dsetname = trim('/xcut/vy_m2')
		call HdfReadReal2D_X(filename,dsetname,vy_m2_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)
		dsetname = trim('/xcut/vz_m2')
		call HdfReadReal2D_X(filename,dsetname,vz_m2_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)
		dsetname = trim('/xcut/temp_m2')
		call HdfReadReal2D_X(filename,dsetname,temp_m2_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)
		dsetname = trim('/xcut/co2_m2')
		call HdfReadReal2D_X(filename,dsetname,co2_m2_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)
		dsetname = trim('/xcut/h2o_m2')
		call HdfReadReal2D_X(filename,dsetname,h2o_m2_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)

		dsetname = trim('/xcut/vx_m3')
		call HdfReadReal2D_X(filename,dsetname,vx_m3_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)
		dsetname = trim('/xcut/vy_m3')
		call HdfReadReal2D_X(filename,dsetname,vy_m3_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)
		dsetname = trim('/xcut/vz_m3')
		call HdfReadReal2D_X(filename,dsetname,vz_m3_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)
		dsetname = trim('/xcut/temp_m3')
		call HdfReadReal2D_X(filename,dsetname,temp_m3_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)
		dsetname = trim('/xcut/co2_m3')
		call HdfReadReal2D_X(filename,dsetname,co2_m3_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)
		dsetname = trim('/xcut/h2o_m3')
		call HdfReadReal2D_X(filename,dsetname,h2o_m3_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)

		dsetname = trim('/xcut/vx_m4')
		call HdfReadReal2D_X(filename,dsetname,vx_m4_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)
		dsetname = trim('/xcut/vy_m4')
		call HdfReadReal2D_X(filename,dsetname,vy_m4_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)
		dsetname = trim('/xcut/vz_m4')
		call HdfReadReal2D_X(filename,dsetname,vz_m4_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)
		dsetname = trim('/xcut/temp_m4')
		call HdfReadReal2D_X(filename,dsetname,temp_m4_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)
		dsetname = trim('/xcut/co2_m4')
		call HdfReadReal2D_X(filename,dsetname,co2_m4_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)
		dsetname = trim('/xcut/h2o_m4')
		call HdfReadReal2D_X(filename,dsetname,h2o_m4_xcut,xstart(2),xend(2),xstart(3),xend(3),ny,nz)

		! Y-CUT

		if ((yc(xstart(2)).le.stats2Dy).and.(yc(xend(2)+lvlhalo).gt.stats2Dy)) then

			dsetname = trim('/ycut/vx_m1')
			call HdfReadReal2D_Y(filename,dsetname,vx_m1_ycut,xstart(3),xend(3),nx,nz)
			dsetname = trim('/ycut/vy_m1')
			call HdfReadReal2D_Y(filename,dsetname,vy_m1_ycut,xstart(3),xend(3),nx,nz)
			dsetname = trim('/ycut/vz_m1')
			call HdfReadReal2D_Y(filename,dsetname,vz_m1_ycut,xstart(3),xend(3),nx,nz)
			dsetname = trim('/ycut/temp_m1')
			call HdfReadReal2D_Y(filename,dsetname,temp_m1_ycut,xstart(3),xend(3),nx,nz)
			dsetname = trim('/ycut/co2_m1')
			call HdfReadReal2D_Y(filename,dsetname,co2_m1_ycut,xstart(3),xend(3),nx,nz)
			dsetname = trim('/ycut/h2o_m1')
			call HdfReadReal2D_Y(filename,dsetname,h2o_m1_ycut,xstart(3),xend(3),nx,nz)

			dsetname = trim('/ycut/vx_m2')
			call HdfReadReal2D_Y(filename,dsetname,vx_m2_ycut,xstart(3),xend(3),nx,nz)
			dsetname = trim('/ycut/vy_m2')
			call HdfReadReal2D_Y(filename,dsetname,vy_m2_ycut,xstart(3),xend(3),nx,nz)
			dsetname = trim('/ycut/vz_m2')
			call HdfReadReal2D_Y(filename,dsetname,vz_m2_ycut,xstart(3),xend(3),nx,nz)
			dsetname = trim('/ycut/temp_m2')
			call HdfReadReal2D_Y(filename,dsetname,temp_m2_ycut,xstart(3),xend(3),nx,nz)
			dsetname = trim('/ycut/co2_m2')
			call HdfReadReal2D_Y(filename,dsetname,co2_m2_ycut,xstart(3),xend(3),nx,nz)
			dsetname = trim('/ycut/h2o_m2')
			call HdfReadReal2D_Y(filename,dsetname,h2o_m2_ycut,xstart(3),xend(3),nx,nz)

			dsetname = trim('/ycut/vx_m3')
			call HdfReadReal2D_Y(filename,dsetname,vx_m3_ycut,xstart(3),xend(3),nx,nz)
			dsetname = trim('/ycut/vy_m3')
			call HdfReadReal2D_Y(filename,dsetname,vy_m3_ycut,xstart(3),xend(3),nx,nz)
			dsetname = trim('/ycut/vz_m3')
			call HdfReadReal2D_Y(filename,dsetname,vz_m3_ycut,xstart(3),xend(3),nx,nz)
			dsetname = trim('/ycut/temp_m3')
			call HdfReadReal2D_Y(filename,dsetname,temp_m3_ycut,xstart(3),xend(3),nx,nz)
			dsetname = trim('/ycut/co2_m3')
			call HdfReadReal2D_Y(filename,dsetname,co2_m3_ycut,xstart(3),xend(3),nx,nz)
			dsetname = trim('/ycut/h2o_m3')
			call HdfReadReal2D_Y(filename,dsetname,h2o_m3_ycut,xstart(3),xend(3),nx,nz)

			dsetname = trim('/ycut/vx_m4')
			call HdfReadReal2D_Y(filename,dsetname,vx_m4_ycut,xstart(3),xend(3),nx,nz)
			dsetname = trim('/ycut/vy_m4')
			call HdfReadReal2D_Y(filename,dsetname,vy_m4_ycut,xstart(3),xend(3),nx,nz)
			dsetname = trim('/ycut/vz_m4')
			call HdfReadReal2D_Y(filename,dsetname,vz_m4_ycut,xstart(3),xend(3),nx,nz)
			dsetname = trim('/ycut/temp_m4')
			call HdfReadReal2D_Y(filename,dsetname,temp_m4_ycut,xstart(3),xend(3),nx,nz)
			dsetname = trim('/ycut/co2_m4')
			call HdfReadReal2D_Y(filename,dsetname,co2_m4_ycut,xstart(3),xend(3),nx,nz)
			dsetname = trim('/ycut/h2o_m4')
			call HdfReadReal2D_Y(filename,dsetname,h2o_m4_ycut,xstart(3),xend(3),nx,nz)

		end if

		! Z-CUT

		if ((zc(xstart(3)).le.stats2Dz).and.(zc(xend(3)+lvlhalo).gt.stats2Dz)) then

			dsetname = trim('/zcut/vx_m1')
			call HdfReadReal2D_Z(filename,dsetname,vx_m1_zcut,xstart(2),xend(2),nx,ny)
			dsetname = trim('/zcut/vy_m1')
			call HdfReadReal2D_Z(filename,dsetname,vy_m1_zcut,xstart(2),xend(2),nx,ny)
			dsetname = trim('/zcut/vz_m1')
			call HdfReadReal2D_Z(filename,dsetname,vz_m1_zcut,xstart(2),xend(2),nx,ny)
			dsetname = trim('/zcut/temp_m1')
			call HdfReadReal2D_Z(filename,dsetname,temp_m1_zcut,xstart(2),xend(2),nx,ny)
			dsetname = trim('/zcut/co2_m1')
			call HdfReadReal2D_Z(filename,dsetname,co2_m1_zcut,xstart(2),xend(2),nx,ny)
			dsetname = trim('/zcut/h2o_m1')
			call HdfReadReal2D_Z(filename,dsetname,h2o_m1_zcut,xstart(2),xend(2),nx,ny)

			dsetname = trim('/zcut/vx_m2')
			call HdfReadReal2D_Z(filename,dsetname,vx_m2_zcut,xstart(2),xend(2),nx,ny)
			dsetname = trim('/zcut/vy_m2')
			call HdfReadReal2D_Z(filename,dsetname,vy_m2_zcut,xstart(2),xend(2),nx,ny)
			dsetname = trim('/zcut/vz_m2')
			call HdfReadReal2D_Z(filename,dsetname,vz_m2_zcut,xstart(2),xend(2),nx,ny)
			dsetname = trim('/zcut/temp_m2')
			call HdfReadReal2D_Z(filename,dsetname,temp_m2_zcut,xstart(2),xend(2),nx,ny)
			dsetname = trim('/zcut/co2_m2')
			call HdfReadReal2D_Z(filename,dsetname,co2_m2_zcut,xstart(2),xend(2),nx,ny)
			dsetname = trim('/zcut/h2o_m2')
			call HdfReadReal2D_Z(filename,dsetname,h2o_m2_zcut,xstart(2),xend(2),nx,ny)

			dsetname = trim('/zcut/vx_m3')
			call HdfReadReal2D_Z(filename,dsetname,vx_m3_zcut,xstart(2),xend(2),nx,ny)
			dsetname = trim('/zcut/vy_m3')
			call HdfReadReal2D_Z(filename,dsetname,vy_m3_zcut,xstart(2),xend(2),nx,ny)
			dsetname = trim('/zcut/vz_m3')
			call HdfReadReal2D_Z(filename,dsetname,vz_m3_zcut,xstart(2),xend(2),nx,ny)
			dsetname = trim('/zcut/temp_m3')
			call HdfReadReal2D_Z(filename,dsetname,temp_m3_zcut,xstart(2),xend(2),nx,ny)
			dsetname = trim('/zcut/co2_m3')
			call HdfReadReal2D_Z(filename,dsetname,co2_m3_zcut,xstart(2),xend(2),nx,ny)
			dsetname = trim('/zcut/h2o_m3')
			call HdfReadReal2D_Z(filename,dsetname,h2o_m3_zcut,xstart(2),xend(2),nx,ny)

			dsetname = trim('/zcut/vx_m4')
			call HdfReadReal2D_Z(filename,dsetname,vx_m4_zcut,xstart(2),xend(2),nx,ny)
			dsetname = trim('/zcut/vy_m4')
			call HdfReadReal2D_Z(filename,dsetname,vy_m4_zcut,xstart(2),xend(2),nx,ny)
			dsetname = trim('/zcut/vz_m4')
			call HdfReadReal2D_Z(filename,dsetname,vz_m4_zcut,xstart(2),xend(2),nx,ny)
			dsetname = trim('/zcut/temp_m4')
			call HdfReadReal2D_Z(filename,dsetname,temp_m4_zcut,xstart(2),xend(2),nx,ny)
			dsetname = trim('/zcut/co2_m4')
			call HdfReadReal2D_Z(filename,dsetname,co2_m4_zcut,xstart(2),xend(2),nx,ny)
			dsetname = trim('/zcut/h2o_m4')
			call HdfReadReal2D_Z(filename,dsetname,h2o_m4_zcut,xstart(2),xend(2),nx,ny)

		end if
		
	else
		if(ismaster) write(6,*) 'Unable to read statistical files'
		if(ismaster) write(6,*) 'Restarting statistics from zero'
	end if

end subroutine ReadStats

subroutine CalcStats

	use param
	use local_arrays, only: vz,vy,vx,temp,co2,h2o
	use decomp_2d, only: xstart,xend
	use stat_arrays
	use mpih

    implicit none

	integer :: i,j,k
	real    :: cpx,cmx,cpy,cmy,cpz,cmz

	nstatsamples = nstatsamples + 1
	tstat        = tstat        + dt

	! For X-CUT

	cpx = (stats2Dx      - xc(stat_xk))/(xc(stat_xk+1) - xc(stat_xk))
	cmx = (xc(stat_xk+1) - stats2Dx)   /(xc(stat_xk+1) - xc(stat_xk))
	cpy = (stats2Dx      - xm(stat_yk))/(xm(stat_yk+1) - xm(stat_yk))
	cmy = (xm(stat_yk+1) - stats2Dx)   /(xm(stat_yk+1) - xm(stat_yk))
	cpz = (stats2Dx      - xm(stat_zk))/(xm(stat_zk+1) - xm(stat_zk))
	cmz = (xm(stat_zk+1) - stats2Dx)   /(xm(stat_zk+1) - xm(stat_zk))

	do j=xstart(2),xend(2)
		do i=xstart(3),xend(3)
			vx_m1_xcut(j,i)   = vx_m1_xcut(j,i)   + dt*(cmx*vx(stat_xk,j,i)      + cpx*vx(stat_xk+1,j,i))
			vy_m1_xcut(j,i)   = vy_m1_xcut(j,i)   + dt*(cmy*vy(stat_yk,j,i)      + cpy*vy(stat_yk+1,j,i))
			vz_m1_xcut(j,i)   = vz_m1_xcut(j,i)   + dt*(cmz*vz(stat_zk,j,i)      + cpz*vz(stat_zk+1,j,i))
			temp_m1_xcut(j,i) = temp_m1_xcut(j,i) + dt*(cmx*temp(stat_xk,j,i)    + cpx*temp(stat_xk+1,j,i))
			co2_m1_xcut(j,i)  = co2_m1_xcut(j,i)  + dt*(cmx*co2(stat_xk,j,i)     + cpx*co2(stat_xk+1,j,i))
			h2o_m1_xcut(j,i)  = h2o_m1_xcut(j,i)  + dt*(cmx*h2o(stat_xk,j,i)     + cpx*h2o(stat_xk+1,j,i))

			vx_m2_xcut(j,i)   = vx_m2_xcut(j,i)   + dt*(cmx*vx(stat_xk,j,i)**2   + cpx*vx(stat_xk+1,j,i)**2)
			vy_m2_xcut(j,i)   = vy_m2_xcut(j,i)   + dt*(cmy*vy(stat_yk,j,i)**2   + cpy*vy(stat_yk+1,j,i)**2)
			vz_m2_xcut(j,i)   = vz_m2_xcut(j,i)   + dt*(cmz*vz(stat_zk,j,i)**2   + cpz*vz(stat_zk+1,j,i)**2)
			temp_m2_xcut(j,i) = temp_m2_xcut(j,i) + dt*(cmx*temp(stat_xk,j,i)**2 + cpx*temp(stat_xk+1,j,i)**2)
			co2_m2_xcut(j,i)  = co2_m2_xcut(j,i)  + dt*(cmx*co2(stat_xk,j,i)**2  + cpx*co2(stat_xk+1,j,i)**2)
			h2o_m2_xcut(j,i)  = h2o_m2_xcut(j,i)  + dt*(cmx*h2o(stat_xk,j,i)**2  + cpx*h2o(stat_xk+1,j,i)**2)

			vx_m3_xcut(j,i)   = vx_m3_xcut(j,i)   + dt*(cmx*vx(stat_xk,j,i)**3   + cpx*vx(stat_xk+1,j,i)**3)
			vy_m3_xcut(j,i)   = vy_m3_xcut(j,i)   + dt*(cmy*vy(stat_yk,j,i)**3   + cpy*vy(stat_yk+1,j,i)**3)
			vz_m3_xcut(j,i)   = vz_m3_xcut(j,i)   + dt*(cmz*vz(stat_zk,j,i)**3   + cpz*vz(stat_zk+1,j,i)**3)
			temp_m3_xcut(j,i) = temp_m3_xcut(j,i) + dt*(cmx*temp(stat_xk,j,i)**3 + cpx*temp(stat_xk+1,j,i)**3)
			co2_m3_xcut(j,i)  = co2_m3_xcut(j,i)  + dt*(cmx*co2(stat_xk,j,i)**3  + cpx*co2(stat_xk+1,j,i)**3)
			h2o_m3_xcut(j,i)  = h2o_m3_xcut(j,i)  + dt*(cmx*h2o(stat_xk,j,i)**3  + cpx*h2o(stat_xk+1,j,i)**3)

			vx_m4_xcut(j,i)   = vx_m4_xcut(j,i)   + dt*(cmx*vx(stat_xk,j,i)**4   + cpx*vx(stat_xk+1,j,i)**4)
			vy_m4_xcut(j,i)   = vy_m4_xcut(j,i)   + dt*(cmy*vy(stat_yk,j,i)**4   + cpy*vy(stat_yk+1,j,i)**4)
			vz_m4_xcut(j,i)   = vz_m4_xcut(j,i)   + dt*(cmz*vz(stat_zk,j,i)**4   + cpz*vz(stat_zk+1,j,i)**4)
			temp_m4_xcut(j,i) = temp_m4_xcut(j,i) + dt*(cmx*temp(stat_xk,j,i)**4 + cpx*temp(stat_xk+1,j,i)**4)
			co2_m4_xcut(j,i)  = co2_m4_xcut(j,i)  + dt*(cmx*co2(stat_xk,j,i)**4  + cpx*co2(stat_xk+1,j,i)**4)
			h2o_m4_xcut(j,i)  = h2o_m4_xcut(j,i)  + dt*(cmx*h2o(stat_xk,j,i)**4  + cpx*h2o(stat_xk+1,j,i)**4)
		end do
	end do

	! For Y-CUT

	if ((yc(xstart(2)).le.stats2Dy).and.(yc(xend(2)+lvlhalo).gt.stats2Dy)) then

		cpx = (stats2Dy      - ym(stat_xj))/(ym(stat_xj+1) - ym(stat_xj))
		cmx = (ym(stat_xj+1) - stats2Dy)   /(ym(stat_xj+1) - ym(stat_xj))
		cpy = (stats2Dy      - yc(stat_yj))/(yc(stat_yj+1) - yc(stat_yj))
		cmy = (yc(stat_yj+1) - stats2Dy)   /(yc(stat_yj+1) - yc(stat_yj))
		cpz = (stats2Dy      - ym(stat_zj))/(ym(stat_zj+1) - ym(stat_zj))
		cmz = (ym(stat_zj+1) - stats2Dy)   /(ym(stat_zj+1) - ym(stat_zj))

		do k=1,nx
			do i=xstart(3),xend(3)
				vx_m1_ycut(k,i)   = vx_m1_ycut(k,i)   + dt*(cmx*vx(k,stat_xj,i)      + cpx*vx(k,stat_xj+1,i))
				vy_m1_ycut(k,i)   = vy_m1_ycut(k,i)   + dt*(cmy*vy(k,stat_yj,i)      + cpy*vy(k,stat_yj+1,i))
				vz_m1_ycut(k,i)   = vz_m1_ycut(k,i)   + dt*(cmz*vz(k,stat_zj,i)      + cpz*vz(k,stat_zj+1,i))
				temp_m1_ycut(k,i) = temp_m1_ycut(k,i) + dt*(cmx*temp(k,stat_xj,i)    + cpx*temp(k,stat_xj+1,i))
				co2_m1_ycut(k,i)  = co2_m1_ycut(k,i)  + dt*(cmx*co2(k,stat_xj,i)     + cpx*co2(k,stat_xj+1,i))
				h2o_m1_ycut(k,i)  = h2o_m1_ycut(k,i)  + dt*(cmx*h2o(k,stat_xj,i)     + cpx*h2o(k,stat_xj+1,i))

				vx_m2_ycut(k,i)   = vx_m2_ycut(k,i)   + dt*(cmx*vx(k,stat_xj,i)**2   + cpx*vx(k,stat_xj+1,i)**2)
				vy_m2_ycut(k,i)   = vy_m2_ycut(k,i)   + dt*(cmy*vy(k,stat_yj,i)**2   + cpy*vy(k,stat_yj+1,i)**2)
				vz_m2_ycut(k,i)   = vz_m2_ycut(k,i)   + dt*(cmz*vz(k,stat_zj,i)**2   + cpz*vz(k,stat_zj+1,i)**2)
				temp_m2_ycut(k,i) = temp_m2_ycut(k,i) + dt*(cmx*temp(k,stat_xj,i)**2 + cpx*temp(k,stat_xj+1,i)**2)
				co2_m2_ycut(k,i)  = co2_m2_ycut(k,i)  + dt*(cmx*co2(k,stat_xj,i)**2  + cpx*co2(k,stat_xj+1,i)**2)
				h2o_m2_ycut(k,i)  = h2o_m2_ycut(k,i)  + dt*(cmx*h2o(k,stat_xj,i)**2  + cpx*h2o(k,stat_xj+1,i)**2)

				vx_m3_ycut(k,i)   = vx_m3_ycut(k,i)   + dt*(cmx*vx(k,stat_xj,i)**3   + cpx*vx(k,stat_xj+1,i)**3)
				vy_m3_ycut(k,i)   = vy_m3_ycut(k,i)   + dt*(cmy*vy(k,stat_yj,i)**3   + cpy*vy(k,stat_yj+1,i)**3)
				vz_m3_ycut(k,i)   = vz_m3_ycut(k,i)   + dt*(cmz*vz(k,stat_zj,i)**3   + cpz*vz(k,stat_zj+1,i)**3)
				temp_m3_ycut(k,i) = temp_m3_ycut(k,i) + dt*(cmx*temp(k,stat_xj,i)**3 + cpx*temp(k,stat_xj+1,i)**3)
				co2_m3_ycut(k,i)  = co2_m3_ycut(k,i)  + dt*(cmx*co2(k,stat_xj,i)**3  + cpx*co2(k,stat_xj+1,i)**3)
				h2o_m3_ycut(k,i)  = h2o_m3_ycut(k,i)  + dt*(cmx*h2o(k,stat_xj,i)**3  + cpx*h2o(k,stat_xj+1,i)**3)

				vx_m4_ycut(k,i)   = vx_m4_ycut(k,i)   + dt*(cmx*vx(k,stat_xj,i)**4   + cpx*vx(k,stat_xj+1,i)**4)
				vy_m4_ycut(k,i)   = vy_m4_ycut(k,i)   + dt*(cmy*vy(k,stat_yj,i)**4   + cpy*vy(k,stat_yj+1,i)**4)
				vz_m4_ycut(k,i)   = vz_m4_ycut(k,i)   + dt*(cmz*vz(k,stat_zj,i)**4   + cpz*vz(k,stat_zj+1,i)**4)
				temp_m4_ycut(k,i) = temp_m4_ycut(k,i) + dt*(cmx*temp(k,stat_xj,i)**4 + cpx*temp(k,stat_xj+1,i)**4)
				co2_m4_ycut(k,i)  = co2_m4_ycut(k,i)  + dt*(cmx*co2(k,stat_xj,i)**4  + cpx*co2(k,stat_xj+1,i)**4)
				h2o_m4_ycut(k,i)  = h2o_m4_ycut(k,i)  + dt*(cmx*h2o(k,stat_xj,i)**4  + cpx*h2o(k,stat_xj+1,i)**4)
			end do
		end do
	end if
	
	! For Z-CUT

	if ((zc(xstart(3)).le.stats2Dz).and.(zc(xend(3)+lvlhalo).gt.stats2Dz)) then

		cpx = (stats2Dz      - zm(stat_xi))/(zm(stat_xi+1) - zm(stat_xi))
		cmx = (zm(stat_xi+1) - stats2Dz)   /(zm(stat_xi+1) - zm(stat_xi))
		cpy = (stats2Dz      - zm(stat_yi))/(zm(stat_yi+1) - zm(stat_yi))
		cmy = (zm(stat_yi+1) - stats2Dz)   /(zm(stat_yi+1) - zm(stat_yi))
		cpz = (stats2Dz      - zc(stat_zi))/(zc(stat_zi+1) - zc(stat_zi))
		cmz = (zc(stat_zi+1) - stats2Dz)   /(zc(stat_zi+1) - zc(stat_zi))

		do k=1,nx
			do j=xstart(2),xend(2)
				vx_m1_zcut(k,j)   = vx_m1_zcut(k,j)   + dt*(cmx*vx(k,j,stat_xi)      + cpx*vx(k,j,stat_xi+1))
				vy_m1_zcut(k,j)   = vy_m1_zcut(k,j)   + dt*(cmy*vy(k,j,stat_yi)      + cpy*vy(k,j,stat_yi+1))
				vz_m1_zcut(k,j)   = vz_m1_zcut(k,j)   + dt*(cmz*vz(k,j,stat_zi)      + cpz*vz(k,j,stat_zi+1))
				temp_m1_zcut(k,j) = temp_m1_zcut(k,j) + dt*(cmx*temp(k,j,stat_xi)    + cpx*temp(k,j,stat_xi+1))
				co2_m1_zcut(k,j)  = co2_m1_zcut(k,j)  + dt*(cmx*co2(k,j,stat_xi)     + cpx*co2(k,j,stat_xi+1))
				h2o_m1_zcut(k,j)  = h2o_m1_zcut(k,j)  + dt*(cmx*h2o(k,j,stat_xi)     + cpx*h2o(k,j,stat_xi+1))

				vx_m2_zcut(k,j)   = vx_m2_zcut(k,j)   + dt*(cmx*vx(k,j,stat_xi)**2   + cpx*vx(k,j,stat_xi+1)**2)
				vy_m2_zcut(k,j)   = vy_m2_zcut(k,j)   + dt*(cmy*vy(k,j,stat_yi)**2   + cpy*vy(k,j,stat_yi+1)**2)
				vz_m2_zcut(k,j)   = vz_m2_zcut(k,j)   + dt*(cmz*vz(k,j,stat_zi)**2   + cpz*vz(k,j,stat_zi+1)**2)
				temp_m2_zcut(k,j) = temp_m2_zcut(k,j) + dt*(cmx*temp(k,j,stat_xi)**2 + cpx*temp(k,j,stat_xi+1)**2)
				co2_m2_zcut(k,j)  = co2_m2_zcut(k,j)  + dt*(cmx*co2(k,j,stat_xi)**2  + cpx*co2(k,j,stat_xi+1)**2)
				h2o_m2_zcut(k,j)  = h2o_m2_zcut(k,j)  + dt*(cmx*h2o(k,j,stat_xi)**2  + cpx*h2o(k,j,stat_xi+1)**2)

				vx_m3_zcut(k,j)   = vx_m3_zcut(k,j)   + dt*(cmx*vx(k,j,stat_xi)**3   + cpx*vx(k,j,stat_xi+1)**3)
				vy_m3_zcut(k,j)   = vy_m3_zcut(k,j)   + dt*(cmy*vy(k,j,stat_yi)**3   + cpy*vy(k,j,stat_yi+1)**3)
				vz_m3_zcut(k,j)   = vz_m3_zcut(k,j)   + dt*(cmz*vz(k,j,stat_zi)**3   + cpz*vz(k,j,stat_zi+1)**3)
				temp_m3_zcut(k,j) = temp_m3_zcut(k,j) + dt*(cmx*temp(k,j,stat_xi)**3 + cpx*temp(k,j,stat_xi+1)**3)
				co2_m3_zcut(k,j)  = co2_m3_zcut(k,j)  + dt*(cmx*co2(k,j,stat_xi)**3  + cpx*co2(k,j,stat_xi+1)**3)
				h2o_m3_zcut(k,j)  = h2o_m3_zcut(k,j)  + dt*(cmx*h2o(k,j,stat_xi)**3  + cpx*h2o(k,j,stat_xi+1)**3)

				vx_m4_zcut(k,j)   = vx_m4_zcut(k,j)   + dt*(cmx*vx(k,j,stat_xi)**4   + cpx*vx(k,j,stat_xi+1)**4)
				vy_m4_zcut(k,j)   = vy_m4_zcut(k,j)   + dt*(cmy*vy(k,j,stat_yi)**4   + cpy*vy(k,j,stat_yi+1)**4)
				vz_m4_zcut(k,j)   = vz_m4_zcut(k,j)   + dt*(cmz*vz(k,j,stat_zi)**4   + cpz*vz(k,j,stat_zi+1)**4)
				temp_m4_zcut(k,j) = temp_m4_zcut(k,j) + dt*(cmx*temp(k,j,stat_xi)**4 + cpx*temp(k,j,stat_xi+1)**4)
				co2_m4_zcut(k,j)  = co2_m4_zcut(k,j)  + dt*(cmx*co2(k,j,stat_xi)**4  + cpx*co2(k,j,stat_xi+1)**4)
				h2o_m4_zcut(k,j)  = h2o_m4_zcut(k,j)  + dt*(cmx*h2o(k,j,stat_xi)**4  + cpx*h2o(k,j,stat_xi+1)**4)
			end do
		end do
	end if

	return

end subroutine CalcStats

subroutine WriteStats(filename)

	use mpih
	use param
	use decomp_2d, only: xstart,xend
	use stat_arrays
	use hdf5
	
	implicit none

	character*50, intent(inout) :: filename
	character*50				:: linkname,dsetname
	
	if (ismaster) then

		dsetname = trim('averaging_samples')
		call HdfSerialWriteIntScalar(dsetname,filename,nstatsamples)
		dsetname = trim('averaging_time')
		call HdfSerialWriteRealScalar(dsetname,filename,tstat)
		dsetname = trim('averaging_interval')
		call HdfSerialWriteRealScalar(dsetname,filename,tinterval+time-tsta)
		dsetname = trim('sample_interval')
		call HdfSerialWriteIntScalar(dsetname,filename,nout)

		dsetname = trim('xc')
		call HdfSerialWriteReal1D(dsetname,filename,xc,nx)
		dsetname = trim('yc')
		call HdfSerialWriteReal1D(dsetname,filename,yc,ny)
		dsetname = trim('zc')
		call HdfSerialWriteReal1D(dsetname,filename,zc,nz)

		dsetname = trim('xm')
		call HdfSerialWriteReal1D(dsetname,filename,xm,nxm)
		dsetname = trim('ym')
		call HdfSerialWriteReal1D(dsetname,filename,ym,nym)
		dsetname = trim('zm')
		call HdfSerialWriteReal1D(dsetname,filename,zm,nzm)

		dsetname = trim('Rayleigh Number')
		call HdfSerialWriteRealScalar(dsetname,filename,ray)
		dsetname = trim('Prandtl Number')
		call HdfSerialWriteRealScalar(dsetname,filename,pra)
		dsetname = trim('CO2 Expansion Ratio')
		call HdfSerialWriteRealScalar(dsetname,filename,lambda_co2)
		dsetname = trim('H2O Expansion Ratio')
		call HdfSerialWriteRealScalar(dsetname,filename,lambda_h2o)
		dsetname = trim('Inlet Velocity')
		call HdfSerialWriteRealScalar(dsetname,filename,ivel)
		dsetname = trim('Inlet Dimension')
		call HdfSerialWriteRealScalar(dsetname,filename,ilen)
		dsetname = trim('Inlet Height')
		call HdfSerialWriteRealScalar(dsetname,filename,iheight)
		dsetname = trim('Outlet Dimension')
		call HdfSerialWriteRealScalar(dsetname,filename,olen)
		dsetname = trim('Outlet Height')
		call HdfSerialWriteRealScalar(dsetname,filename,oheight)

	end if
	
	! X-CUT

	linkname = trim('/xcut')
	call HdfParallelCreateGroup(linkname,filename,comm_xcut)

	dsetname = trim('/xcut/vx_m1')
	call HdfWriteReal2D_X(dsetname,filename,vx_m1_xcut)
	dsetname = trim('/xcut/vy_m1')
	call HdfWriteReal2D_X(dsetname,filename,vy_m1_xcut)
	dsetname = trim('/xcut/vz_m1')
	call HdfWriteReal2D_X(dsetname,filename,vz_m1_xcut)
	dsetname = trim('/xcut/temp_m1')
	call HdfWriteReal2D_X(dsetname,filename,temp_m1_xcut)
	dsetname = trim('/xcut/co2_m1')
	call HdfWriteReal2D_X(dsetname,filename,co2_m1_xcut)
	dsetname = trim('/xcut/h2o_m1')
	call HdfWriteReal2D_X(dsetname,filename,h2o_m1_xcut)

	dsetname = trim('/xcut/vx_m2')
	call HdfWriteReal2D_X(dsetname,filename,vx_m2_xcut)
	dsetname = trim('/xcut/vy_m2')
	call HdfWriteReal2D_X(dsetname,filename,vy_m2_xcut)
	dsetname = trim('/xcut/vz_m2')
	call HdfWriteReal2D_X(dsetname,filename,vz_m2_xcut)
	dsetname = trim('/xcut/temp_m2')
	call HdfWriteReal2D_X(dsetname,filename,temp_m2_xcut)
	dsetname = trim('/xcut/co2_m2')
	call HdfWriteReal2D_X(dsetname,filename,co2_m2_xcut)
	dsetname = trim('/xcut/h2o_m2')
	call HdfWriteReal2D_X(dsetname,filename,h2o_m2_xcut)

	dsetname = trim('/xcut/vx_m3')
	call HdfWriteReal2D_X(dsetname,filename,vx_m3_xcut)
	dsetname = trim('/xcut/vy_m3')
	call HdfWriteReal2D_X(dsetname,filename,vy_m3_xcut)
	dsetname = trim('/xcut/vz_m3')
	call HdfWriteReal2D_X(dsetname,filename,vz_m3_xcut)
	dsetname = trim('/xcut/temp_m3')
	call HdfWriteReal2D_X(dsetname,filename,temp_m3_xcut)
	dsetname = trim('/xcut/co2_m3')
	call HdfWriteReal2D_X(dsetname,filename,co2_m3_xcut)
	dsetname = trim('/xcut/h2o_m3')
	call HdfWriteReal2D_X(dsetname,filename,h2o_m3_xcut)

	dsetname = trim('/xcut/vx_m4')
	call HdfWriteReal2D_X(dsetname,filename,vx_m4_xcut)
	dsetname = trim('/xcut/vy_m4')
	call HdfWriteReal2D_X(dsetname,filename,vy_m4_xcut)
	dsetname = trim('/xcut/vz_m4')
	call HdfWriteReal2D_X(dsetname,filename,vz_m4_xcut)
	dsetname = trim('/xcut/temp_m4')
	call HdfWriteReal2D_X(dsetname,filename,temp_m4_xcut)
	dsetname = trim('/xcut/co2_m4')
	call HdfWriteReal2D_X(dsetname,filename,co2_m4_xcut)
	dsetname = trim('/xcut/h2o_m4')
	call HdfWriteReal2D_X(dsetname,filename,h2o_m4_xcut)

	! Y-CUT

	if ((yc(xstart(2)).le.stats2Dy).and.(yc(xend(2)+lvlhalo).gt.stats2Dy)) then

		linkname = trim('/ycut/')
		call HdfParallelCreateGroup(linkname,filename,comm_ycut)

		dsetname = trim('/ycut/vx_m1')
		call HdfWriteReal2D_Y(dsetname,filename,vx_m1_ycut)
		dsetname = trim('/ycut/vy_m1')
		call HdfWriteReal2D_Y(dsetname,filename,vy_m1_ycut)
		dsetname = trim('/ycut/vz_m1')
		call HdfWriteReal2D_Y(dsetname,filename,vz_m1_ycut)
		dsetname = trim('/ycut/temp_m1')
		call HdfWriteReal2D_Y(dsetname,filename,temp_m1_ycut)
		dsetname = trim('/ycut/co2_m1')
		call HdfWriteReal2D_Y(dsetname,filename,co2_m1_ycut)
		dsetname = trim('/ycut/h2o_m1')
		call HdfWriteReal2D_Y(dsetname,filename,h2o_m1_ycut)

		dsetname = trim('/ycut/vx_m2')
		call HdfWriteReal2D_Y(dsetname,filename,vx_m2_ycut)
		dsetname = trim('/ycut/vy_m2')
		call HdfWriteReal2D_Y(dsetname,filename,vy_m2_ycut)
		dsetname = trim('/ycut/vz_m2')
		call HdfWriteReal2D_Y(dsetname,filename,vz_m2_ycut)
		dsetname = trim('/ycut/temp_m2')
		call HdfWriteReal2D_Y(dsetname,filename,temp_m2_ycut)
		dsetname = trim('/ycut/co2_m2')
		call HdfWriteReal2D_Y(dsetname,filename,co2_m2_ycut)
		dsetname = trim('/ycut/h2o_m2')
		call HdfWriteReal2D_Y(dsetname,filename,h2o_m2_ycut)

		dsetname = trim('/ycut/vx_m3')
		call HdfWriteReal2D_Y(dsetname,filename,vx_m3_ycut)
		dsetname = trim('/ycut/vy_m3')
		call HdfWriteReal2D_Y(dsetname,filename,vy_m3_ycut)
		dsetname = trim('/ycut/vz_m3')
		call HdfWriteReal2D_Y(dsetname,filename,vz_m3_ycut)
		dsetname = trim('/ycut/temp_m3')
		call HdfWriteReal2D_Y(dsetname,filename,temp_m3_ycut)
		dsetname = trim('/ycut/co2_m3')
		call HdfWriteReal2D_Y(dsetname,filename,co2_m3_ycut)
		dsetname = trim('/ycut/h2o_m3')
		call HdfWriteReal2D_Y(dsetname,filename,h2o_m3_ycut)

		dsetname = trim('/ycut/vx_m4')
		call HdfWriteReal2D_Y(dsetname,filename,vx_m4_ycut)
		dsetname = trim('/ycut/vy_m4')
		call HdfWriteReal2D_Y(dsetname,filename,vy_m4_ycut)
		dsetname = trim('/ycut/vz_m4')
		call HdfWriteReal2D_Y(dsetname,filename,vz_m4_ycut)
		dsetname = trim('/ycut/temp_m4')
		call HdfWriteReal2D_Y(dsetname,filename,temp_m4_ycut)
		dsetname = trim('/ycut/co2_m4')
		call HdfWriteReal2D_Y(dsetname,filename,co2_m4_ycut)
		dsetname = trim('/ycut/h2o_m4')
		call HdfWriteReal2D_Y(dsetname,filename,h2o_m4_ycut)

	end if

	! Z-CUT

	if ((zc(xstart(3)).le.stats2Dz).and.(zc(xend(3)+lvlhalo).gt.stats2Dz)) then

		linkname = trim('/zcut')
		call HdfParallelCreateGroup(linkname,filename,comm_zcut)

		dsetname = trim('/zcut/vx_m1')
		call HdfWriteReal2D_Z(dsetname,filename,vx_m1_zcut)
		dsetname = trim('/zcut/vy_m1')
		call HdfWriteReal2D_Z(dsetname,filename,vy_m1_zcut)
		dsetname = trim('/zcut/vz_m1')
		call HdfWriteReal2D_Z(dsetname,filename,vz_m1_zcut)
		dsetname = trim('/zcut/temp_m1')
		call HdfWriteReal2D_Z(dsetname,filename,temp_m1_zcut)
		dsetname = trim('/zcut/co2_m1')
		call HdfWriteReal2D_Z(dsetname,filename,co2_m1_zcut)
		dsetname = trim('/zcut/h2o_m1')
		call HdfWriteReal2D_Z(dsetname,filename,h2o_m1_zcut)

		dsetname = trim('/zcut/vx_m2')
		call HdfWriteReal2D_Z(dsetname,filename,vx_m2_zcut)
		dsetname = trim('/zcut/vy_m2')
		call HdfWriteReal2D_Z(dsetname,filename,vy_m2_zcut)
		dsetname = trim('/zcut/vz_m2')
		call HdfWriteReal2D_Z(dsetname,filename,vz_m2_zcut)
		dsetname = trim('/zcut/temp_m2')
		call HdfWriteReal2D_Z(dsetname,filename,temp_m2_zcut)
		dsetname = trim('/zcut/co2_m2')
		call HdfWriteReal2D_Z(dsetname,filename,co2_m2_zcut)
		dsetname = trim('/zcut/h2o_m2')
		call HdfWriteReal2D_Z(dsetname,filename,h2o_m2_zcut)

		dsetname = trim('/zcut/vx_m3')
		call HdfWriteReal2D_Z(dsetname,filename,vx_m3_zcut)
		dsetname = trim('/zcut/vy_m3')
		call HdfWriteReal2D_Z(dsetname,filename,vy_m3_zcut)
		dsetname = trim('/zcut/vz_m3')
		call HdfWriteReal2D_Z(dsetname,filename,vz_m3_zcut)
		dsetname = trim('/zcut/temp_m3')
		call HdfWriteReal2D_Z(dsetname,filename,temp_m3_zcut)
		dsetname = trim('/zcut/co2_m3')
		call HdfWriteReal2D_Z(dsetname,filename,co2_m3_zcut)
		dsetname = trim('/zcut/h2o_m3')
		call HdfWriteReal2D_Z(dsetname,filename,h2o_m3_zcut)

		dsetname = trim('/zcut/vx_m4')
		call HdfWriteReal2D_Z(dsetname,filename,vx_m4_zcut)
		dsetname = trim('/zcut/vy_m4')
		call HdfWriteReal2D_Z(dsetname,filename,vy_m4_zcut)
		dsetname = trim('/zcut/vz_m4')
		call HdfWriteReal2D_Z(dsetname,filename,vz_m4_zcut)
		dsetname = trim('/zcut/temp_m4')
		call HdfWriteReal2D_Z(dsetname,filename,temp_m4_zcut)
		dsetname = trim('/zcut/co2_m4')
		call HdfWriteReal2D_Z(dsetname,filename,co2_m4_zcut)
		dsetname = trim('/zcut/h2o_m4')
		call HdfWriteReal2D_Z(dsetname,filename,h2o_m4_zcut)

	end if

	return

end subroutine WriteStats

subroutine WriteStatsEnd
    implicit none

    CHARACTER*50 :: fname
    
    fname = trim('Results/Stats/stafield_master.h5')
    call WriteStats(fname)
    
    return
end subroutine WriteStatsEnd

subroutine WriteStatsSnap
    use param, only: time
    
    implicit none
    
    CHARACTER*50 :: fname
    CHARACTER*8 :: citime
    
    write(citime,"(I8.8)") nint(time)
    fname = trim('Results/Stats/stafield_master'//trim(citime)//'.h5')
    call WriteStats(fname)
    
    return
end subroutine WriteStatsSnap