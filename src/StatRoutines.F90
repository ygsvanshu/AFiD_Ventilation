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

	nstatsamples 	= 0
	tstat 			= 0.0d0
	tinterval 		= 0.0d0
	
	stat3d_vx_m1(:,:,:)   = 0.d0
	stat3d_vy_m1(:,:,:)   = 0.d0
	stat3d_vz_m1(:,:,:)   = 0.d0
	stat3d_pr_m1(:,:,:)   = 0.d0
	stat3d_temp_m1(:,:,:) = 0.d0
	stat3d_co2_m1(:,:,:)  = 0.d0
	stat3d_h2o_m1(:,:,:)  = 0.d0

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

	filename = trim('stafield_master.h5')
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

	else
		if(ismaster) write(6,*) 'Unable to read statistical files'
		if(ismaster) write(6,*) 'Restarting statistics from zero'
	end if

	filename = trim('stafield_vx_m1.h5')
	call HdfReadReal3D(filename,stat3d_vx_m1,xstart(2),xend(2),xstart(3),xend(3),nx,ny,nz)

	filename = trim('stafield_vy_m1.h5')
	call HdfReadReal3D(filename,stat3d_vy_m1,xstart(2),xend(2),xstart(3),xend(3),nx,ny,nz)

	filename = trim('stafield_vz_m1.h5')
	call HdfReadReal3D(filename,stat3d_vz_m1,xstart(2),xend(2),xstart(3),xend(3),nx,ny,nz)

	filename = trim('stafield_pr_m1.h5')
	call HdfReadReal3D(filename,stat3d_pr_m1,xstart(2),xend(2),xstart(3),xend(3),nx,ny,nz)

	filename = trim('stafield_temp_m1.h5')
	call HdfReadReal3D(filename,stat3d_temp_m1,xstart(2),xend(2),xstart(3),xend(3),nx,ny,nz)

	filename = trim('stafield_co2_m1.h5')
	call HdfReadReal3D(filename,stat3d_co2_m1,xstart(2),xend(2),xstart(3),xend(3),nx,ny,nz)

	filename = trim('stafield_h2o_m1.h5')
	call HdfReadReal3D(filename,stat3d_h2o_m1,xstart(2),xend(2),xstart(3),xend(3),nx,ny,nz)

end subroutine ReadStats

subroutine CalcStats

	use param
	use local_arrays, only: vx,vy,vz,pr,temp,co2,h2o
	use decomp_2d, only: xstart,xend
	use stat_arrays
	use mpih

    implicit none

	integer :: i,j,k

	nstatsamples = nstatsamples + 1
	tstat        = tstat        + dt

	! For X-CUT

	do k=1,nx
		do j=xstart(2),xend(2)
			do i=xstart(3),xend(3)
				stat3d_vx_m1(k,j,i)   = stat3d_vx_m1(k,j,i)   + dt*vx(k,j,i)
				stat3d_vy_m1(k,j,i)   = stat3d_vy_m1(k,j,i)   + dt*vy(k,j,i)
				stat3d_vz_m1(k,j,i)   = stat3d_vz_m1(k,j,i)   + dt*vz(k,j,i)
				stat3d_pr_m1(k,j,i)   = stat3d_pr_m1(k,j,i)   + dt*pr(k,j,i)
				stat3d_temp_m1(k,j,i) = stat3d_temp_m1(k,j,i) + dt*temp(k,j,i)
				stat3d_co2_m1(k,j,i)  = stat3d_co2_m1(k,j,i)  + dt*co2(k,j,i)
				stat3d_h2o_m1(k,j,i)  = stat3d_h2o_m1(k,j,i)  + dt*h2o(k,j,i)
			end do
		end do
	end do

	return

end subroutine CalcStats

subroutine WriteStatsSnap

	use mpih
	use param
	use decomp_2d, only: xstart,xend
	use stat_arrays
	use hdf5
	
	implicit none

	character*50				:: filename,linkname,dsetname
    character*8 				:: citime
    
    write(citime,"(I8.8)") nint(time)

	if (ismaster) then

		filename = trim('Results/Grid/stafield_master'//trim(citime)//'.h5')

		dsetname = trim('averaging_samples')
		call HdfSerialWriteIntScalar(dsetname,filename,nstatsamples)
		dsetname = trim('averaging_time')
		call HdfSerialWriteRealScalar(dsetname,filename,tstat)
		dsetname = trim('averaging_interval')
		call HdfSerialWriteRealScalar(dsetname,filename,tinterval+time-tsta)
		dsetname = trim('sample_interval')
		call HdfSerialWriteIntScalar(dsetname,filename,nout)

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
	
	filename = trim('Results/Vx/stafield_vx_m1'//trim(citime)//'.h5')
	call HdfWriteReal3D(filename,stat3d_vx_m1,xstart(2),xend(2),xstart(3),xend(3),nx,ny,nz)

	filename = trim('Results/Vy/stafield_vy_m1'//trim(citime)//'.h5')
	call HdfWriteReal3D(filename,stat3d_vy_m1,xstart(2),xend(2),xstart(3),xend(3),nx,ny,nz)

	filename = trim('Results/Vz/stafield_vz_m1'//trim(citime)//'.h5')
	call HdfWriteReal3D(filename,stat3d_vz_m1,xstart(2),xend(2),xstart(3),xend(3),nx,ny,nz)

	filename = trim('Results/Pr/stafield_pr_m1'//trim(citime)//'.h5')
	call HdfWriteReal3D(filename,stat3d_pr_m1,xstart(2),xend(2),xstart(3),xend(3),nx,ny,nz)

	filename = trim('Results/Temp/stafield_temp_m1'//trim(citime)//'.h5')
	call HdfWriteReal3D(filename,stat3d_temp_m1,xstart(2),xend(2),xstart(3),xend(3),nx,ny,nz)

	filename = trim('Results/CO2/stafield_co2_m1'//trim(citime)//'.h5')
	call HdfWriteReal3D(filename,stat3d_co2_m1,xstart(2),xend(2),xstart(3),xend(3),nx,ny,nz)

	filename = trim('Results/H2O/stafield_h2o_m1'//trim(citime)//'.h5')
	call HdfWriteReal3D(filename,stat3d_h2o_m1,xstart(2),xend(2),xstart(3),xend(3),nx,ny,nz)

	return

end subroutine WriteStatsSnap

subroutine WriteStatsEnd

	use mpih
	use param
	use decomp_2d, only: xstart,xend
	use stat_arrays
	use hdf5
	
	implicit none

	character*50				:: filename,linkname,dsetname
    
	if (ismaster) then

		filename = trim('stafield_master.h5')

		dsetname = trim('averaging_samples')
		call HdfSerialWriteIntScalar(dsetname,filename,nstatsamples)
		dsetname = trim('averaging_time')
		call HdfSerialWriteRealScalar(dsetname,filename,tstat)
		dsetname = trim('averaging_interval')
		call HdfSerialWriteRealScalar(dsetname,filename,tinterval+time-tsta)
		dsetname = trim('sample_interval')
		call HdfSerialWriteIntScalar(dsetname,filename,nout)

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
	
	filename = trim('stafield_vx_m1.h5')
	call HdfWriteReal3D(filename,stat3d_vx_m1,xstart(2),xend(2),xstart(3),xend(3),nx,ny,nz)

	filename = trim('stafield_vy_m1.h5')
	call HdfWriteReal3D(filename,stat3d_vy_m1,xstart(2),xend(2),xstart(3),xend(3),nx,ny,nz)

	filename = trim('stafield_vz_m1.h5')
	call HdfWriteReal3D(filename,stat3d_vz_m1,xstart(2),xend(2),xstart(3),xend(3),nx,ny,nz)

	filename = trim('stafield_pr_m1.h5')
	call HdfWriteReal3D(filename,stat3d_pr_m1,xstart(2),xend(2),xstart(3),xend(3),nx,ny,nz)

	filename = trim('stafield_temp_m1.h5')
	call HdfWriteReal3D(filename,stat3d_temp_m1,xstart(2),xend(2),xstart(3),xend(3),nx,ny,nz)

	filename = trim('stafield_co2_m1.h5')
	call HdfWriteReal3D(filename,stat3d_co2_m1,xstart(2),xend(2),xstart(3),xend(3),nx,ny,nz)

	filename = trim('stafield_h2o_m1.h5')
	call HdfWriteReal3D(filename,stat3d_h2o_m1,xstart(2),xend(2),xstart(3),xend(3),nx,ny,nz)

	return

end subroutine WriteStatsEnd