!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: WriteFlowField.F90                             !
!    CONTAINS: subroutine WriteFlowField                  !
!                                                         ! 
!    PURPOSE: Write down the full flow snapshot for       !
!     restarting the simulation at a later date           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Modified on 10/02/2020 by Vanshu [Modification #16]

subroutine WriteFlowField
	
	use param
	use local_arrays, only: vz,vy,vx,temp,co2,h2o
	
	implicit none
	
	character*50 :: filnam1,dsetname			!! <- character should be of size 50 because of edits to the HdfRoutines.F90

	filnam1 = trim('continua_vx.h5')
	call HdfWriteRealHalo3D(filnam1,vx)
	filnam1 = trim('continua_vy.h5')
	call HdfWriteRealHalo3D(filnam1,vy)
	filnam1 = trim('continua_vz.h5')
	call HdfWriteRealHalo3D(filnam1,vz)
	filnam1 = trim('continua_temp.h5')
	call HdfWriteRealHalo3D(filnam1,temp)
	filnam1 = trim('continua_co2.h5')
	call HdfWriteRealHalo3D(filnam1,co2)
	filnam1 = trim('continua_h2o.h5')
	call HdfWriteRealHalo3D(filnam1,h2o)

	if (ismaster) then !EP only write once

		filnam1 = trim('continua_master.h5')
		call HdfCreateBlankFile(filnam1)

		dsetname = trim('nx')
		call HdfSerialWriteIntScalar(dsetname,filnam1,nx)
		dsetname = trim('ny')
		call HdfSerialWriteIntScalar(dsetname,filnam1,ny)
		dsetname = trim('nz')
		call HdfSerialWriteIntScalar(dsetname,filnam1,nz)
		dsetname = trim('ylen')
		call HdfSerialWriteRealScalar(dsetname,filnam1,ylen)
		dsetname = trim('zlen')
		call HdfSerialWriteRealScalar(dsetname,filnam1,zlen)
		dsetname = trim('time')
		call HdfSerialWriteRealScalar(dsetname,filnam1,time)
		dsetname = trim('istr3')
		call HdfSerialWriteIntScalar(dsetname,filnam1,istr3)
		dsetname = trim('str3')
		call HdfSerialWriteRealScalar(dsetname,filnam1,str3)

	end if

end subroutine WriteFlowField

!! End [Modification #16]