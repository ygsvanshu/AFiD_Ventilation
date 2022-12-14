!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: WriteFlowFieldSnapshot.F90                     !
!    CONTAINS: subroutine WriteFlowFieldSnapshot          !
!                                                         !
!    PURPOSE: Write down the full snapshot of temp,       !
!    vx,vy,vz for post processing purposes                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Modified on 05/02/2020 by Vanshu [Modification #14]
!! Modified on 30/12/2019 by Vanshu [Modification #1]

subroutine WriteFlowFieldSnapshot

    use param
    use local_arrays, only:vz,vy,vx,temp,co2,h2o

    implicit none

    character*50 :: filnam1,namfile
    ! ==============================================================================
    !  ModR05 Robert 2020-09-09
    !     Write Gridinfo snapshots
    ! ==============================================================================
    character*50 :: dsetname
    ! ==============================================================================
    !  End of ModR05
    ! ==============================================================================
    integer :: int_time
    character*8 :: citime

    int_time = nint(time)
    write(citime,"(I8.8)") int_time

    filnam1 = trim('Results/Vx/continua_vx')
    namfile=trim(trim(filnam1)//trim(citime)//'.h5')
    call HdfWriteRealHalo3D(namfile,vx)
    	
    filnam1 = trim('Results/Vy/continua_vy')
    namfile=trim(trim(filnam1)//trim(citime)//'.h5')
    call HdfWriteRealHalo3D(namfile,vy)
    
    filnam1 = trim('Results/Vz/continua_vz')
    namfile=trim(trim(filnam1)//trim(citime)//'.h5')
    call HdfWriteRealHalo3D(namfile,vz)

    filnam1 = trim('Results/Temp/continua_temp')
    namfile=trim(trim(filnam1)//trim(citime)//'.h5')
    call HdfWriteRealHalo3D(namfile,temp)

    filnam1 = trim('Results/H2O/continua_h2o')
    namfile=trim(trim(filnam1)//trim(citime)//'.h5')
    call HdfWriteRealHalo3D(namfile,co2)

    filnam1 = trim('Results/CO2/continua_co2')
    namfile=trim(trim(filnam1)//trim(citime)//'.h5')
    call HdfWriteRealHalo3D(namfile,h2o)

    ! ==============================================================================
    !  ModR05 Robert 2020-09-09
    !     Write Gridinfo snapshots
    ! ==============================================================================
    if (ismaster) then !EP only write once

		filnam1 = trim('Results/Grid/continua_master')
                namfile=trim(trim(filnam1)//trim(citime)//'.h5')
		call HdfCreateBlankFile(namfile)

		dsetname = trim('nx')
		call HdfSerialWriteIntScalar(dsetname,namfile,nx)
		dsetname = trim('ny')
		call HdfSerialWriteIntScalar(dsetname,namfile,ny)
		dsetname = trim('nz')
		call HdfSerialWriteIntScalar(dsetname,namfile,nz)
		dsetname = trim('ylen')
		call HdfSerialWriteRealScalar(dsetname,namfile,ylen)
		dsetname = trim('zlen')
		call HdfSerialWriteRealScalar(dsetname,namfile,zlen)
		dsetname = trim('time')
		call HdfSerialWriteRealScalar(dsetname,namfile,time)
		dsetname = trim('istr3')
		call HdfSerialWriteIntScalar(dsetname,namfile,istr3)
		dsetname = trim('str3')
		call HdfSerialWriteRealScalar(dsetname,namfile,str3)

	end if
    ! ==============================================================================
    !  End of ModR05
    ! ==============================================================================

end subroutine WriteFlowFieldSnapshot

!! End [Modification #1]
!! End [Modification #14]