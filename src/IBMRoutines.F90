!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                           ! 
!    FILE: IBMRoutines.F90                                  !
!    CONTAINS: subroutines CreateBodyIBM, AddBodyIBM,       !
!              CreateBreathIBM and AddBreathIBM             !
!                                                           ! 
!    PURPOSE: IBM subroutines to add the person geometry    !
!                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CreateBodyIBM

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: ibm_body
    use mpih

    implicit none

    integer         :: kc,jc,ic
    character*50    :: filnam
    logical         :: inside
    real            :: p(3),chksum,res_dummy,val_max,val_min,val_max_all,val_min_all
    real            :: xmax,ymax,zmax, xmin,ymin,zmin, xdiff,ydiff,zdiff, alldiff

    integer ( kind = 4 ), allocatable, dimension(:,:) :: face_node
    integer ( kind = 4 ) face_num
    integer ( kind = 4 ), allocatable, dimension(:) :: face_order
    integer ( kind = 4 ) ierror
    integer ( kind = 4 ) node_num
    real    ( kind = 8 ), allocatable, dimension(:,:) :: node_xyz
    real    ( kind = 8 ), allocatable, dimension(:,:) :: normal_vector
    integer ( kind = 4 ) normal_num
    integer ( kind = 4 ) order_max
    integer ( kind = 4 ), allocatable, dimension(:,:) :: vertex_normal

    call obj_size ( person_objfile, node_num, face_num, normal_num, order_max )

    allocate ( face_node(order_max,face_num) )
    allocate ( face_order(face_num) )
    allocate ( node_xyz(3,node_num) )
    allocate ( normal_vector(3,normal_num) )
    allocate ( vertex_normal(order_max,face_num) )

    call obj_read ( person_objfile, node_num, face_num, normal_num, order_max, node_xyz, face_order, face_node, normal_vector, vertex_normal )

    ! Get the dimensions of the bounding box of the person geometry
    xmax=maxval(node_xyz(1,:))
    xmin=minval(node_xyz(1,:))
    ymax=maxval(node_xyz(2,:))
    ymin=minval(node_xyz(2,:))
    zmax=maxval(node_xyz(3,:))
    zmin=minval(node_xyz(3,:))
    
    xdiff=xmax-xmin
    ydiff=ymax-ymin
    zdiff=zmax-zmin

    ! Find the largest dimension of the bounding box for rescaling
    if(xdiff.ge.ydiff .and. xdiff.ge.zdiff) then 
        alldiff=xdiff
    elseif(ydiff.ge.xdiff .and. ydiff.ge.zdiff) then
        alldiff=ydiff
    else
        alldiff=zdiff
    endif

    ! Scale down such that the largest dimension of bounding box is unity
    node_xyz(1,:)=(node_xyz(1,:)-(xmax+xmin)/2.d0)/alldiff
    node_xyz(2,:)=(node_xyz(2,:)-(ymax+ymin)/2.d0)/alldiff
    node_xyz(3,:)=(node_xyz(3,:)-(zmax+zmin)/2.d0)/alldiff

    ! Scale up the largest dimension of bounding box by the scale factor "sclf"
    node_xyz(1,:)=node_xyz(1,:)*sclf + personz
    node_xyz(2,:)=node_xyz(2,:)*sclf + persony
    node_xyz(3,:)=node_xyz(3,:)*sclf + personx !+sclf/2.d0 !alx3/4.d0

    ! Check if any part of the person is jutting out of the domain
    xmax=maxval(node_xyz(1,:))
    xmin=minval(node_xyz(1,:))
    ymax=maxval(node_xyz(2,:))
    ymin=minval(node_xyz(2,:))
    zmax=maxval(node_xyz(3,:))
    zmin=minval(node_xyz(3,:))
    ! If true, then abort
    if ((xmin.lt.0.0d0).or.(xmax.gt.alx3).or.(ymin.lt.0.0d0).or.(ymax.gt.ylen).or.(zmin.lt.0.0d0).or.(zmax.gt.zlen)) then
        call NotifyError(666)
        call MPI_Abort(MPI_COMM_WORLD,1) 
    end if

    ! Some debugging routines that can be un-commented
    ! call obj_face_node_print ( face_num, order_max, face_order, face_node )
    ! call obj_normal_vector_print ( normal_num, normal_vector )
    ! call obj_node_xyz_print ( node_num, node_xyz )

    ! Initialize the ibm_body 
    ibm_body(:,:,:) = 0.d0
    chksum = 0
    ! Check for points inside the geometry and set the ibm_body to 1.0
    do kc = 1,nxm
        do jc = xstart(2),xend(2)
            do ic = xstart(3),xend(3)
                p = (/ym(jc),zm(ic),xm(kc)/)    ! This is such that the person faces the -z direction i.e. the inlet vent
                call polyhedron_contains_point_3d ( node_num, face_num, order_max, node_xyz, face_order, face_node, p, inside )
                if(inside) then
                    ibm_body(kc,jc,ic) = 1.d0
                    chksum = chksum + 1.0d0
                endif
            enddo
        enddo
    enddo

    call MpiSumRealScalar(chksum,res_dummy)
    chksum = res_dummy
    if(ismaster) write(6,*) 'Number of grid points in IBM person geometry = ',int(chksum)

    deallocate ( face_node )
    deallocate ( face_order )
    deallocate ( node_xyz )
    deallocate ( normal_vector )
    deallocate ( vertex_normal )
    
    filnam = trim('ibm_body.h5')
    call HdfWriteRealHalo3D(filnam,ibm_body)

    return

end subroutine CreateBodyIBM

subroutine AddBreathIBM

    use param
    use decomp_2d, only: xstart,xend
    ! use local_arrays, only: vx,vy,vz,temp,co2,h2o
    use local_arrays, only: qcap,dph,dq,hro,qco2,qh2o
    use mpih

    implicit none

    integer :: kc,jc,ic
    real    :: time_signal,space_signal  
    real    :: tprefactor,sprefactor,qprefactor
    real    :: injectedvol,injectmeanvx,injectmeanvy,injectmeanvz,injectmeantemp,injectmeanco2,injectmeanh2o

    ! Check if the breath forcing location is out of the domain volume
    if ((breathx.lt.0.0d0).or.(breathx.gt.alx3).or.(breathy.lt.0.0d0).or.(breathy.gt.ylen).or.(breathz.lt.0.0d0).or.(breathz.gt.zlen)) then
        call NotifyError(667)
        call MPI_Abort(MPI_COMM_WORLD,1) 
    end if

    ! Compute injection quantities
    ! injectedvol     = 5e-4  /3.0/3.0/3.0       ! normalized 0.5L (by length scale 3m)
    injectedvol     = breath_volume
    ! injectmeanvx    = -0.5*dsin(pi/3.0)/0.71   ! normalized 0.5m/s with angle (by free fall vel 0.71m/s)
    injectmeanvx    = breath_velocity*dsin(breath_angle*pi/180.0d0)
    injectmeanvy    = 0.0
    ! injectmeanvz    = -0.5*dcos(pi/3.0)/0.71   ! normalized 0.5m/s with angle (by free fall vel 0.71m/s)
    injectmeanvz    = breath_velocity*dcos(breath_angle*pi/180.0d0)
    injectmeantemp  = 1.0
    injectmeanco2   = 1.0
    injectmeanh2o   = 1.0

    ! Set temporal Gaussian func
    ! breath_offset   = 2.0/4.25             ! normalized 2s (by free fall time 4.25s)
    ! breath_interval = 4.25/4.25            ! normalized 4.25s (by free fall time 4.25s)

    tprefactor      = (2.0/(2.0*pi)**0.5)/kernel_width_time
    time_signal     = tprefactor*exp(-0.5*( 2.0*(modulo(time,breath_interval)-breath_offset)/kernel_width_time)**2)
    sprefactor      = (2.0/(2.0*pi)**0.5)**3.0/kernel_width_space/kernel_width_space/kernel_width_space

    ! Set breath
    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nx     
                space_signal    = exp(-0.5*((2.0*(xc(kc)-breathx)/kernel_width_space)**2  + (2.0*(ym(jc)-breathy)/kernel_width_space)**2 + (2.0*(zm(ic)-breathz)/kernel_width_space)**2))
                qprefactor      = (sprefactor*space_signal*time_signal*injectedvol)!*ga*dt)
                ! vx(kc,jc,ic)    = vx(kc,jc,ic)   + injectmeanvx*qprefactor
                ! temp(kc,jc,ic)  = temp(kc,jc,ic) + injectmeantemp*qprefactor
                ! co2(kc,jc,ic)   = co2(kc,jc,ic)  + injectmeanco2*qprefactor
                ! h2o(kc,jc,ic)   = h2o(kc,jc,ic)  + injectmeanh2o*qprefactor
                qcap(kc,jc,ic)  = qcap(kc,jc,ic) + (injectmeanvx*qprefactor)
                hro(kc,jc,ic)   = hro(kc,jc,ic)  + (injectmeantemp*qprefactor)
                qco2(kc,jc,ic)  = qco2(kc,jc,ic) + (injectmeanco2*qprefactor)
                qh2o(kc,jc,ic)  = qh2o(kc,jc,ic) + (injectmeanh2o*qprefactor)

            end do
            do kc=1,nxm
                space_signal    = exp(-0.5*((2.0*(xm(kc)-breathx)/kernel_width_space)**2  + (2.0*(yc(jc)-breathy)/kernel_width_space)**2 + (2.0*(zm(ic)-breathz)/kernel_width_space)**2))
                qprefactor      = (sprefactor*space_signal*time_signal*injectedvol)!*ga*dt)
                ! vy(kc,jc,ic)    = vy(kc,jc,ic)   + injectmeanvy*qprefactor
                dph(kc,jc,ic)   = dph(kc,jc,ic)  + (injectmeanvy*qprefactor)
                space_signal    = exp(-0.5*((2.0*(xm(kc)-breathx)/kernel_width_space)**2  + (2.0*(ym(jc)-breathy)/kernel_width_space)**2 + (2.0*(zc(ic)-breathz)/kernel_width_space)**2))
                qprefactor      = (sprefactor*space_signal*time_signal*injectedvol)!*ga*dt)
                ! vz(kc,jc,ic)    = vz(kc,jc,ic)   + injectmeanvz*qprefactor
                dq(kc,jc,ic)    = dq(kc,jc,ic)   + (injectmeanvz*qprefactor)
            enddo
        enddo
    enddo

    return

end subroutine AddBreathIBM

! ********** FOR DEBUG PURPOSES ********** !
! Use a sphere as a test IBM body
subroutine CreateDebugBodyIBM

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: ibm_body
    use mpih

    implicit none

    integer         :: kc,jc,ic
    character*50    :: filnam

    real    :: radius, distance
    radius = 0.1d0
    do kc=1,nxm
        do jc=xstart(2),xend(2)
            do ic=xstart(3),xend(3)
                distance = 0.0d0
                distance = distance + ((xm(kc)-personx)**2)
                distance = distance + ((ym(jc)-persony)**2)
                distance = distance + ((zm(ic)-personz)**2)
                distance = distance**0.5
                if (distance.le.radius) ibm_body(kc,jc,ic) = 1.0d0
            end do
        end do
    end do
    

    filnam = trim('ibm_body.h5')
    call HdfWriteRealHalo3D(filnam,ibm_body)

    return

end subroutine CreateDebugBodyIBM
! ********** FOR DEBUG PURPOSES ********** !
