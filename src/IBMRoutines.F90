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

    use mpih
    use param
    use decomp_2d, only: xstart,xend,ystart,yend,zstart,zend,transpose_x_to_y,transpose_y_to_z
    use local_arrays, only: rhsx,rhsy,rhsz
    use ibm_arrays

    implicit none

    integer         :: kc,kp,jc,jp,ic,ip
    character*50    :: filnam
    real            :: chksum,res_dummy,val_max,val_min,val_max_all,val_min_all
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

    ! Initialize the ibm body 
    rhsx(:,:,:) = 0.0d0
    rhsy(:,:,:) = 0.0d0
    rhsz(:,:,:) = 0.0d0
    chksum = 0
    ! Check for points inside the geometry and set the ibm body to 1.0
    do ic = xstart(3),xend(3)
        do jc = xstart(2),xend(2)
            do kc = 1,nxm
                ! This is such that the person faces the -z direction i.e. the inlet vent
                call polyhedron_contains_point_3d ( node_num, face_num, order_max, node_xyz, face_order, face_node, (/ym(jc),zm(ic),xc(kc)/), ibm_gx_px(kc,jc,ic))
                call polyhedron_contains_point_3d ( node_num, face_num, order_max, node_xyz, face_order, face_node, (/yc(jc),zm(ic),xm(kc)/), ibm_gy_px(kc,jc,ic))
                call polyhedron_contains_point_3d ( node_num, face_num, order_max, node_xyz, face_order, face_node, (/ym(jc),zc(ic),xm(kc)/), ibm_gz_px(kc,jc,ic))
            enddo
        enddo
    enddo

    do ic = xstart(3),xend(3)
        ip = min(ic+1,xend(3))
        do jc = xstart(2),xend(2)
            jp = min(jc+1,xend(2))
            do kc = 1,nxm
                kp = kc+1
                ! Check if all face velocity components are non inside the IBM body
                if(ibm_gx_px(kc,jc,ic).or.ibm_gy_px(kc,jc,ic).or.ibm_gz_px(kc,jc,ic).or.ibm_gx_px(kp,jc,ic).or.ibm_gy_px(kc,jp,ic).or.ibm_gz_px(kc,jc,ip)) then
                    ibm_gc_px(kc,jc,ic) = .true.
                    rhsx(kc,jc,ic) = 1.0d0
                    chksum = chksum + 1.0d0
                endif
            enddo
        enddo
    enddo

    call MpiSumRealScalar(chksum,res_dummy)
    chksum = res_dummy
    if(ismaster) write(6,*) 'Number of cells intersecting IBM person geometry = ',int(chksum)

    filnam = trim('ibm_body.h5')
    call HdfWriteRealHalo3D(filnam,rhsx)

    ! Transforms for the x-grid
    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                if (ibm_gx_px(kc,jc,ic)) rhsx = 1.0d0
            enddo
        enddo
    enddo
    call transpose_x_to_y(rhsx,rhsy)
    do ic=ystart(3),yend(3)
        do kc=ystart(1),yend(1)
            do jc=1,nym
                if (rhsy(kc,jc,ic).gt.0.5d0) ibm_gx_py = .true.
            enddo
        enddo
    enddo
    call transpose_y_to_z(rhsy,rhsz)
    do jc=zstart(2),zend(2)
        do kc=zstart(1),zend(1)
            do ic=1,nzm
                if (rhsz(kc,jc,ic).gt.0.5d0) ibm_gx_pz = .true.
            enddo
        enddo
    enddo
    rhsx(nx,:,:) = 0.0d0

    ! Transforms for the y-grid
    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                if (ibm_gy_px(kc,jc,ic)) rhsx = 1.0d0
            enddo
        enddo
    enddo
    call transpose_x_to_y(rhsx,rhsy)
    do ic=ystart(3),yend(3)
        do kc=ystart(1),yend(1)
            do jc=1,nym
                if (rhsy(kc,jc,ic).gt.0.5d0) ibm_gy_py = .true.
            enddo
        enddo
    enddo
    call transpose_y_to_z(rhsy,rhsz)
    do jc=zstart(2),zend(2)
        do kc=zstart(1),zend(1)
            do ic=1,nzm
                if (rhsz(kc,jc,ic).gt.0.5d0) ibm_gy_pz = .true.
            enddo
        enddo
    enddo

    ! Transforms for the z-grid
    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                if (ibm_gz_px(kc,jc,ic)) rhsx = 1.0d0
            enddo
        enddo
    enddo
    call transpose_x_to_y(rhsx,rhsy)
    do ic=ystart(3),yend(3)
        do kc=ystart(1),yend(1)
            do jc=1,nym
                if (rhsy(kc,jc,ic).gt.0.5d0) ibm_gz_py = .true.
            enddo
        enddo
    enddo
    call transpose_y_to_z(rhsy,rhsz)
    do jc=zstart(2),zend(2)
        do kc=zstart(1),zend(1)
            do ic=1,nzm
                if (rhsz(kc,jc,ic).gt.0.5d0) ibm_gz_pz = .true.
            enddo
        enddo
    enddo

    deallocate ( face_node )
    deallocate ( face_order )
    deallocate ( node_xyz )
    deallocate ( normal_vector )
    deallocate ( vertex_normal )
    
    return

end subroutine CreateBodyIBM

subroutine AddBodyIBM(qua,ibm,val)

    use mpih
    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: rhsx
    
    implicit none

    real,    intent(in), dimension(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo) :: qua
    logical, intent(in), dimension(1:nx,xstart(2):xend(2),xstart(3):xend(3)) :: ibm
    real,    intent(in) :: val
    integer :: kc,jc,ic

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nx
                if (ibm(kc,jc,ic)) rhsx(kc,jc,ic) = val - qua(kc,jc,ic)
            end do
        end do
    end do

end subroutine AddBodyIBM

subroutine AddBreathIBM

    use param
    use decomp_2d, only: xstart,xend
    ! use local_arrays, only: vx,vy,vz,temp,co2,h2o
    use local_arrays, only: qcap,dph,dq,hro,qco2,qh2o
    use mpih

    implicit none

    integer :: kc,jc,ic
    real    :: time_shift,breath_interval,time_signal,vel_peak,time_signal_exp,space_signal  
    real    :: tprefactor,sprefactor,qprefactor
    real    :: injectedvol,injectmeanvx,injectmeanvy,injectmeanvz,injectmeantemp,injectmeanco2,injectmeanh2o

    ! Check if the breath forcing location is out of the domain volume
    if ((breathx.lt.0.0d0).or.(breathx.gt.alx3).or.(breathy.lt.0.0d0).or.(breathy.gt.ylen).or.(breathz.lt.0.0d0).or.(breathz.gt.zlen)) then
        call NotifyError(667)
        call MPI_Abort(MPI_COMM_WORLD,1) 
    end if

    ! Compute injection quantities
    injectedvol     = 5e-4  /3.0/3.0/3.0       ! normalized 0.5L (by length scale 3m)
    injectmeanvx    = -0.5*dsin(pi/3.0)/0.71   ! normalized 0.5m/s with angle (by free fall vel 0.71m/s)
    injectmeanvy    = 0.0
    injectmeanvz    = -0.5*dcos(pi/3.0)/0.71   ! normalized 0.5m/s with angle (by free fall vel 0.71m/s)
    injectmeantemp  = 1.0
    injectmeanco2   = 1.0
    injectmeanh2o   = 1.0

    ! Set temporal Gaussian func
    time_shift      = 2.0/4.25             ! normalized 2s (by free fall time 4.25s)
    breath_interval = 4.25/4.25            ! normalized 4.25s (by free fall time 4.25s)

    time_signal_exp = exp(-0.5*( 2.0*(modulo(time,breath_interval)-time_shift)/kernel_width_time)**2)
    if(time_signal_exp.le.1d-8) time_signal_exp=0.d0
    tprefactor      = (2.0/(2.0*pi)**0.5)/kernel_width_time
    time_signal     = tprefactor*time_signal_exp
    sprefactor      = (2.0/(2.0*pi)**0.5)**3.0/kernel_width_space/kernel_width_space/kernel_width_space

    ! Set breath
    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nx     
                space_signal    = exp(-0.5*((2.0*(xc(kc)-breathx)/kernel_width_space)**2  + (2.0*(ym(jc)-breathy)/kernel_width_space)**2 + (2.0*(zm(ic)-breathz)/kernel_width_space)**2))
                qprefactor      = (sprefactor*space_signal*time_signal*injectedvol)!*ga*dt)
                qcap(kc,jc,ic)  = qcap(kc,jc,ic) + (injectmeanvx*qprefactor)
                hro(kc,jc,ic)   = hro(kc,jc,ic)  + (injectmeantemp*qprefactor)
                qco2(kc,jc,ic)  = qco2(kc,jc,ic) + (injectmeanco2*qprefactor)
                qh2o(kc,jc,ic)  = qh2o(kc,jc,ic) + (injectmeanh2o*qprefactor)

            end do
            do kc=1,nxm
                space_signal    = exp(-0.5*((2.0*(xm(kc)-breathx)/kernel_width_space)**2  + (2.0*(yc(jc)-breathy)/kernel_width_space)**2 + (2.0*(zm(ic)-breathz)/kernel_width_space)**2))
                qprefactor      = (sprefactor*space_signal*time_signal*injectedvol)!*ga*dt)
                dph(kc,jc,ic)   = dph(kc,jc,ic)  + (injectmeanvy*qprefactor)
                space_signal    = exp(-0.5*((2.0*(xm(kc)-breathx)/kernel_width_space)**2  + (2.0*(ym(jc)-breathy)/kernel_width_space)**2 + (2.0*(zc(ic)-breathz)/kernel_width_space)**2))
                qprefactor      = (sprefactor*space_signal*time_signal*injectedvol)!*ga*dt)
                dq(kc,jc,ic)    = dq(kc,jc,ic)   + (injectmeanvz*qprefactor)
            enddo
        enddo
    enddo

    return

end subroutine AddBreathIBM

! ********** FOR DEBUG PURPOSES ********** !
! Use a sphere as a test IBM body
subroutine CreateDebugBodyIBM

    use mpih
    use param
    use decomp_2d, only: xstart,xend,ystart,yend,zstart,zend,transpose_x_to_y,transpose_y_to_z
    use local_arrays, only: rhsx,rhsy,rhsz
    use ibm_arrays

    implicit none

    integer         :: kc,kp,jc,jp,ic,ip
    character*50    :: filnam
    logical         :: inside_xc,inside_xp,inside_yc,inside_yp,inside_zc,inside_zp
    real            :: sqradius,sqdistance,chksum,res_dummy

    rhsx(:,:,:) = 0.0d0
    rhsy(:,:,:) = 0.0d0
    rhsz(:,:,:) = 0.0d0

    sqradius = 0.1d0**2.0d0

    do ic=xstart(3),xend(3)
        ip = ic + 1
        do jc=xstart(2),xend(2)
            jp = jc + 1
            do kc=1,nxm
                kp = kc + 1

                sqdistance = ((xc(kc)-personx)**2) + ((ym(jc)-persony)**2) + ((zm(ic)-personz)**2)
                if (sqdistance.le.sqradius) inside_xc = .true.
                sqdistance = ((xc(kp)-personx)**2) + ((ym(jc)-persony)**2) + ((zm(ic)-personz)**2)
                if (sqdistance.le.sqradius) inside_xp = .true.
                sqdistance = ((xm(kc)-personx)**2) + ((yc(jc)-persony)**2) + ((zm(ic)-personz)**2)
                if (sqdistance.le.sqradius) inside_yc = .true.
                sqdistance = ((xm(kc)-personx)**2) + ((yc(jp)-persony)**2) + ((zm(ic)-personz)**2)
                if (sqdistance.le.sqradius) inside_yp = .true.
                sqdistance = ((xm(kc)-personx)**2) + ((ym(jc)-persony)**2) + ((zc(ic)-personz)**2)
                if (sqdistance.le.sqradius) inside_zc = .true.
                sqdistance = ((xm(kc)-personx)**2) + ((ym(jc)-persony)**2) + ((zc(ip)-personz)**2)
                if (sqdistance.le.sqradius) inside_zp = .true.
                ! Check if all face velocity components are non inside the IBM body
                if(inside_xc.and.inside_xp.and.inside_yc.and.inside_yp.and.inside_zc.and.inside_zp) then
                    ibm_gc_px(kc,jc,ic) = .true.
                    rhsx(kc,jc,ic) = 1.0d0
                    chksum = chksum + 1.0d0
                endif
                ibm_gx_px(kc,jc,ic) = inside_xc
                ibm_gy_px(kc,jc,ic) = inside_yc
                ibm_gz_px(kc,jc,ic) = inside_zc
            enddo
        enddo
    enddo

    call MpiSumRealScalar(chksum,res_dummy)
    chksum = res_dummy
    if(ismaster) write(6,*) 'Number of cells intersecting IBM person geometry = ',int(chksum)

    filnam = trim('ibm_body.h5')
    call HdfWriteRealHalo3D(filnam,rhsx)

    ! Transforms for the x-grid
    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nx
                if (ibm_gx_px(kc,jc,ic)) rhsx = 1.0d0
            enddo
        enddo
    enddo
    call transpose_x_to_y(rhsx,rhsy)
    do ic=ystart(3),yend(3)
        do kc=ystart(1),yend(1)
            do jc=1,nym
                if (rhsy(kc,jc,ic).gt.0.5d0) ibm_gx_py = .true.
            enddo
        enddo
    enddo
    call transpose_y_to_z(rhsy,rhsz)
    do jc=zstart(2),zend(2)
        do kc=zstart(1),zend(1)
            do ic=1,nzm
                if (rhsz(kc,jc,ic).gt.0.5d0) ibm_gx_pz = .true.
            enddo
        enddo
    enddo
    rhsx(nx,:,:) = 0.0d0

    ! Transforms for the y-grid
    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                if (ibm_gy_px(kc,jc,ic)) rhsx = 1.0d0
            enddo
        enddo
    enddo
    call transpose_x_to_y(rhsx,rhsy)
    do ic=ystart(3),yend(3)
        do kc=ystart(1),yend(1)
            do jc=1,ny
                if (rhsy(kc,jc,ic).gt.0.5d0) ibm_gy_py = .true.
            enddo
        enddo
    enddo
    call transpose_y_to_z(rhsy,rhsz)
    do jc=zstart(2),zend(2)
        do kc=zstart(1),zend(1)
            do ic=1,nzm
                if (rhsz(kc,jc,ic).gt.0.5d0) ibm_gy_pz = .true.
            enddo
        enddo
    enddo

    ! Transforms for the z-grid
    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                if (ibm_gz_px(kc,jc,ic)) rhsx = 1.0d0
            enddo
        enddo
    enddo
    call transpose_x_to_y(rhsx,rhsy)
    do ic=ystart(3),yend(3)
        do kc=ystart(1),yend(1)
            do jc=1,nym
                if (rhsy(kc,jc,ic).gt.0.5d0) ibm_gz_py = .true.
            enddo
        enddo
    enddo
    call transpose_y_to_z(rhsy,rhsz)
    do jc=zstart(2),zend(2)
        do kc=zstart(1),zend(1)
            do ic=1,nz
                if (rhsz(kc,jc,ic).gt.0.5d0) ibm_gz_pz = .true.
            enddo
        enddo
    enddo

    return

end subroutine CreateDebugBodyIBM
! ********** FOR DEBUG PURPOSES ********** !
