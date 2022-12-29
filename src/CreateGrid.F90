!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CreateGrid.F90                                 !
!    CONTAINS: subroutine CreateGrid                      !
!                                                         ! 
!    PURPOSE: Compute the indices, grid, grid metrics     !
!     and coefficients for differentiation                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CreateGrid

    use param
    use AuxiliaryRoutines

    implicit none

    real :: x1,x2,x3
    real :: a33, a33m, a33p
    real :: delet, etain, tstr3
    real :: z2dp

    real :: alpha
    real :: x_temporary
    logical :: fexist

    real, allocatable, dimension(:) :: etaz, etazm

    integer :: i, j, kc, km, kp
    integer :: nxmo, nclip

    do kc=1,nxm
        kmv(kc)=kc-1
        kpv(kc)=kc+1
        if(kc.eq.1) kmv(kc)=kc
        if(kc.eq.nxm) kpv(kc)=kc
    end do

    do kc=1,nxm
        kpc(kc)=kpv(kc)-kc
        kmc(kc)=kc-kmv(kc)
    end do

    ! UNIFORM (HORIZONTAL DIRECTIONS) GRID
    
    do  i=1,nz
        x1=real(i-1)/real(nzm)
        zc(i)= zlen*x1
    end do

    do i=1,nzm
        zm(i)=(zc(i)+zc(i+1))*0.5d0
    end do

    do j=1,ny
        x2=real(j-1)/real(nym)
        yc(j)= ylen*x2
    end do

    do j=1,nym
        ym(j)=(yc(j)+yc(j+1))*0.5d0
    end do

    ! VERTICAL COORDINATE DEFINITION

    ! OPTION 0: UNIFORM CLUSTERING
    
    call AllocateReal1DArray(etaz,1,nx+500)
    call AllocateReal1DArray(etazm,1,nx+500)

    if (istr3.eq.0) then
        do kc=1,nx
            x3=real(kc-1)/real(nxm)
            etaz(kc)=alx3*x3
            xc(kc)=etaz(kc)
        enddo
    end if

    ! OPTION 4: HYPERBOLIC TANGENT-TYPE CLUSTERING

    tstr3=tanh(str3)

    if (istr3.eq.4) then
        xc(1)=0.0d0
        do kc=2,nx
            z2dp=float(2*kc-nx-1)/float(nxm)
            xc(kc)=(1+tanh(str3*z2dp)/tstr3)*0.5*alx3
            if(xc(kc).lt.0.or.xc(kc).gt.alx3)then
                write(*,*)'Grid is too streched: ','zc(',kc,')=',xc(kc)
                stop
            endif
        end do
    end if

    ! OPTION 6: CLIPPED CHEBYCHEV-TYPE CLUSTERING

    if(istr3.eq.6) then
        nclip = int(str3)
        nxmo = nx+nclip+nclip
        do kc=1,nxmo
            etazm(kc)=+cos(pi*(float(kc)-0.5)/float(nxmo))
        end do
        do kc=1,nx
            etaz(kc)=etazm(kc+nclip)
        end do
        delet = etaz(1)-etaz(nx)
        etain = etaz(1)
        do kc=1,nx
            etaz(kc)=etaz(kc)/(0.5*delet)
        end do
        xc(1) = 0.
        do kc=2,nxm
            xc(kc) = alx3*(1.-etaz(kc))*0.5
        end do
        xc(nx) = alx3
    end if

    call DestroyReal1DArray(etaz)
    call DestroyReal1DArray(etazm)

    ! OPTION 7: ERROR FUNCTION FOR VALIDATION OF PURE COUETTE FLOW WITH PIROZZOLI.S. et.al. "TURBULENCE STATISTICS IN COUETTE FLOW AT HIGH REYNOLDS NUMBERS" (2014)

    alpha=str3
	if(istr3.eq.7) then
		do kc=1,nx
			x_temporary	= (real(kc-1)/real(nxm)) - 0.5
			x_temporary = erf(alpha*x_temporary)/erf(alpha/2)
			xc(kc) 		= 0.5*(1+x_temporary)*alx3
		end do   
	end if	

    ! OPTION 8: READ GRID FROM GRID INPUT FILE

	if(istr3.eq.8) then
        inquire(file='./wall_normal_grid.in', exist=fexist) 
        if(fexist) then
            if (ismaster) write(6,*) 'Reading custom grid from wall_normal_grid.in'
        else
            if (ismaster) write(6,*) 'Warning: wall_normal_grid.in not found!'
            call MpiAbort
        end if
        open(unit=78,file='wall_normal_grid.in',status='old')
        do kc=1,nx
            read(78,*) xc(kc)
        end do
        close(78)
	endif
      
    ! METRIC FOR UNIFORM DIRECTIONS

    dx=real(nxm)/alx3
    dy=real(nym)/ylen
    dz=real(nzm)/zlen

    dxq=dx*dx                                                      
    dyq=dy*dy                                                      
    dzq=dz*dz                                                      

    ! STAGGERED COORDINATES AND
    ! METRIC QUANTITIES FOR NON-UNIFORM 
    ! DIRECTIONS

    do kc=1,nxm
        xm(kc)=(xc(kc)+xc(kc+1))*0.5d0
        dx3c(kc)=(xc(kc+1)-xc(kc))*dx
    enddo
      
    do kc=1,nxm-1
        dx3m(kc)=(xm(kc+1)-xm(kc))*dx
    enddo
      
    do kc=2,nxm
        g3rc(kc)=(xc(kc+1)-xc(kc-1))*dx*0.5d0
    enddo
    g3rc(1)=(xc(2)-xc(1))*dx
    g3rc(nx)= (xc(nx)-xc(nxm))*dx
      
    do kc=2,nxm-1
        g3rm(kc)=(xm(kc+1)-xm(kc-1))*dx*0.5d0
    enddo
    g3rm(1)   = (xm(2)  - xc(1))    *dx*0.5d0
    g3rm(nxm) = (xc(nx) - xm(nxm-1))*dx*0.5d0

    ! WRITE GRID INFORMATION

    do kc=1,nxm
        udx3m(kc) = dx/dx3c(kc)
        udx3c(kc) = dx/g3rc(kc)
    end do
    udx3c(nx) = dx/g3rc(nx)

    if(ismaster) then
        open(unit=78,file='Results/axicor.out',status='unknown')
        write(78,'(A8,2X,4(A24,2X))') 'kc','xc(kc)','xm(kc)','g3rc(kc)','dx3c(kc)'
        do kc=1,nx
            write(78,'(I8,2X,4(E24.16,2X))') kc,xc(kc),xm(kc),g3rc(kc),dx3c(kc) ! Nomenclature g3rm->dx3c
        end do
        close(78)

        ! QUANTITIES FOR DERIVATIVES

        open(unit=78,file='Results/fact3.out',status='unknown')
        write(78,'(A8,2X,2(A24,2X))') 'kc','udx3m(kc)','udx3c(kc)'
        do kc=1,nxm
            write(78,'(I8,2X,2(E24.16,2X))') kc,udx3m(kc),udx3c(kc)
        end do
            write(78,'(I8,2X,2(E24.16,2X))') nx,udx3m(nxm),udx3c(nx)
        close(78)
    endif
    
    ! VX DIFFERENTIATION (XC VARIABLE)
    
    do kc=2,nxm
        km=kc-1
        kp=kc+1
        ap3ck(kc) = dxq/g3rc(kc)/dx3c(kc)
        am3ck(kc) = dxq/g3rc(kc)/dx3c(km)
        ac3ck(kc) =-(ap3ck(kc)+am3ck(kc))
    enddo
    
    ! Boundary Conditions: Dummy values, not used at all

    am3ck(1)=dxq/g3rc(1)/dx3c(1)        ! 0.0d0
    ap3ck(1)=dxq/g3rc(1)/dx3c(1)        ! 0.0d0
    ac3ck(1)=-(ap3ck(1)+am3ck(1))       ! 1.0d0
    am3ck(nx)=dxq/g3rc(nx)/dx3c(nxm)    ! 0.0d0
    ap3ck(nx)=dxq/g3rc(nx)/dx3c(nxm)    ! 0.0d0
    ac3ck(nx)=-(ap3ck(nx)+am3ck(nx))    ! 1.0d0

    ! VY,VZ DIFFERENTIATION (XM VARIABLE)

    do kc=2,nxm-1
        kp=kc+1
        km=kc-1
        ap3sk(kc) = dxq/g3rm(kc)/dx3m(kc)
        am3sk(kc) = dxq/g3rm(kc)/dx3m(km)
        ac3sk(kc) =-(ap3sk(kc)+am3sk(kc))
    enddo

    ! LOWER WALL BOUNDARY CONDITIONS (INSLWS SETS NO-SLIP vs STRESS-FREE WALL)
    ! kc=1
    ! Note: dx3m(0)=g3rc(1)

    ap3sk(1) = dxq/g3rm(1)/dx3m(1)
    am3sk(1) = dxq/g3rm(1)*2.d0*inslws/g3rc(1)
    ac3sk(1) =-dxq/g3rm(1)*(1.d0/dx3m(1) + 2.d0*inslws/g3rc(1))
    
    ! UPPER WALL BOUNDARY CONDITIONS (INSLWN SETS NO-SLIP vs STRESS-FREE WALL)
    ! kc=nxm
    ! Note: dx3m(nxm)=g3rc(nx)

    ap3sk(nxm) = dxq/g3rm(nxm)*2.d0*inslwn/g3rc(nx)
    am3sk(nxm) = dxq/g3rm(nxm)/dx3m(nxm-1)
    ac3sk(nxm) =-dxq/g3rm(nxm)*(2.d0*inslwn/g3rc(nx) + 1.d0/dx3m(nxm-1))

    ! TEMPERATURE DIFFERENTIATION (XC VARIABLE)
    ! Note: As long Temp defined on XC coefficients

    do kc=2,nxm
        ap3ssk(kc) = ap3ck(kc)
        am3ssk(kc) = am3ck(kc)
        ac3ssk(kc) = ac3ck(kc)
    enddo

    ! Boundary Conditions: Dummy value, not used at all

    am3ssk(1)=am3ck(1)      ! 0.0d0
    ap3ssk(1)=ap3ck(1)      ! 0.0d0
    ac3ssk(1)=ac3ck(1)      ! 1.0d0
    am3ssk(nx)=am3ck(nx)    ! 0.0d0
    ap3ssk(nx)=ap3ck(nx)    ! 0.0d0
    ac3ssk(nx)=ac3ck(nx)    ! 1.0d0

    ! These arrays can be used for first order numerical derviative of
    ! quantities
    ! in wall normal direction with second order accuracy

    ! d     :   differentiation
    ! c     :   coefficient
    ! m/c   :   middle/corner
    ! t/b   :   top/bottom
    ! 1/2/3 :   node number counting from boundary (i.e. 1 is on boundary, 2
    ! is next point from boundary, 3 is second point after boundary)

    ! COEFFICIENTS FOR ONE SIDED 2nd ORDER DIFFERENTIATION AT TOP AND BOTTOM
    ! BOUNDARIES
    
    ! FOR CELL CENTER VARIABLES (i.e. vz and vy)

    dcmb1 = (xc(1)+xc(1)-xm(1)-xm(2))/((xm(1)-xc(1))*(xm(2)-xc(1)))
    dcmb2 = (xc(1)-xm(2))/((xm(1)-xc(1))*(xm(1)-xm(2)))
    dcmb3 = (xm(1)-xc(1))/((xm(2)-xc(1))*(xm(1)-xm(2)))

    dcmt1 = (xc(nx)+xc(nx)-xm(nxm)-xm(nxm-1))/((xm(nxm)-xc(nx))*(xm(nxm-1)-xc(nx)))
    dcmt2 = (xc(nx)-xm(nxm-1))/((xm(nxm)-xc(nx))*(xm(nxm)-xm(nxm-1)))
    dcmt3 = (xm(nxm)-xc(nx))/((xm(nxm-1)-xc(nx))*(xm(nxm)-xm(nxm-1)))

    ! FOR CELL VERTEX OR NODE VARIABLES (i.e. vx and temp)

    dccb1 = (xc(1)+xc(1)-xc(2)-xc(3))/((xc(2)-xc(1))*(xc(3)-xc(1)))
    dccb2 = (xc(1)-xc(3))/((xc(2)-xc(1))*(xc(2)-xc(3)))
    dccb3 = (xc(2)-xc(1))/((xc(3)-xc(1))*(xc(2)-xc(3)))

    dcct1 = (xc(nx)+xc(nx)-xc(nxm)-xc(nxm-1))/((xc(nxm)-xc(nx))*(xc(nxm-1)-xc(nx)))
    dcct2 = (xc(nx)-xc(nxm-1))/((xc(nxm)-xc(nx))*(xc(nxm)-xc(nxm-1)))
    dcct3 = (xc(nxm)-xc(nx))/((xc(nxm-1)-xc(nx))*(xc(nxm)-xc(nxm-1)))

    !! [TESTING]
    if(ismaster) then
        open(unit=78,file='Results/diffc_surf.out',status='unknown')
        write(78,'(6(A24,2X))') 'dcmb1', 'dcmb2', 'dcmb3', 'dcmt1', 'dcmt2', 'dcmt3'
        write(78,'(6(E24.16,2X))') dcmb1, dcmb2, dcmb3, dcmt1, dcmt2, dcmt3
        write(78,'(6(A24,2X))') 'dccb1', 'dccb2', 'dccb3', 'dcct1', 'dcct2', 'dcct3'
        write(78,'(6(E24.16,2X))') dccb1, dccb2, dccb3, dcct1, dcct2, dcct3
        close(78)
    endif
    !! [TESTING] 

    ! d     :   differentiation
    ! c     :   coefficient
    ! p/c/m :   plus/center/minus
    ! dx    :   wall normal direction
    ! m/c   :   middle/corner

    ! COEFFICIENTS FOR ONE 2nd ORDER CENTRAL DIFFERENCE IN WALL NORMAL
    ! DIRECTION (FOR POST PROCESSING)

    do kc=2,nxm

        kp = kc + 1
        km = kc - 1

        dcmdxc(kc) = (xc(kp)-xc(kc))/((xc(km)-xc(kc))*(xc(kp)-xc(km)))
        dccdxc(kc) = ((xc(kc)+xc(kc))-(xc(kp)+xc(km)))/((xc(km)-xc(kc))*(xc(kp)-xc(kc)))
        dcpdxc(kc) = -(xc(km)-xc(kc))/((xc(kp)-xc(kc))*(xc(kp)-xc(km)))

        dcmdxm(kc) = (xm(kp)-xm(kc))/((xm(km)-xm(kc))*(xm(kp)-xm(km)))
        dccdxm(kc) = ((xm(kc)+xm(kc))-(xm(kp)+xm(km)))/((xm(km)-xm(kc))*(xm(kp)-xm(kc)))
        dcpdxm(kc) = -(xm(km)-xm(kc))/((xm(kp)-xm(kc))*(xm(kp)-xm(km)))

    enddo

    ! COEFFICIENTS FOR ONE SIDED 2nd ORDER DIFFERENTIATION AT FIRST AND LAST
    ! GRID POINT

    ! FOR CELL CENTER VARIABLES (i.e. vz and vy)

    dcmdxm(1) = (xm(2)-xm(1))/((xc(1)-xm(1))*(xm(2)-xc(1)))
    dccdxm(1) = ((xm(1)+xm(1))-(xm(2)+xc(1)))/((xc(1)-xm(1))*(xm(2)-xm(1)))
    dcpdxm(1) = -(xc(1)-xm(1))/((xm(2)-xm(1))*(xm(2)-xc(1)))

    dcpdxm(nxm) = -(xm(nxm-1)-xm(nxm))/((xc(nx)-xm(nxm))*(xc(nx)-xm(nxm-1)))
    dccdxm(nxm) = ((xm(nxm)+xm(nxm))-(xc(nx)+xm(nxm-1)))/((xm(nxm-1)-xm(nxm))*(xc(nx)-xm(nxm)))
    dcmdxm(nxm) = (xc(nx)-xm(nxm))/((xm(nxm-1)-xm(nxm))*(xc(nx)-xm(nxm-1)))

    ! FOR CELL VERTEX OR NODE VARIABLES (i.e. vx and temp)

    dcmdxc(1) = dccb1
    dccdxc(1) = dccb2
    dcpdxc(1) = dccb3

    dcpdxc(nx) = dcct1
    dccdxc(nx) = dcct2
    dcmdxc(nx) = dcct3

    !! [TESTING]
    if(ismaster) then
        open(unit=78,file='Results/diffc.out',status='unknown')
        write(78,'(A8,4X,3(A24,2X),2X,3(A24,2X))') 'kc', 'dcpdxc(kc)', 'dccdxc(kc)', 'dcmdxc(kc)', 'dcpdxm(kc)', 'dccdxm(kc)', 'dcmdxm(kc)'
        do kc=1,nx
            write(78,'(I8,4X,3(E24.16,2X),2X,3(E24.16,2X))') kc,dcpdxc(kc),dccdxc(kc),dcmdxc(kc),dcpdxm(kc),dccdxm(kc),dcmdxm(kc)
        end do
        close(78)
    endif
    !! [TESTING]

    return

    end