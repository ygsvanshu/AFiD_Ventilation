
! Declaration of global variables
!***********************************************************  

module param

    implicit none

    !==========================================================			
    !       read from input file bou.in
    !==========================================================

    logical         :: readflow=.false.
    logical         :: statread=.false.

    integer         :: nx, ny, nz
    real            :: alx3,ylen,zlen
    integer         :: istr3
    real            :: str3

    integer         :: ntst
    real            :: walltimemax,tmax
    
    integer         :: nsst
    logical         :: variabletstep=.true.
    real            :: dt,dtmin,dtmax,limitCFL,resid,limitVel,eps
    
    real            :: ray,pra,lambda_co2,lambda_h2o
    real            :: iheight,ilen,oheight,olen
    real            :: ivel,tvel
    real            :: ocou,ovsc,odst
    logical         :: person_on=.false.
    logical         :: breath_on=.false.
    character*50    :: person_objfile
    real            :: personx,persony,personz,sclf
    real            :: breathx,breathy,breathz
    real            :: kernel_width_space,kernel_width_time
    integer         :: inslws,inslwn

    logical         :: statcalc=.false.
    real            :: tsta
    integer         :: nout
    real            :: stats2Dx,stats2Dy,stats2Dz

    logical         :: savesnap=.false.
    real            :: startsnap,tsnap

    logical         :: savemovie=.false.
    real            :: tframe
    real            :: movie2Dx,movie2Dy,movie2Dz
    
    real            :: tcontinua

    !=================================================
    !       end of input file
    !=================================================

    real :: time

    !******* Grid parameters**************************
    real :: dx,dy,dz,dxq,dyq,dzq        
    real, allocatable, dimension(:) :: xc,xm
    real, allocatable, dimension(:) :: yc,ym
    real, allocatable, dimension(:) :: zc,zm
    real, allocatable, dimension(:) :: g3rc,g3rm
    real, allocatable, dimension(:) :: dx3c,dx3m ! Additional Grid parameters
    real, allocatable, dimension(:) :: udx3c,udx3m
    integer, allocatable, dimension(:) :: kmc,kpc,kmv,kpv
    real, allocatable, dimension(:) :: ap3ck,ac3ck,am3ck
    real, allocatable, dimension(:) :: ap3sk,ac3sk,am3sk
    real, allocatable, dimension(:) :: ap3ssk,ac3ssk,am3ssk   

    !******* Variables for FFTW and Poisson solver****************
    real, allocatable, dimension(:) :: ak2,ap
    real, allocatable, dimension(:) :: ak1,ao
    real, allocatable, dimension(:) :: amphk,acphk,apphk
    
    !******* Other variables ***********************************
    integer                 :: nxm, nym, nzm
    real                    :: ren, pec
    real                    :: pi
    real                    :: al,ga,ro
    real                    :: beta
    real                    :: qqmax,qqtot
    real                    :: re
    real                    :: tempmax,tempmin,tempm
    real                    :: co2max,co2min,co2m
    real                    :: h2omax,h2omin,h2om
    integer                 :: ntime
    integer, parameter      :: ndv=3
    real, dimension(1:ndv)  :: vmax
    real, dimension(1:3)    :: gam,rom,alm
    logical                 :: ismaster=.false.
    integer                 :: lvlhalo=1
    integer                 :: tsteps
    integer                 :: comm_xcut,comm_ycut,comm_zcut
        
end module param
      
!************* End of param module******************************
!===============================================================
!******* 2D arrays, dynamically allocated by each process*******

module local_arrays
    use param
    implicit none
    real,allocatable,dimension(:,:,:)   :: vx,vy,vz
    real,allocatable,dimension(:,:,:)   :: pr,temp,co2,h2o,rhs
    real,allocatable,dimension(:,:,:)   :: rux,ruy,ruz,rutemp,ruco2,ruh2o
    real,allocatable,dimension(:,:,:)   :: dph,qcap,dq,hro,dphhalo,qco2,qh2o
    real,allocatable,dimension(:,:,:)   :: ibm_body
end module local_arrays

module stat_arrays

    implicit none

    integer                             :: nstatsamples
    real                                :: tstat
    real                                :: tinterval

    integer                             :: stat_xi,stat_xj,stat_xk
    integer                             :: stat_yi,stat_yj,stat_yk
    integer                             :: stat_zi,stat_zj,stat_zk

    real,allocatable,dimension(:,:)     :: vx_m1_xcut,vy_m1_xcut,vz_m1_xcut,temp_m1_xcut,co2_m1_xcut,h2o_m1_xcut
    real,allocatable,dimension(:,:)     :: vx_m2_xcut,vy_m2_xcut,vz_m2_xcut,temp_m2_xcut,co2_m2_xcut,h2o_m2_xcut
    real,allocatable,dimension(:,:)     :: vx_m3_xcut,vy_m3_xcut,vz_m3_xcut,temp_m3_xcut,co2_m3_xcut,h2o_m3_xcut
    real,allocatable,dimension(:,:)     :: vx_m4_xcut,vy_m4_xcut,vz_m4_xcut,temp_m4_xcut,co2_m4_xcut,h2o_m4_xcut

    real,allocatable,dimension(:,:)     :: vx_m1_ycut,vy_m1_ycut,vz_m1_ycut,temp_m1_ycut,co2_m1_ycut,h2o_m1_ycut
    real,allocatable,dimension(:,:)     :: vx_m2_ycut,vy_m2_ycut,vz_m2_ycut,temp_m2_ycut,co2_m2_ycut,h2o_m2_ycut
    real,allocatable,dimension(:,:)     :: vx_m3_ycut,vy_m3_ycut,vz_m3_ycut,temp_m3_ycut,co2_m3_ycut,h2o_m3_ycut
    real,allocatable,dimension(:,:)     :: vx_m4_ycut,vy_m4_ycut,vz_m4_ycut,temp_m4_ycut,co2_m4_ycut,h2o_m4_ycut

    real,allocatable,dimension(:,:)     :: vx_m1_zcut,vy_m1_zcut,vz_m1_zcut,temp_m1_zcut,co2_m1_zcut,h2o_m1_zcut
    real,allocatable,dimension(:,:)     :: vx_m2_zcut,vy_m2_zcut,vz_m2_zcut,temp_m2_zcut,co2_m2_zcut,h2o_m2_zcut
    real,allocatable,dimension(:,:)     :: vx_m3_zcut,vy_m3_zcut,vz_m3_zcut,temp_m3_zcut,co2_m3_zcut,h2o_m3_zcut
    real,allocatable,dimension(:,:)     :: vx_m4_zcut,vy_m4_zcut,vz_m4_zcut,temp_m4_zcut,co2_m4_zcut,h2o_m4_zcut

end module stat_arrays

module movie_indices
    implicit none
    integer                             :: mov_xi,mov_xj,mov_xk
    integer                             :: mov_yi,mov_yj,mov_yk
    integer                             :: mov_zi,mov_zj,mov_zk
    integer                             :: mov_ci,mov_cj,mov_ck
end module movie_indices

module vent_arrays
    implicit none
    real                                :: isvel,tsvel
    integer                             :: ixcst,ixcen
    integer                             :: ixmst,ixmen
    integer                             :: ixfst,ixfen
    integer                             :: oxcst,oxcen
    integer                             :: oxmst,oxmen
    integer                             :: oxfst,oxfen
    real                                :: iflux,oflux
    real                                :: iarea,oarea
    real                                :: varea,harea
    real,allocatable,dimension(:)       :: icell,ocell
    real,allocatable,dimension(:)       :: igrid,ogrid
    real,allocatable,dimension(:,:)     :: outvx,outvy,outvz
    real,allocatable,dimension(:,:)     :: outtemp,outco2,outh2o
    real,allocatable,dimension(:,:)     :: outvscx,outvscy,outvscz
end module vent_arrays

!=====================================================    

module stat3_param
    implicit none
    integer :: kslab(1:9)
    real    :: xslab(1:9)
end module stat3_param

!=====================================================    

module mpih
    implicit none
    include 'mpif.h'
    integer :: ierr
    integer, parameter :: master=0
    integer :: MDP = MPI_DOUBLE_PRECISION
end module mpih

!====================================================
module fftw_params
    ! use param, only: m2m,m2mh,m1m
    use iso_c_binding

    type, bind(C) :: fftw_iodim
        integer(C_INT) n, is, os
    end type fftw_iodim

    integer, parameter :: C_FFTW_R2R_KIND = C_INT32_T

    interface

        type(C_PTR) function fftw_plan_guru_dft(rank,dims,howmany_rank,howmany_dims,in,out,sign,flags) bind(C, name='fftw_plan_guru_dft')
            import
            integer(C_INT), value :: rank
            type(fftw_iodim), dimension(*), intent(in) :: dims
            integer(C_INT), value :: howmany_rank
            type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
            complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
            complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
            integer(C_INT), value :: sign
            integer(C_INT), value :: flags
        end function fftw_plan_guru_dft

        type(C_PTR) function fftw_plan_guru_dft_r2c(rank,dims,howmany_rank,howmany_dims,in,out,flags) bind(C, name='fftw_plan_guru_dft_r2c')
            import
            integer(C_INT), value :: rank
            type(fftw_iodim), dimension(*), intent(in) :: dims
            integer(C_INT), value :: howmany_rank
            type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
            real(C_DOUBLE), dimension(*), intent(out) :: in
            complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
            integer(C_INT), value :: flags
        end function fftw_plan_guru_dft_r2c
        
        type(C_PTR) function fftw_plan_guru_dft_c2r(rank,dims,howmany_rank,howmany_dims,in,out,flags) bind(C, name='fftw_plan_guru_dft_c2r')
            import
            integer(C_INT), value :: rank
            type(fftw_iodim), dimension(*), intent(in) :: dims
            integer(C_INT), value :: howmany_rank
            type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
            complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
            real(C_DOUBLE), dimension(*), intent(out) :: out
            integer(C_INT), value :: flags
        end function fftw_plan_guru_dft_c2r

        type(C_PTR) function fftw_plan_guru_r2r(rank,dims,howmany_rank,howmany_dims,in,out,kind,flags) bind(C, name='fftw_plan_guru_r2r')
            import
            integer(C_INT), value :: rank
            type(fftw_iodim), dimension(*), intent(in) :: dims
            integer(C_INT), value :: howmany_rank
            type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
            real(C_DOUBLE), dimension(*), intent(out) :: in
            real(C_DOUBLE), dimension(*), intent(out) :: out
            integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kind
            integer(C_INT), value :: flags
        end function fftw_plan_guru_r2r

    end interface

    integer FFTW_PATIENT, FFTW_FORWARD, FFTW_BACKWARD,FFTW_ESTIMATE
    integer FFTW_REDFT01, FFTW_REDFT10
    parameter (FFTW_PATIENT=32)
    parameter (FFTW_MEASURE=0)
    parameter (FFTW_ESTIMATE=64)   
    parameter (FFTW_FORWARD=-1)   
    parameter (FFTW_BACKWARD=1)
    parameter (FFTW_REDFT01=4)
    parameter (FFTW_REDFT10=5)

    type(C_PTR) :: fwd_guruplan_y,bwd_guruplan_y 
    type(C_PTR) :: fwd_guruplan_z,bwd_guruplan_z
    logical :: planned=.false.

    real,allocatable,dimension(:,:,:) :: ry1,rz1
    complex,allocatable,dimension(:,:,:) :: cy1,cz1,dphc

    complex,allocatable,dimension(:,:,:) :: fouvar1,fouvar2

    real,allocatable,dimension(:,:,:) :: ry2,rz2
    real,allocatable,dimension(:,:,:) :: dphr 

end module fftw_params