MODULE GLOBAL
       USE PARAM
       IMPLICIT NONE

       SAVE





! MPI variables
       INTEGER :: myid,ier
       INTEGER :: comm2d
       INTEGER :: n_west, n_east, n_suth, n_nrth
       INTEGER :: npx,npy
       INTEGER :: ndims=2
       INTEGER :: NumberProcessor
       INTEGER, DIMENSION(2) :: dims, coords
       LOGICAL, DIMENSION(2) :: periods
       LOGICAL :: reorder = .true.

       INTEGER :: px,py

! station data
       INTEGER :: NumberStations
       INTEGER,DIMENSION(:),ALLOCATABLE :: ista,jsta,nsta
       REAL(SP):: PLOT_INTV_STATION,PLOT_COUNT_STATION


! define parameters
       CHARACTER(LEN=80) TITLE
       CHARACTER(LEN=80) ReadType
       CHARACTER(LEN=80) DEPTH_TYPE
       CHARACTER(LEN=80) DEPTH_FILE
       CHARACTER(LEN=80) ETA_FILE
       CHARACTER(LEN=80) U_FILE
       CHARACTER(LEN=80) V_FILE
       CHARACTER(LEN=80) DX_FILE
       CHARACTER(LEN=80) DY_FILE

       CHARACTER(LEN=80) X_FILE
       CHARACTER(LEN=80) Y_FILE

       CHARACTER(LEN=80) LATITUDE_FILE
       CHARACTER(LEN=80) OBSTACLE_FILE
       CHARACTER(LEN=80) RESULT_FOLDER
       CHARACTER(LEN=80) STATIONS_FILE
       CHARACTER(LEN=80) WaveMaker
       CHARACTER(LEN=80) Time_Scheme
       CHARACTER(LEN=80) CONSTR
       CHARACTER(LEN=80) HIGH_ORDER

       REAL(SP),PARAMETER,DIMENSION(3)::alpha=(/0.0_SP,3.0_SP/4.0_SP,1.0_SP/3.0_SP/)
       REAL(SP),PARAMETER,DIMENSION(3)::beta=(/1.0_SP,1.0_SP/4.0_SP,2.0_SP/3.0_SP/)

! some global variables
! Mloc1=Mloc+1, Nloc1=Nloc+1
       INTEGER :: Mglob,Nglob,Mloc,Nloc,Mloc1,Nloc1
       INTEGER, PARAMETER :: Nghost = 3
       INTEGER :: Ibeg,Iend,Jbeg,Jend,Iend1,Jend1

       LOGICAL :: CORI_CONSTANT

       REAL(SP):: DX,DY,LATITUDE
       REAL(SP)::  DT,TIME,TOTAL_TIME,PLOT_INTV,PLOT_COUNT,&
                  SCREEN_INTV,SCREEN_COUNT
       REAL(SP) :: HOTSTART_INTV,HOTSTART_COUNT
       INTEGER :: icount=0  ! for output file number
       INTEGER :: icount_riverdata = 1
       INTEGER :: icount_hotstart=0
       INTEGER :: FileNumber_HOTSTART
       
       REAL(SP) :: MinDepth=0.001
       REAL(SP) :: MinDepthFrc=0.5
       REAL(SP) :: CFL=0.15_SP
       REAL(SP) :: FroudeCap=10.0_SP
       REAL(SP) :: SWE_ETA_DEP=0.7_SP
       REAL(SP) :: Cd=0.005
       REAL(SP) :: Manning=0.0_SP
       REAL(SP),DIMENSION(:,:),ALLOCATABLE :: nu_total 
                         ! eddy viscosity (nu_0+nu_b+nu_smg)
       REAL(SP) :: nu_0=0.0_SP
       REAL(SP) :: nu_bkgd=0.0_SP
! curvilinear

       REAL(SP),DIMENSION(:,:),ALLOCATABLE :: &
       X,Y,Jaco,JL11x,JL12x,JL21x,JL22x, JS1,JS2,JacoX,JacoY,&
       JL11y,JL12y,JL21y,JL22y,L11,L12,L21,L22, &
       Uc,Vc,UnL,UnR,VnL,VnR,UbarXL,UbarXR,UbarYL,UbarYR,&
       VbarXL,VbarXR,VbarYL,VbarYR,DelxUbar,DelxVbar,&
       DelyUbar,DelyVbar

! stationary mode for shorecirc




! some local variables
       REAL(SP),DIMENSION(:,:),ALLOCATABLE :: &
          Coriolis,&
          DelxU,DelxHU,DelxV,DelxEtar,&
          DelxHV, DelyHU, &
          DelyU,DelyHV,DelyV,DelyEtar,&
          UxL,UxR,VxL,VxR,&
          HUxL,HUxR,HUyL,HUyR,HxL,HxR, &
          EtaRxL,EtaRxR,&
          UyL,UyR,VyL,VyR,&
          HVxL,HVxR,HVyL,HVyR,HyL,HyR, &
          EtaRyL,EtaRyR, &
          PL,PR,QL,QR, &
          FxL,FxR,FyL,FyR, &
          GxL,GxR,GyL,GyR, &
          SxL,SxR,SyL,SyR, &
! original variables
          Fx,Fy,U,V,HU,HV,&
          Gx,Gy,P,Q,SourceX,SourceY,Int2Flo, &
          tmp4preview,HeightMax
! quasi-3d stuff
       REAL(SP),DIMENSION(:,:),ALLOCATABLE :: dd11,ee11,ff11,dd12,ee12,ff12
       REAL(SP),DIMENSION(:,:),ALLOCATABLE :: Ubott,Vbott,Usurf,Vsurf

! friction manning
       REAL(SP),DIMENSION(:,:),ALLOCATABLE :: CD_2D
       LOGICAL :: FRC_MANNING_DATA = .FALSE.
       CHARACTER(LEN=80) FRC_FILE       




! open boundary variables
       LOGICAL :: ETA_CLAMPED = .FALSE.,FLUX_TIDE=.FALSE.
       REAL(SP),DIMENSION(:),ALLOCATABLE :: ETA_bd
       INTEGER,DIMENSION(:),ALLOCATABLE :: I_eta_bd,J_eta_bd
       CHARACTER(LEN=80) :: TIDE_FILE, FLUX_FILE,FLUX_TIDE_FILE
       INTEGER :: NumConstituent,NumEtaPoint
       REAL(SP),DIMENSION(:,:),ALLOCATABLE :: TidePeriod,TideAmp,TidePha
       REAL(SP),DIMENSION(:,:),ALLOCATABLE :: FluxU_tide,FluxPeriod,FluxPha
       REAL(SP),DIMENSION(:),ALLOCATABLE :: FluxAngle_Tide
       REAL(SP),DIMENSION(:),ALLOCATABLE :: TideFac,Tideu0
       REAL(SP) :: TideStartDate
       CHARACTER(LEN=1),DIMENSION(:),ALLOCATABLE  :: FLUX_ORIENTATION
  ! flux file
       LOGICAL :: FLUX_CLAMPED = .FALSE.
       INTEGER :: NumFluxPoint,NumTimeData  
       REAL(SP),DIMENSION(:,:),ALLOCATABLE :: FluxU_bd,FluxAngle
       CHARACTER(LEN=1),DIMENSION(:),ALLOCATABLE  :: RIVER_ORIENTATION
       INTEGER,DIMENSION(:),ALLOCATABLE :: I_flux_bd,J_flux_bd 
       REAL(SP),DIMENSION(:),ALLOCATABLE :: TimeFlux
   ! wind file
       LOGICAL :: WindForce = .FALSE.
       CHARACTER(LEN=80) WIND_FILE
       INTEGER :: NumTimeWindData
       REAL(SP), DIMENSION(:,:),ALLOCATABLE :: WindU2D,WindV2D
       REAL(SP),DIMENSION(:),ALLOCATABLE :: TimeWind
       REAL(SP),DIMENSION(:),ALLOCATABLE :: WU,WV
       REAL(SP) :: Cdw   
       INTEGER :: icount_winddata = 1       
! wetting and drying
        INTEGER,DIMENSION(:,:),ALLOCATABLE :: MASK,MASK_STRUC,MASK9
! wave maker
        REAL(SP)::AMP_SOLI,DEP_SOLI,LAG_SOLI, CPH_SOLI,XWAVEMAKER, &
                  Xc,Yc, WID

        REAL(SP)::x1_Nwave = 5.0, &
                  x2_Nwave = 5.0, &
                  a0_Nwave = 1.0, &
                  gamma_Nwave = -3.0, &
                  dep_Nwave = 1.0
! sponge
        REAL(SP),DIMENSION(:,:),ALLOCATABLE :: SPONGE
        REAL(SP)::Sponge_west_width,Sponge_east_width, &
                  Sponge_south_width,Sponge_north_width, &
                  R_sponge,A_sponge


        LOGICAL :: INI_UVZ=.FALSE.

! smagorinsky and mixing

      REAL(SP),DIMENSION(:,:),ALLOCATABLE :: Umean,Vmean,&
                  UmeanL,VmeanL,UsumL,VsumL, &
                  UmeanS,VmeanS,UsumS,VsumS, &
                  UmeanB,VmeanB,UsumB,VsumB, &
                  ETAmean,Usum,Vsum,ETAsum, nu_smg, &
                  Tau_bx,Tau_by,Tau_sx,Tau_sy, &
                  D_11,D_22,D_12
      REAL(SP) :: Disp     
      LOGICAL :: DISP_EMPIRICAL = .TRUE.
            
      REAL(SP)::T_INTV_mean = 20.0,T_sum=0.0,C_smg=0.25

! friction related
      REAL(SP) :: beta1,beta2


! depth H=Eta+Depth, 
       REAL(SP),DIMENSION(:,:),ALLOCATABLE :: Depth,H,&
             Depthx,Depthy
       REAL(SP)::Depth_Flat, SLP,Xslp

! updating variables
       REAL(SP),DIMENSION(:,:),ALLOCATABLE::Ubar0,Vbar0,Eta0,&
                               Ubar,Vbar,Eta
! output logical parameters
       LOGICAL :: OUT_U=.FALSE., OUT_V=.FALSE.,OUT_ETA=.FALSE., &
                  OUT_MASK=.FALSE., &
                  OUT_SourceX=.FALSE., OUT_SourceY=.FALSE., &
                  OUT_MASK9=.FALSE., OUT_DEPTH=.FALSE., &
                  OUT_TMP=.FALSE., OUT_AGE=.FALSE., &
                  OUT_HS=.FALSE.,OUT_PER=.FALSE., &
                  OUT_WFC=.FALSE.,OUT_WBV=.FALSE.,&
                  OUT_WDIR=.FALSE.,OUT_3D=.FALSE., &
                  OUT_Qw=.FALSE.,OUT_Wdis=.FALSE.

! periodic boundary conditions
       LOGICAL :: PERIODIC_X=.FALSE.,PERIODIC_Y=.FALSE.
       INTEGER :: Num_Transit = 0
! dispersion control
       LOGICAL :: OBSTACLE=.FALSE., HOT_START=.FALSE.
! sponge
       LOGICAL :: SPONGE_ON=.FALSE.
! breaking
!       LOGICAL :: SHOW_BREAKING=.TRUE.




END MODULE GLOBAL

