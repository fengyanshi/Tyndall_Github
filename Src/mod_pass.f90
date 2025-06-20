MODULE PASS
       USE PARAM
       IMPLICIT NONE

       INTEGER :: IX,IY,INDX
! pass variables
       REAL(SP),DIMENSION(:,:),ALLOCATABLE :: &
         WaveHeightSW,WaveHeightSC,WaveHeightGL,&
         PeakPeriodSW,PeakPeriodSC,PeakPeriodGL, &
         AvePeriodSW,AvePeriodSC,AvePeriodGL, &
         WaveAngleSW, WaveAngleSC, WaveAngleGL, &
         PeakAngleSW, PeakAngleSC, PeakAngleGL, &
         WaveUbottSW, WaveUbottSC, WaveUbottGL, &
         WaveDissSW, WaveDissSC, WaveDissGL, &
         WaveDissFrcSW,WaveDissFrcSC,WaveDissFrcGL,&
         WaveDissWcapSW,WaveDissWcapSC,WaveDissWcapGL,&
         WaveFluxXSW,WaveFluxXSC,WaveFluxXGL, &
         WaveFluxYSW,WaveFluxYSC,WaveFluxYGL, &
         WaveFxSW,WaveFxSC,WaveFxGL, &
         WaveFySW,WaveFySC,WaveFyGL, &
         WaveBrFraSW, WaveBrFraGL,WaveBrFraSC,&
         Uswan,Vswan, &
         EtaSwan

! other variable needed
       REAL(SP),DIMENSION(:,:),ALLOCATABLE :: &
         GRID_ANGLE_SW2SC,MASK_WC_INTERACT

       REAL(SP) :: WC_LAG

! number of grid points for wave-currnt interaction region
       INTEGER :: WC_BOUND_WEST,WC_BOUND_EAST, &
                  WC_BOUND_SOUTH,WC_BOUND_NORTH

!       INTEGER,DIMENSION(:),ALLOCATABLE :: VOQR
!       REAL(SP),DIMENSION(:),ALLOCATABLE :: VOQ
         

      INTEGER,DIMENSION(:),ALLOCATABLE :: MXFs,MYFs,& 
             IXTRIM_Ls,IXTRIM_Rs,IYTRIM_Ls,IYTRIM_Rs     
      INTEGER :: IXTRIM_L,IXTRIM_R,IYTRIM_L,IYTRIM_R  
      LOGICAL :: FIRST_CALL_SW2SC = .TRUE.
      LOGICAL :: FIRST_CALL_SC2SW = .TRUE.


      LOGICAL :: SWAN_RUN = .TRUE.
      LOGICAL :: SHORECIRC_RUN = .TRUE.
      LOGICAL :: FIRST_CALL_PASS_SW = .TRUE.
      LOGICAL :: FIRST_CALL_PASS_SC = .TRUE.
     REAL(SP) :: DT_SC

END MODULE PASS


