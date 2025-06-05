SUBROUTINE GET_COUPLING_NEEDS (COMPDA_C,KGRPNT_C,XCGRID_C,YCGRID_C,  &
                 MXK,MYK,VOQR_C,VOQ_C)
      USE OUTP_DATA                                                       
      USE SWAN_COMMON
      IMPLICIT NONE
      REAL COMPDA_C(MCGRD,MCMVAR),XCGRID_C(MXC,MYC),YCGRID_C(MXC,MYC)
      integer KGRPNT_C(MXC,MYC),MXK,MYK
      integer VOQR_C(*)
      real VOQ_C(MXK*MYK,*)

       write(*,*) MXC,MYC,MXK,MYK

END SUBROUTINE GET_COUPLING_NEEDS

SUBROUTINE SHORECIRC2SWAN
!      USE SWAN_COMMON, ONLY : IX,IY,INODE,MXC,MYC,FIRST_CALL_PASS_SC
      USE SWAN_COMMON
      USE PARAM
      USE GLOBAL, ONLY : Mloc,Nloc,U,V,ETA,Mglob,Nglob,TIME,MASK
      IMPLICIT NONE
      REAL(SP),DIMENSION(Mglob,Nglob)::VaGlob
      REAL(SP),DIMENSION(Mloc,Nloc)::ETA_MASKED
      REAL(SP),DIMENSION(MXC,MYC)::EtaSW,USW,VSW

        IF(FIRST_CALL_PASS_SC)THEN
          FIRST_CALL_PASS_SC = .FALSE.
          IF(.NOT.ALLOCATED(GRID_ANGLE_SW2SC)) &
             ALLOCATE(GRID_ANGLE_SW2SC(MXC,MYC))
          GRID_ANGLE_SW2SC = ZERO
          DO IY=1,MYC
          DO IX=1,MXC-1
           GRID_ANGLE_SW2SC(IX,IY)= &
           ATAN2(YCGRID(IX+1,IY)-YCGRID(IX,IY),XCGRID(IX+1,IY)-XCGRID(IX,IY))
          ENDDO
          ENDDO
          DO IY=1,MYC
           GRID_ANGLE_SW2SC(MXC,IY)=GRID_ANGLE_SW2SC(MXC-1,IY)
          ENDDO

          VaGlob = ZERO
          EtaSW = ZERO
          USW = ZERO
          VSW = ZERO
        ENDIF

       IF(TIME.GE.WC_LAG)THEN 

! there's a bug here, dry eta can pass to swan, fixed 01/22/2012
! it's not a good idea not to pass eta outside mask_wc_interact 04/12/2013
! so I remove it
        ETA_MASKED=ETA*MASK
        CALL GATHER_SC2GLOBAL(ETA_MASKED,VaGlob)

        VaGlob=VaGlob    ! I removed *MASK_WC_INTERACT 

        CALL DISTRIBUTE2SWAN(VaGlob,EtaSW)

        CALL GATHER_SC2GLOBAL(U,VaGlob)


        VaGlob=VaGlob*MASK_WC_INTERACT 

        CALL DISTRIBUTE2SWAN(VaGlob,USW)

        CALL GATHER_SC2GLOBAL(V,VaGlob)

        VaGlob=VaGlob*MASK_WC_INTERACT 

        CALL DISTRIBUTE2SWAN(VaGlob,VSW)

! end parallel

! assign to previous time level
        DO INDX = 2, MCGRD
          COMPDA(INDX,JWLV1)=COMPDA(INDX,JWLV2)
          COMPDA(INDX,JVX1)=COMPDA(INDX,JVX2)
          COMPDA(INDX,JVY1)=COMPDA(INDX,JVY2)
        ENDDO
         DO IX=1,MXC
         DO IY=1,MYC
          INDX=KGRPNT(IX,IY)
          IF(INDX.GT.1)THEN
           COMPDA(INDX,JWLV2)=EtaSW(IX,IY)
           COMPDA(INDX,JVX2)=USW(IX,IY)*COS(GRID_ANGLE_SW2SC(IX,IY)) &
                            +VSW(IX,IY)*SIN(GRID_ANGLE_SW2SC(IX,IY))
           COMPDA(INDX,JVY2)=VSW(IX,IY)*COS(GRID_ANGLE_SW2SC(IX,IY)) &
                            -USW(IX,IY)*SIN(GRID_ANGLE_SW2SC(IX,IY))
          ENDIF
         ENDDO
         ENDDO

       ENDIF ! WC_LAG

!      IF(INODE.EQ.1)THEN
!        OPEN(2,FILE='tmp1.txt')
!         DO IY=1,MYC
!          WRITE(2,1221)(EtaSW(IX,IY),IX=1,MXC)
!         ENDDO
!        CLOSE(2)
!      ENDIF
1221     format(500f12.6)

END SUBROUTINE SHORECIRC2SWAN


SUBROUTINE GATHER_SC2GLOBAL(PHI,PHIGLOB)
      USE PARAM
      USE GLOBAL,ONLY : ier,NumberProcessor,MGlob,NGlob,Nghost,Mloc,Nloc, &
                       npx,npy,myid,px,py,Ibeg,Iend,Jbeg,Jend
      USE PASS
     IMPLICIT NONE
     INTEGER :: l
     INTEGER,DIMENSION(NumberProcessor) :: npxs,npys
     REAL(SP),DIMENSION(NumberProcessor) :: xx
     REAL(SP),DIMENSION(MGlob,NGlob),INTENT(OUT) :: PHIGLOB
     REAL(SP),DIMENSION(MGlob+2*Nghost,NGlob+2*Nghost) :: PHIGLOB_GHOST
     REAL(SP),DIMENSION(Mloc,Nloc),INTENT(IN) :: PHI

     call MPI_Gather(npx,1,MPI_INTEGER,npxs,1,MPI_INTEGER,&
          0,MPI_COMM_WORLD,ier)
     call MPI_Gather(npy,1,MPI_INTEGER,npys,1,MPI_INTEGER,&
          0,MPI_COMM_WORLD,ier)

     do i=1,Mloc
     do j=1,Nloc
        call MPI_Gather(PHI(i,j),1,MPI_SP,&
             xx,1,MPI_SP,0,MPI_COMM_WORLD,ier)

        if (j.eq.1) call MPI_Barrier(MPI_COMM_WORLD,ier)

        if (myid.eq.0) then
           do l=1,px*py
              PHIGLOB_GHOST(i+npxs(l)*(Iend-Ibeg+1),&
                   j+npys(l)*(Jend-Jbeg+1)) = xx(l)
           enddo
        endif
     enddo
     enddo

     DO J=1,NGlob
     DO I=1,MGlob
       PHIGLOB(I,J)=PHIGLOB_GHOST(I+Nghost,J+Nghost)
     ENDDO
     ENDDO

END SUBROUTINE GATHER_SC2GLOBAL


SUBROUTINE SWAN2SHORECIRC
      USE SWAN_COMMON
      USE OUTP_DATA
      USE GLOBAL, ONLY : Mloc,Nloc,TIME,H,MinDepthFrc,tmp4preview
      IMPLICIT NONE

        IF(FIRST_CALL_PASS_SW)THEN
         FIRST_CALL_PASS_SW = .FALSE.
! SHORECIRC INITIALIZATION

        IF(.NOT.ALLOCATED(WaveHeightSW)) ALLOCATE(WaveHeightSW(MXC,MYC))
        IF(.NOT.ALLOCATED(WaveHeightSC)) ALLOCATE(WaveHeightSC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(PeakPeriodSW)) ALLOCATE(PeakPeriodSW(MXC,MYC))
        IF(.NOT.ALLOCATED(PeakPeriodSC)) ALLOCATE(PeakPeriodSC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(AvePeriodSW)) ALLOCATE(AvePeriodSW(MXC,MYC))
        IF(.NOT.ALLOCATED(AvePeriodSC)) ALLOCATE(AvePeriodSC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(WaveAngleSW)) ALLOCATE(WaveAngleSW(MXC,MYC))
        IF(.NOT.ALLOCATED(WaveAngleSC)) ALLOCATE(WaveAngleSC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(PeakAngleSW)) ALLOCATE(PeakAngleSW(MXC,MYC))
        IF(.NOT.ALLOCATED(PeakAngleSC)) ALLOCATE(PeakAngleSC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(WaveUbottSW)) ALLOCATE(WaveUbottSW(MXC,MYC))
        IF(.NOT.ALLOCATED(WaveUbottSC)) ALLOCATE(WaveUbottSC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(WaveFxSW)) ALLOCATE(WaveFxSW(MXC,MYC))
        IF(.NOT.ALLOCATED(WaveFxSC)) ALLOCATE(WaveFxSC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(WaveFySW)) ALLOCATE(WaveFySW(MXC,MYC))
        IF(.NOT.ALLOCATED(WaveFySC)) ALLOCATE(WaveFySC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(WaveFluxXSW)) ALLOCATE(WaveFluxXSW(MXC,MYC))
        IF(.NOT.ALLOCATED(WaveFluxXSC)) ALLOCATE(WaveFluxXSC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(WaveFluxYSW)) ALLOCATE(WaveFluxYSW(MXC,MYC))
        IF(.NOT.ALLOCATED(WaveFluxYSC)) ALLOCATE(WaveFluxYSC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(WaveDissSW)) ALLOCATE(WaveDissSW(MXC,MYC))
        IF(.NOT.ALLOCATED(WaveDissSC)) ALLOCATE(WaveDissSC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(WaveBrFraSW)) ALLOCATE(WaveBrFraSW(MXC,MYC))
        IF(.NOT.ALLOCATED(WaveBrFraSC)) ALLOCATE(WaveBrFraSC(Mloc,Nloc))

        WaveHeightSW = ZERO
        WaveHeightSC = ZERO
        PeakPeriodSW = ZERO
        PeakPeriodSC = ZERO
        WaveAngleSW  = ZERO
        WaveAngleSC  = ZERO
        PeakAngleSW  = ZERO
        PeakAngleSC  = ZERO
        WaveUbottSW  = ZERO
        WaveUbottSC  = ZERO
        WaveFxSW = ZERO
        WaveFxSC = ZERO
        WaveFySW = ZERO
        WaveFySC = ZERO
        WaveFluxXSW = ZERO
        WaveFluxYSW = ZERO
        WaveFluxXSC = ZERO
        WaveFluxYSC = ZERO
        WaveDissSW = ZERO
        WaveDissSC = ZERO
        WaveBrFraSW = ZERO
        WaveBrFraSC = ZERO


        IF(.NOT.ALLOCATED(WaveHeightGL))ALLOCATE(WaveHeightGL(MXCGL,MYCGL))
        IF(.NOT.ALLOCATED(PeakPeriodGL))ALLOCATE(PeakPeriodGL(MXCGL,MYCGL))
        IF(.NOT.ALLOCATED(AvePeriodGL))ALLOCATE(AvePeriodGL(MXCGL,MYCGL))
        IF(.NOT.ALLOCATED(WaveAngleGL))ALLOCATE(WaveAngleGL(MXCGL,MYCGL))
        IF(.NOT.ALLOCATED(PeakAngleGL))ALLOCATE(PeakAngleGL(MXCGL,MYCGL))
        IF(.NOT.ALLOCATED(WaveUbottGL))ALLOCATE(WaveUbottGL(MXCGL,MYCGL))
        IF(.NOT.ALLOCATED(WaveFxGL))ALLOCATE(WaveFxGL(MXCGL,MYCGL))
        IF(.NOT.ALLOCATED(WaveFyGL))ALLOCATE(WaveFyGL(MXCGL,MYCGL))
        IF(.NOT.ALLOCATED(WaveFluxXGL))ALLOCATE(WaveFluxXGL(MXCGL,MYCGL))
        IF(.NOT.ALLOCATED(WaveFluxYGL))ALLOCATE(WaveFluxYGL(MXCGL,MYCGL))
        IF(.NOT.ALLOCATED(WaveDissGL))ALLOCATE(WaveDissGL(MXCGL,MYCGL))
        IF(.NOT.ALLOCATED(WaveBrFraGL))ALLOCATE(WaveBrFraGL(MXCGL,MYCGL))

        WaveHeightGL = ZERO
        PeakPeriodGL = ZERO
        WaveAngleGL  = ZERO
        WaveUbottGL  = ZERO
        WaveFxGL  = ZERO
        WaveFyGL  = ZERO
        WaveFluxXGL = ZERO
        WaveFluxYGL = ZERO
        WaveDissGL = ZERO
        WaveBrFraGL = ZERO



        ENDIF ! end first_call_pass

! wave height and bottom orbital velocity, dissipation
! note: for dissipation in SWAN, 
! JDISS - total dissipation
! JDSXB - bottom friction dissipation
! JDSXS - surfzone breaking
! JDSXW - whitecapping dissipation
! JQB   - wave breaking fractoin
         DO IX=1,MXC
         DO IY=1,MYC
          INDX=KGRPNT(IX,IY)
          WaveHeightSW(IX,IY)=COMPDA(INDX,JHS)
          WaveUbottSW(IX,IY) =COMPDA(INDX,JUBOT)
          WaveDissSW(IX,IY) = COMPDA(INDX,JDSXS)
          WaveBrFraSW(IX,IY) = COMPDA(INDX,JQB)
         ENDDO
         ENDDO



! wave height
        CALL GATHER_SWAN2GLOBAL(WaveHeightSW,WaveHeightGL)



        CALL DISTRIBUTE2SHORECIRC(WaveHeightGL,WaveHeightSC)

! wave angle
        CALL GATHER_SWAN2GLOBAL(WaveAngleSW,WaveAngleGL)

        CALL DISTRIBUTE2SHORECIRC(WaveAngleGL,WaveAngleSC)


! bottom u
        CALL GATHER_SWAN2GLOBAL(WaveUbottSW,WaveUbottGL)

        CALL DISTRIBUTE2SHORECIRC(WaveUbottGL,WaveUbottSC)

! wave fluxX
        CALL GATHER_SWAN2GLOBAL(WaveFluxXSW,WaveFluxXGL)

        CALL DISTRIBUTE2SHORECIRC(WaveFluxXGL,WaveFluxXSC)

! wave fluxY
        CALL GATHER_SWAN2GLOBAL(WaveFluxYSW,WaveFluxYGL)

        CALL DISTRIBUTE2SHORECIRC(WaveFluxYGL,WaveFluxYSC)

! meanperiod
        CALL GATHER_SWAN2GLOBAL(AvePeriodSW,AvePeriodGL)

        CALL DISTRIBUTE2SHORECIRC(AvePeriodGL,AvePeriodSC)

! peakperiod
        CALL GATHER_SWAN2GLOBAL(PeakPeriodSW,PeakPeriodGL)

        CALL DISTRIBUTE2SHORECIRC(PeakPeriodGL,PeakPeriodSC)

! wave dissi
        CALL GATHER_SWAN2GLOBAL(WaveDissSW,WaveDissGL)

        CALL DISTRIBUTE2SHORECIRC(WaveDissGL,WaveDissSC)

! wave breaking fraction
        CALL GATHER_SWAN2GLOBAL(WaveBrFraSW,WaveBrFraGL)

        CALL DISTRIBUTE2SHORECIRC(WaveBrFraGL,WaveBrFraSC)



! end parallel


! wave force
     IF(SHORECIRC_RUN)THEN
       IF(TIME.GE.WC_LAG)THEN 


 
        CALL GATHER_SWAN2GLOBAL(WaveFxSW,WaveFxGL)


        CALL GATHER_SWAN2GLOBAL(WaveFySW,WaveFyGL)


        WaveFxGL=WaveFxGL*MASK_WC_INTERACT
        WaveFyGL=WaveFyGL*MASK_WC_INTERACT

        CALL DISTRIBUTE2SHORECIRC(WaveFxGL,WaveFxSC)
        CALL DISTRIBUTE2SHORECIRC(WaveFyGL,WaveFySC)


       ENDIF
     ENDIF ! end shorecirc_run = .true.

END SUBROUTINE SWAN2SHORECIRC


SUBROUTINE DISTRIBUTE2SHORECIRC(PHIGLOB_IN,PHI)
     USE GLOBAL
     IMPLICIT NONE
     INTEGER :: l
     INTEGER,DIMENSION(NumberProcessor) :: npxs,npys
     REAL(SP),DIMENSION(NumberProcessor) :: xx
     REAL(SP),DIMENSION(MGlob,NGlob),INTENT(IN) :: PHIGLOB_IN
    REAL(SP),DIMENSION(MGlob,NGlob) :: PHIGLOB
     REAL(SP),DIMENSION(MGlob+2*Nghost,NGlob+2*Nghost) :: PHIGL_GHOST
     REAL(SP),DIMENSION(Mloc,Nloc),INTENT(OUT) :: PHI
     REAL(SP) :: TMP_L, TMP_R

    IF(myid.eq.0)THEN

     PHIGLOB=PHIGLOB_IN

! periodic
     IF(PERIODIC_X)THEN
       DO J=1,NGlob
        TMP_L=PHIGLOB(MGlob-Num_Transit+1,J)
        TMP_R=PHIGLOB(Num_Transit,J)
        DO I=1,Num_Transit
          PHIGLOB(I,J)=TMP_L+(TMP_R-TMP_L)*(I+Num_Transit-1)/(2.0*Num_Transit)
        ENDDO
        DO I=MGlob-Num_Transit,MGlob
          PHIGLOB(I,J)=TMP_L+(TMP_R-TMP_L)*(I-MGlob+Num_Transit)/(2.0*Num_Transit)
        ENDDO
       ENDDO
     ENDIF

     IF(PERIODIC_Y)THEN
       DO I=1,MGlob
        TMP_L=PHIGLOB(I,NGlob-Num_Transit+1)
        TMP_R=PHIGLOB(I, Num_Transit)
        DO J=1,Num_Transit
          PHIGLOB(I,J)=TMP_L+(TMP_R-TMP_L)*(J+Num_Transit-1)/(2.0*Num_Transit)
        ENDDO
        DO J=NGlob-Num_Transit,NGlob
          PHIGLOB(I,J)=TMP_L+(TMP_R-TMP_L)*(J-NGlob+Num_Transit)/(2.0*Num_Transit)
        ENDDO
       ENDDO
     ENDIF

! end periodic

     DO J=Nghost+1,NGlob+NGhost
     DO I=Nghost+1,MGlob+Nghost
        PHIGL_GHOST(I,J) = PHIGLOB(I-Nghost,J-Nghost)
     ENDDO
     ENDDO
! ghost cell
        DO I=Nghost+1,MGlob+Nghost
           DO J=1,Nghost
              PHIGL_GHOST(I,J)=PHIGL_GHOST(I,Nghost+1)
           ENDDO
           DO J=NGlob+Nghost+1,NGlob+2*Nghost
              PHIGL_GHOST(I,J)=PHIGL_GHOST(I,NGlob+Nghost)
           ENDDO
        ENDDO
        DO J=1,NGlob+2*Nghost
           DO I=1,Nghost
              PHIGL_GHOST(I,J)=PHIGL_GHOST(Nghost+1,J)
           ENDDO
           DO I=MGlob+Nghost+1,MGlob+2*Nghost
              PHIGL_GHOST(I,J)=PHIGL_GHOST(MGlob+Nghost,J)
           ENDDO
        ENDDO
     ENDIF ! end of myid=0

! distribute

     call MPI_Gather(npx,1,MPI_INTEGER,npxs,1,MPI_INTEGER,&
          0,MPI_COMM_WORLD,ier)
     call MPI_Gather(npy,1,MPI_INTEGER,npys,1,MPI_INTEGER,&
          0,MPI_COMM_WORLD,ier)

       DO J=1,Nloc
       DO I=1,Mloc   
        if (myid.eq.0) then
           do l=1,px*py
              xx(l) = PHIGL_GHOST(i+npxs(l)*(Iend-Ibeg+1),&
                   j+npys(l)*(Jend-Jbeg+1))
           enddo
        endif
        call MPI_Scatter(xx,1,MPI_SP,&
             PHI(i,j),1,MPI_SP,0,MPI_COMM_WORLD,ier)
       ENDDO
       ENDDO


END SUBROUTINE DISTRIBUTE2SHORECIRC




SUBROUTINE GATHER_SWAN2GLOBAL(PHI,PHIGLOB)

      USE PARAM
      USE GLOBAL
      USE SWCOMM3
      USE M_PARALL
      USE PASS
      IMPLICIT NONE
      INTEGER,DIMENSION(:),ALLOCATABLE :: MXCs,MYCs
      
      INTEGER :: l,IDIR

      REAL(SP),DIMENSION(NPROC) :: xx
      REAL(SP),DIMENSION(MXCGL,MYCGL),INTENT(OUT) :: PHIGLOB
      REAL(SP),DIMENSION(MXC,MYC),INTENT(IN) :: PHI

      IF(FIRST_CALL_SW2SC)THEN
        FIRST_CALL_SW2SC=.FALSE.
      IF(.NOT.ALLOCATED(MXFs)) ALLOCATE(MXFs(NPROC))
      IF(.NOT.ALLOCATED(MYFs)) ALLOCATE(MYFs(NPROC))
      IF(.NOT.ALLOCATED(IXTRIM_Ls)) ALLOCATE(IXTRIM_Ls(NPROC))
      IF(.NOT.ALLOCATED(IYTRIM_Ls)) ALLOCATE(IYTRIM_Ls(NPROC))
      IF(.NOT.ALLOCATED(IXTRIM_Rs)) ALLOCATE(IXTRIM_Rs(NPROC))
      IF(.NOT.ALLOCATED(IYTRIM_Rs)) ALLOCATE(IYTRIM_Rs(NPROC))


! find relation between different processes

      IF ( MXCGL.GT.MYCGL ) THEN
         IDIR = 2  ! cut in x direction
      ELSE
         IDIR = 1
      END IF 

       call MPI_Gather(MXF,1,MPI_INTEGER,MXFs,1,MPI_INTEGER,&
          0,MPI_COMM_WORLD,ier)
       call MPI_Gather(MYF,1,MPI_INTEGER,MYFs,1,MPI_INTEGER,&
          0,MPI_COMM_WORLD,ier)
      

        IXTRIM_L=0
        IYTRIM_L=0
        IXTRIM_R=0
        IYTRIM_R=0

        IF (IDIR.EQ.2)THEN
           IF(INODE.NE.1)IXTRIM_L=3
           IF(INODE.NE.NPROC)IXTRIM_R=3
        ELSE
           IF(INODE.NE.1)IYTRIM_L=3
           IF(INODE.NE.NPROC)IYTRIM_R=3
        ENDIF

 
        call MPI_Gather(IXTRIM_L,1,MPI_INTEGER,IXTRIM_Ls,1,MPI_INTEGER,&
          0,MPI_COMM_WORLD,ier)  
        call MPI_Gather(IXTRIM_R,1,MPI_INTEGER,IXTRIM_Rs,1,MPI_INTEGER,&
          0,MPI_COMM_WORLD,ier)  
        call MPI_Gather(IYTRIM_L,1,MPI_INTEGER,IYTRIM_Ls,1,MPI_INTEGER,&
          0,MPI_COMM_WORLD,ier)  
        call MPI_Gather(IYTRIM_R,1,MPI_INTEGER,IYTRIM_Rs,1,MPI_INTEGER,&
          0,MPI_COMM_WORLD,ier)      

! this coupling algorithm requires grid number (longer one) be divided by nproc
! I make a warning in case not satisfy the requirement

         ALLOCATE(MXCs(NPROC),MYCs(NPROC))
 
         IF (IDIR.EQ.2)THEN

           call MPI_Gather(MXC,1,MPI_INTEGER,MXCs,1,MPI_INTEGER,&
               0,MPI_COMM_WORLD,ier) 
           IF(INODE.EQ.1)THEN
             DO l=1,NPROC
              IF(MXCs(l)-IXTRIM_Rs(l)-IXTRIM_Ls(l) &
              .NE.MXCs(1)-IXTRIM_Rs(1)-IXTRIM_Ls(1))THEN
                WRITE(*,*)'Grid number(in X direction) can NOT be divided by NPROC'
                WRITE(*,*)'Please STOP the program and reset NPROC'
              ENDIF
!WRITE(*,*)l,MXCs(l)-IXTRIM_Rs(l)-IXTRIM_Ls(l),MXCs(l),IXTRIM_Rs(l),IXTRIM_Ls(l)
             ENDDO
           ENDIF

         ELSE

           call MPI_Gather(MYC,1,MPI_INTEGER,MYCs,1,MPI_INTEGER,&
               0,MPI_COMM_WORLD,ier) 
           IF(INODE.EQ.1)THEN
             DO l=1,NPROC
              IF(MYCs(l)-IYTRIM_Rs(l)-IYTRIM_Ls(l) &
              .NE.MYCs(1)-IYTRIM_Rs(1)-IYTRIM_Ls(1))THEN
               WRITE(*,*)'Grid number(in Y direction) can NOT be divided by NPROC'
               WRITE(*,*)'Please STOP the program and reset NPROC'
              ENDIF
!WRITE(*,*)l,MYCs(l)-IYTRIM_Rs(l)-IYTRIM_Ls(l)
             ENDDO
           ENDIF

         ENDIF
         DEALLOCATE(MXCs,MYCs)

       ENDIF ! end first_call


       do IX=IXTRIM_L+1,MXC-IXTRIM_R
       do IY=IYTRIM_L+1,MYC-IYTRIM_R
          call MPI_Gather(PHI(IX,IY),1,MPI_SP,&
             xx,1,MPI_SP,0,MPI_COMM_WORLD,ier)
       
          if (IY.eq.IYTRIM_L+1) call MPI_Barrier(MPI_COMM_WORLD,ier)

          if (INODE.eq.1) then

             do l=1,NPROC
               PHIGLOB(MXFs(l)+IX-1+IXTRIM_Ls(l),MYFs(l)+IY-1+IYTRIM_Ls(l))=xx(l)
             enddo
          endif

       enddo
       enddo

END SUBROUTINE GATHER_SWAN2GLOBAL



SUBROUTINE DISTRIBUTE2SWAN(PHIGLOB,PHI)
      USE SWAN_COMMON
      IMPLICIT NONE
      INTEGER :: l,IDIR,ier
      REAL(SP),DIMENSION(NPROC) :: xx
      REAL(SP),DIMENSION(MXCGL,MYCGL),INTENT(IN) :: PHIGLOB
      REAL(SP),DIMENSION(MXC,MYC),INTENT(OUT) :: PHI

! the following may be performed in DISTRIBUTE2SWAN

      IF(FIRST_CALL_SW2SC)THEN
        FIRST_CALL_SW2SC=.FALSE.
      IF(.NOT.ALLOCATED(MXFs)) ALLOCATE(MXFs(NPROC))
      IF(.NOT.ALLOCATED(MYFs)) ALLOCATE(MYFs(NPROC))
      IF(.NOT.ALLOCATED(IXTRIM_Ls)) ALLOCATE(IXTRIM_Ls(NPROC))
      IF(.NOT.ALLOCATED(IYTRIM_Ls)) ALLOCATE(IYTRIM_Ls(NPROC))
      IF(.NOT.ALLOCATED(IXTRIM_Rs)) ALLOCATE(IXTRIM_Rs(NPROC))
      IF(.NOT.ALLOCATED(IYTRIM_Rs)) ALLOCATE(IYTRIM_Rs(NPROC))


! find relation between different processes

      IF ( MXCGL.GT.MYCGL ) THEN
         IDIR = 2  ! cut in x direction
      ELSE
         IDIR = 1
      END IF 

       call MPI_Gather(MXF,1,MPI_INTEGER,MXFs,1,MPI_INTEGER,&
          0,MPI_COMM_WORLD,ier)
       call MPI_Gather(MYF,1,MPI_INTEGER,MYFs,1,MPI_INTEGER,&
          0,MPI_COMM_WORLD,ier)
      

        IXTRIM_L=0
        IYTRIM_L=0
        IXTRIM_R=0
        IYTRIM_R=0

        IF (IDIR.EQ.2)THEN
           IF(INODE.NE.1)IXTRIM_L=3
           IF(INODE.NE.NPROC)IXTRIM_R=3
        ELSE
           IF(INODE.NE.1)IYTRIM_L=3
           IF(INODE.NE.NPROC)IYTRIM_R=3
        ENDIF

 
        call MPI_Gather(IXTRIM_L,1,MPI_INTEGER,IXTRIM_Ls,1,MPI_INTEGER,&
          0,MPI_COMM_WORLD,ier)  
        call MPI_Gather(IXTRIM_R,1,MPI_INTEGER,IXTRIM_Rs,1,MPI_INTEGER,&
          0,MPI_COMM_WORLD,ier)  
        call MPI_Gather(IYTRIM_L,1,MPI_INTEGER,IYTRIM_Ls,1,MPI_INTEGER,&
          0,MPI_COMM_WORLD,ier)  
        call MPI_Gather(IYTRIM_R,1,MPI_INTEGER,IYTRIM_Rs,1,MPI_INTEGER,&
          0,MPI_COMM_WORLD,ier)      


       ENDIF ! end first_call

       DO IX=IXTRIM_L+1,MXC-IXTRIM_R
       DO IY=IYTRIM_L+1,MYC-IYTRIM_R 

        if (INODE.eq.1) THEN
           do l=1,NPROC
              xx(l) = PHIGLOB(IX+MXFs(l)+IXTRIM_Ls(l)-1,&
                   MYFs(l)+IY-1+IYTRIM_Ls(l))
           enddo
        endif
        call MPI_Scatter(xx,1,MPI_SP,&
             PHI(IX,IY),1,MPI_SP,0,MPI_COMM_WORLD,ier)
       ENDDO
       ENDDO

! ghost cells
       DO IY=IYTRIM_L+1,MYC-IYTRIM_R
        DO IX=1,IXTRIM_L
         PHI(IX,IY)=PHI(IXTRIM_L+1,IY)
        ENDDO
        DO IX=MXC-IXTRIM_R+1,MXC
         PHI(IX,IY)=PHI(MXC-IXTRIM_R,IY)
        ENDDO
       ENDDO

       DO IX=1,MXC
        DO IY=1,IYTRIM_L
          PHI(IX,IY)=PHI(IX,IYTRIM_L+1)
        ENDDO
        DO IY=MYC-IYTRIM_R+1,MYC
          PHI(IX,IY)=PHI(IX,MYC-IYTRIM_R)
        ENDDO
       ENDDO

END SUBROUTINE DISTRIBUTE2SWAN



