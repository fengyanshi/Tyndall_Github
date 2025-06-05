      PROGRAM COUPLER
      USE SWAN_COMMON, ONLY : IT,IT0,MTC, DT_SC,TFINC,TINIC,&
                          SWAN_RUN,SHORECIRC_RUN



      USE GLOBAL, ONLY : myid,Total_Time,TIME



! end parallel



      USE PARAM
      IMPLICIT NONE
      REAL(SP) :: TOTAL_SEC,DT_SWAN_SEC,SC_INNER_TIME=ZERO




      CALL SWINITMPI                                                     

      CALL SWAN_INITIALIZATION


      CALL SHORECIRC_INITIALIZATION


      TOTAL_SEC = TFINC-TINIC
      DT_SWAN_SEC = INT(TOTAL_SEC/MTC)
      Total_Time = TOTAL_SEC

      DO 500 IT = IT0, MTC

! pass uvz to swan
     IF(SWAN_RUN.AND.SHORECIRC_RUN)THEN

      CALL SHORECIRC2SWAN

     ENDIF

! call swan
     IF(SWAN_RUN)THEN
      CALL SWAN_CYCLE
     ENDIF

! pass variable from swan to shorecirc
! note that when shorecirc_run=.false., we still need this pass for 
! swan output of height etc, not radiation stresses
     IF(SWAN_RUN)THEN
      CALL SWAN2SHORECIRC
     ENDIF

! call shorecirc


!      SC_INNER_TIME = ZERO ! there's small residual, should correct it
      SC_INNER_TIME = MAX(ZERO,SC_INNER_TIME - DT_SWAN_SEC)
      DO WHILE (SC_INNER_TIME<DT_SWAN_SEC)



      IF(SHORECIRC_RUN)THEN
        CALL SHORECIRC_CYCLE
      ELSE

        DT_SC = DT_SWAN_SEC
        TIME=TIME+DT_SC

      ENDIF


        SC_INNER_TIME=SC_INNER_TIME+DT_SC

! dye tracer


! sediment module

! end sediment inside shorecirc

        CALL OUTPUT


! end non-stationary

      END DO




 500      CONTINUE


     if (myid.eq.0) WRITE(*,*)'Normal Termination!'
     if (myid.eq.0) WRITE(3,*)'Normal Termination!'


      CALL SWEXITMPI 



      END
                                    
SUBROUTINE OUTPUT
      USE SWAN_COMMON
      USE GLOBAL
      IMPLICIT NONE

     SCREEN_COUNT=SCREEN_COUNT+DT_SC

     IF(SCREEN_COUNT>=SCREEN_INTV)THEN
      SCREEN_COUNT=SCREEN_COUNT-SCREEN_INTV
      CALL STATISTICS
     ENDIF

! stations
      IF(NumberStations>0)THEN
      PLOT_COUNT_STATION=PLOT_COUNT_STATION+DT_SC
      IF(PLOT_COUNT_STATION>=PLOT_INTV_STATION)THEN
       PLOT_COUNT_STATION=PLOT_COUNT_STATION-PLOT_INTV_STATION
       CALL STATIONS
      ENDIF
      ENDIF
! preview
      PLOT_COUNT=PLOT_COUNT+DT_SC
      IF(PLOT_COUNT>=PLOT_INTV)THEN
       PLOT_COUNT=PLOT_COUNT-PLOT_INTV
       CALL PREVIEW
      ENDIF

END SUBROUTINE OUTPUT



  

