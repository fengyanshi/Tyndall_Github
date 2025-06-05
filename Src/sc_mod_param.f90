MODULE PARAM

       use mpi

       IMPLICIT NONE


       INTEGER, PARAMETER::SP=SELECTED_REAL_KIND(6,30)

       INTEGER, PARAMETER::MPI_SP=MPI_REAL


! define parameters
       REAL(SP), PARAMETER::pi=3.141592653
       REAL(SP), PARAMETER::R_earth = 6371000.0_SP
       REAL(SP), PARAMETER::SMALL=0.000001_SP
       REAL(SP), PARAMETER::LARGE=100000.0_SP
       REAL(SP), PARAMETER:: grav=9.81_SP
       REAL(SP), PARAMETER:: zero = 0.0_SP
       REAL(SP), PARAMETER:: DEG2RAD = 0.0175_SP
       REAL(SP), PARAMETER:: RHO_AW = 0.0012041_SP
       REAL(SP), PARAMETER:: RHO_W = 1025.0_SP

! some global variables
       INTEGER :: I,J,K
       INTEGER :: itmp1,itmp2,itmp3,itmp4,itmp5
       REAL(SP):: tmp1,tmp2,tmp3,tmp4,tmp5

END MODULE PARAM

