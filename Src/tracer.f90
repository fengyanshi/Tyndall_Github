! end tracking



SUBROUTINE ADVECTION_DIFFUSION(Con,DT)
      USE PARAM
      USE GLOBAL,ONLY : U,V,Mloc,Nloc,H,Mloc1,Nloc1,MASK,Ibeg,Iend,Jbeg,Jend, &
                         Jaco,Ubott,Vbott,tmp4preview,MinDepthFrc, &
                         L11,L12,L22,nu_total 
      IMPLICIT NONE
      REAL(SP),DIMENSION(Mloc,Nloc),INTENT(INOUT) :: Con
      REAL(SP), INTENT(IN) :: DT
      REAL(SP),DIMENSION(Mloc,Nloc) :: Q1,Q2
      REAL(SP),DIMENSION(Mloc1,Nloc) :: Q1L,Q1R,Q1c
      REAL(SP),DIMENSION(Mloc,Nloc1) :: Q2L,Q2R,Q2c
      REAL(SP),DIMENSION(Mloc,Nloc) :: DelxQ1
      REAL(SP),DIMENSION(Mloc,Nloc) :: DelyQ2
      REAL(SP) :: HL11_1,HL11_2,HL12_1,HL12_2,HL22_1,HL22_2, &
                  C_11,C_12,C_22,C_1,C_2,HL11,HL12,HL22

      DO J=1,Nloc
      DO I=1,Mloc
        Q1(I,J)=(U(I,J)*L11(I,J)+V(I,J)*L12(I,J))*Jaco(I,J)*Con(I,J)
        Q2(I,J)=(U(I,J)*L12(I,J)+V(I,J)*L22(I,J))*Jaco(I,J)*Con(I,J)
      ENDDO
      ENDDO

      CALL DelxFun_minmod(1.0_SP,Mloc,Nloc,Q1,DelxQ1)
      CALL DelyFun_minmod(1.0_SP,Mloc,Nloc,Q2,DelyQ2)
     
      CALL CONSTRUCT_X(Mloc,Mloc1,Nloc,1.0_SP,Q1,DelxQ1,Q1L,Q1R)
      CALL CONSTRUCT_Y(Mloc,Nloc,Nloc1,1.0_SP,Q2,DelyQ2,Q2L,Q2R)

! 1st order TVD
       DO J=1,Nloc
       DO I=1,Mloc1
        IF((Q1L(I,J)+Q1R(I,J)).GE.ZERO) THEN
          Q1c(I,J)=Q1L(I,J)
        ELSE
          Q1c(I,J)=Q1R(I,J)
        ENDIF
       ENDDO
       ENDDO

       DO J=1,Nloc1
       DO I=1,Mloc
        IF((Q2R(I,J)+Q2L(I,J)).LE.ZERO) THEN
          Q2c(I,J)=Q2R(I,J)
        ELSE
          Q2c(I,J)=Q2L(I,J)
        ENDIF        
       ENDDO
       ENDDO

! deal with masks

      DO J=1,Nloc
      DO I=1,Mloc
         IF(MASK(I,J)<1)THEN
           Q1c(I,J)=ZERO
           Q1c(I+1,J)=ZERO
           Q2c(I,J+1)=ZERO
           Q2c(I,J)=ZERO
         ENDIF
      ENDDO
      ENDDO

      DO J=Jbeg,Jend
      DO I=Ibeg,Iend

       IF(MASK(I,J)>0)THEN
        Con(I,J)=Con(I,J)-1.0_SP/MAX(SMALL,Jaco(I,J))*(Q1c(I+1,J)-Q1c(I,J)  &
                               +Q2c(I,J+1)-Q2c(I,J))*DT                               
       ELSE
        Con(I,J)=ZERO
       ENDIF

      ENDDO
      ENDDO

! diffusion

      DO J=Jbeg,Jend
      DO I=Ibeg,Iend

       IF(MASK(I,J)>0)THEN
         HL11_1=0.5_SP*(H(I+1,J)*nu_total(I+1,J)*L11(I+1,J)-H(I-1,J)*nu_total(I-1,J)*L11(I-1,J))
         HL12_1=0.5_SP*(H(I+1,J)*nu_total(I+1,J)*L12(I+1,J)-H(I-1,J)*nu_total(I-1,J)*L12(I-1,J))
         HL22_1=0.5_SP*(H(I+1,J)*nu_total(I+1,J)*L22(I+1,J)-H(I-1,J)*nu_total(I-1,J)*L22(I-1,J))
         HL11_2=0.5_SP*(H(I,J+1)*nu_total(I,J+1)*L11(I,J+1)-H(I,J-1)*nu_total(I,J-1)*L11(I,J-1))
         HL12_2=0.5_SP*(H(I,J+1)*nu_total(I,J+1)*L12(I,J+1)-H(I,J-1)*nu_total(I,J-1)*L12(I,J-1))
         HL22_2=0.5_SP*(H(I,J+1)*nu_total(I,J+1)*L22(I,J+1)-H(I,J-1)*nu_total(I,J-1)*L22(I,J-1))
         HL11=H(I,J)*L11(I,J)*nu_total(I,J)
         HL12=H(I,J)*L12(I,J)*nu_total(I,J)
         HL22=H(I,J)*L22(I,J)*nu_total(I,J)

         C_1=(Con(I+1,J)-Con(I-1,J))*0.5_SP/MAX(SMALL,H(I,J))
         C_2=(Con(I,J+1)-Con(I,J-1))*0.5_SP/MAX(SMALL,H(I,J))
         C_11=Con(I+1,J)-2.0_SP*Con(I,J)+Con(I-1,J)/MAX(SMALL,H(I,J))
         C_22=Con(I,J+1)-2.0_SP*Con(I,J)+Con(I,J-1)/MAX(SMALL,H(I,J))
         C_12=0.25_SP*(Con(I+1,J+1)-Con(I-1,J+1)  &
            -Con(I+1,J-1)+Con(I-1,J-1))/MAX(SMALL,H(I,J))
         Con(I,J)=Con(I,J) + ( &
           HL11*(HL11*C_11+C_1*HL11_1+HL12*C_12+C_2*HL12_1) &
          +HL12*(HL11*C_12+C_1*HL11_2+HL12*C_22+C_2*HL12_2) &
          +HL12*(HL12*C_11+C_1*HL12_1+HL22*C_12+C_2*HL22_1) &
          +HL22*(HL12*C_12+C_1*HL12_2+HL22*C_22+C_2*HL22_2) )
       ELSE
        Con(I,J)=ZERO
       ENDIF

      ENDDO
      ENDDO


END SUBROUTINE ADVECTION_DIFFUSION

