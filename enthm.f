CC
CC      ************************************************************
CC      **              MODULE TO CALCULATE ENTHALPY              **
CC      **      WHEN T P AND THE INITIAL RP ARE KNOWN             **
CC      **     UNITS                                              **
CC      **         P,Pa; RP, kg/m3; T,K.                          **
CC      ************************************************************
CC
        DOUBLE PRECISION FUNCTION ENTHM (P,D,T)
        
	  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DATA S0, H0 /1964.2D0, 15273.D0/
        DATA R/8.31434D-3/
      
        DD1=D/4.0026D0
        TT=T
        PP=P/1.0D6
        CALL PROPS(SD,DD1,TT,2)
        CALL PROPS(UD,DD1,TT,3)
        DD0=0.D0
        CALL PROPS(S0,DD0,TT,2)
        CALL PROPS(U0,DD0,TT,3)
        print *,"In ENTHM CPI = ",CPI(TT,3)
        print *,"In ENTHM SD  = ",SD
        print *,"In ENTHM UD = ",UD
        print *,"In ENTHM S0 = ",S0
        print *,"In ENTHM U0 = ",U0
        ENTHM=TT*(SD-S0)*1000.D0+(UD-U0)*1000.D0+CPI(TT,3)
     &   +(PP/DD1-R*TT)*1.D+3
        ENTHM=ENTHM/4.0026D0*1.0D3
	  ENTHM=ENTHM+H0
	  RETURN
          END

CC
CC      ************************************************************
CC      **   Calculation of the specific heat capacity of helium  **
CC      **    Cp J/(MOL.K)    Temperature:  4-300k                 **
CC      ************************************************************
CC
        DOUBLE PRECISION FUNCTION CPI(T,K)
        
	  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        
        DIMENSION G(11),GS(11)
        DATA RR/8.31434D0/
        DATA G /  .00D+00, .00D+00, .00D+00, .25D+01, .00D+00,
     A  .00D+00, .00D+00, .00D+00, .25D+01, .00D+00, .00D+00/
      
        DO 1 I=1,11
    1   GS(I)=G(I)
        U=G(9)/T
        EU=DEXP(U)
        TS=1.D0/T**4
        GO TO (20,40,55),K
   20   CPI=G(8)*U*U*EU/(EU-1.D0)**2
        DO 25 I=1,7
        TS=TS*T
   25   CPI=CPI+G(I)*TS
        CPI=CPI*RR
        GO TO 60
   40   CPI=G(8)*(U/(EU-1.D0)-DLOG(1.D0-1.D0/EU))
     1   -G(1)*TS*T/3.D0-G(2)*TS*T*T/2.D0-G(3)/T+G(4)*DLOG(T)+G(5)*T
     2   +G(6)*T*T/2.D0+G(7)*T**3/3.D0
        CPI=CPI*RR+G(11)
        GO TO 60
   55   CPI=G(8)*U*T/(EU-1.D0)-G(1)/(2.D0*T*T)-G(2)/T+G(3)*DLOG(T)
     1   +G(4)*T+G(5)*T*T/2.D0+G(6)*T**3/3.D0+G(7)*T**4/4.D0
        CPI=CPI*RR+G(10)
   60   DO 61 I=1,11
   61   G(I)=GS(I)

        RETURN
        END


C THE 32 TERM EQUATION OF STATE, INPUT IS DENSITY(MOL/L),
C TEMPERATURE(K), OUTPUT (PP) IS PRESSURE(MPA),OR DP/DD IN
C L-MPA/MOL OR DP/DT MPA/K OR S,H,OR CV AT ONE LIMIT OF
C INTEGRATION
        SUBROUTINE PROPS(PP,DD,TT,K)  
        
	  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    
        DIMENSION X(33),G(32),B(32)
        DATA GAMMA/-.33033259D-02/
        DATA M/32/       
	  DATA R/0.00831434D0/      
      DATA G /               .4558980227431D-04,  .1260692007853D-02,
     A -.7139657549318D-02,  .9728903861441D-02, -.1589302471562D-01,
     B  .1454229259623D-05, -.4708238429298D-04,  .1132915223587D-02,
     C  .2410763742104D-02, -.5093547838381D-08,  .2699726927900D-05,
     D -.3954146691114D-04,  .1551961438127D-08,  .1050712335785D-07,
     E -.5501158366750D-07, -.1037673478521D-09,  .6446881346448D-12,
     F  .3298960057071D-10, -.3555585738784D-12, -.6885401367690D-02,
     G  .9166109232806D-02, -.6544314242937D-05, -.3315398880031D-04,
     H -.2067693644676D-07,  .3850153114958D-07, -.1399040626999D-10,
     I -.1888462892389D-11, -.4595138561035D-14,  .6872567403738D-14,
     J -.6097223119177D-18, -.7636186157005D-17,  .3848665703556D-17/
       
        D=DD
        P=PP
        T=TT
        GM=GAMMA
        D2=D*D
        D3=D2*D
        D4=D3*D
        D5=D4*D
        D6=D5*D
        D7=D6*D
        D8=D7*D
        D9=D8*D
        D10=D9*D
        D11=D10*D
        D12=D11*D
        D13=D12*D
        TS=DSQRT(T)
        T2=T*T
        T3=T2*T
        T4=T3*T
        T5=T4*T
        F=DEXP(GM*D2)
        GO TO (100,200,300,400,500),K

C     ENTRY DPDD
  400 F1=2.D0*F*GM*D
      F21=3.D0*F*D2 +F1*D3
      F22=5.D0*F*D4 +F1*D5
      F23=7.D0*F*D6 +F1*D7
      F24=9.D0*F*D8 +F1*D9
      F25=11.D0*F*D10+F1*D11
      F26=13.D0*F*D12+F1*D13
      B( 1)=2.D0*D*T
      B( 2)=2.D0*D*TS
      B( 3)=2.D0*D
      B( 4)=2.D0*D/T
      B( 5)=2.D0*D/T2
      B( 6)=3.D0*D2*T
      B( 7)=3.D0*D2
      B( 8)=3.D0*D2/T
      B( 9)=3.D0*D2/T2
      B(10)=4.D0*D3*T
      B(11)=4.D0*D3
      B(12)=4.D0*D3/T
      B(13)=5.D0*D4
      B(14)=6.D0*D5/T
      B(15)=6.D0*D5/T2
      B(16)=7.D0*D6/T
      B(17)=8.D0*D7/T
      B(18)=8.D0*D7/T2
      B(19)=9.D0*D8/T2
      B(20)=F21/T2
      B(21)=F21/T3
      B(22)=F22/T2
      B(23)=F22/T4
      B(24)=F23/T2
      B(25)=F23/T3
      B(26)=F24/T2
      B(27)=F24/T4
      B(28)=F25/T2
      B(29)=F25/T3
      B(30)=F26/T2
      B(31)=F26/T3
      B(32)=F26/T4
      P=0
      DO 201 I=1,M
  201 P=P+B(I)*G(I)
      P=P+R*T
      PP=P
      RETURN
C     ENTRY DPDT
  500 X( 1)=D2
      X( 2)=D2/(2.D0*TS)
      X( 3)=0.D0
      X( 4)=-D2/T2
      X( 5)=-2.D0*D2/T3
      X( 6)=D3
      X( 7)=0.D0
      X( 8)=-D3/T2
      X( 9)=-2.D0*D3/T3
      X(10)=D4
      X(11)=0.D0
      X(12)=-D4/T2
      X(13)=0.D0
      X(14)=-D6/T2
      X(15)=-2.D0*D6/T3
      X(16)=-D7/T2
      X(17)=-D8/T2
      X(18)=-2.D0*D8/T3
      X(19)=-2.D0*D9/T3
      X(20)=-2.D0*D3*F/T3
      X(21)=-3.D0*D3*F/T4
      X(22)=-2.D0*D5*F/T3
      X(23)=-4.D0*D5*F/T5
      X(24)=-2.D0*D7*F/T3
      X(25)=-3.D0*D7*F/T4
      X(26)=-2.D0*D9*F/T3
      X(27)=-4.D0*D9*F/T5
      X(28)=-2.D0*D11*F/T3
      X(29)=-3.D0*D11*F/T4
      X(30)=-2.D0*D13*F/T3
      X(31)=-3.D0*D13*F/T4
      X(32)=-4.D0*D13*F/T5
      P=0
      DO 301 I=1,M
  301 P=P+G(I)*X(I)
      PP=P+R*D
      RETURN
C     ENTRY TDSDT
C     TEMP. TIMES THE PARTIAL OF
C     ENTROPY WITH RESPECT TO TEMP.
C     CV=CV0+(TDSDN(/)-TDSDN(D))*1000.

100     G1=F/(2.D0*GM)
        G2=(F*D2-2.D0*G1)/(2.D0*GM)
        G3=(F*D4-4.D0*G2)/(2.D0*GM)
        G4=(F*D6-6.D0*G3)/(2.D0*GM)
        G5=(F*D8-8.D0*G4)/(2.D0*GM)
        G6=(F*D10-10.D0*G5)/(2.D0*GM)
        X(1)=0.D0
        X( 2)=-D/(4.D0*TS)
        X(3)=0.D0
        X( 4)=2.D0*D/T2
        X( 5)=6.D0*D/T3
        X(6)=0.D0
        X(7)=0.D0
        X( 8)=D2/T2
        X( 9)=3.D0*D2/T3
        X(10)=0.D0
        X(11)=0.D0
        X(12)=(2.D0*D3)/(3.D0*T2)
        X(13)=0.D0
        X(14)=(2.D0*D5)/(5.D0*T2)
        X(15)=(6.D0*D5)/(5.D0*T3)
        X(16)=D6/(3.D0*T2)
        X(17)=(2.D0*D7)/(7.D0*T2)
        X(18)=(6.D0*D7)/(7.D0*T3)
        X(19)=(3.D0*D8)/(4.D0*T3)
        X(20)=6.D0*G1/T3
        X(21)=12.D0*G1/T4
        X(22)=6.D0*G2/T3
        X(23)=20.D0*G2/T5
        X(24)=6.D0*G3/T3
        X(25)=12.D0*G3/T4
        X(26)=6.D0*G4/T3
        X(27)=20.D0*G4/T5
        X(28)=6.D0*G5/T3
        X(29)=12.D0*G5/T4
        X(30)=6.D0*G6/T3
        X(31)=12.D0*G6/T4
        X(32)=20.D0*G6/T5
        P=0
        DO 601 I=1,M
  601   P=P+G(I)*X(I)
        PP=P
        RETURN
C     ENTRY DSDN
C     PARTIAL OF ENTROPY WITH
C     RESPECT TO THE G COEFFICIENTS
C     S=S0-R*LOGF(D*R*T/P0)+(DSDG(D)-DSDG(0))*1000.DO +CPOS(T)
  200   G1=F/(2.D0*GM)
        G2=(F*D2-2.D0*G1)/(2.D0*GM)
        G3=(F*D4-4.D0*G2)/(2.D0*GM)
        G4=(F*D6-6.D0*G3)/(2.D0*GM)
        G5=(F*D8-8.D0*G4)/(2.D0*GM)
        G6=(F*D10-10.D0*G5)/(2.D0*GM)
        X( 1)=-D
        X( 2)=-D/(2.D0*TS)
        X( 3)=0.D0
        X( 4)=+D/T2
        X( 5)=2.D0*D/T3
        X( 6)=-D2/2.D0
        X( 7)=0.D0
        X( 8)=D2/(2.D0*T2)
        X( 9)=D2/T3
        X(10)=-D3/3.D0
        X(11)=0.D0
        X(12)=D3/(3.D0*T2)
        X(13)=0.D0
        X(14)=D5/(5.D0*T2)
        X(15)= 2.D0*D5/(5.D0*T3)
        X(16)=D6/(6.D0*T2)
        X(17)=D7/(7.D0*T2)
        X(18)=2.D0*D7/(7.D0*T3)
        X(19)=D8/(4.D0*T3)
        X(20)=2.D0*G1/T3
        X(21)=3.D0*G1/T4
        X(22)=2.D0*G2/T3
        X(23)=4.D0*G2/T5
        X(24)=2.D0*G3/T3
        X(25)=3.D0*G3/T4
        X(26)=2.D0*G4/T3
        X(27)=4.D0*G4/T5
        X(28)=2.D0*G5/T3
        X(29)=3.D0*G5/T4
        X(30)=2.D0*G6/T3
        X(31)=3.D0*G6/T4
        X(32)=4.D0*G6/T5
        P=0
        DO 401 I=1,M
  401   P=P+G(I)*X(I)
        PP=P
        RETURN
C     ENTRY DUDN
C     TERMS NEEDED FOR ENTHALPY CALCULATION
C     H=H0+(T*DSDG(D)-DSDG(0))*1000.+(DUDG(D-DUDG(0))*1000.+CPOH(T)
C     +(P/D-R*T)*1000.
  300   G1=F/(2.D0*GM)
        G2=(F*D2-2.D0*G1)/(2.D0*GM)
        G3=(F*D4-4.D0*G2)/(2.D0*GM)
        G4=(F*D6-6.D0*G3)/(2.D0*GM)
        G5=(F*D8-8.D0*G4)/(2.D0*GM)
        G6=(F*D10-10.D0*G5)/(2.D0*GM)
        X( 1)=D*T
        X( 2)=D*TS
        X( 3)=D
        X( 4)=D/T
        X( 5)=D/T2
        X( 6)=D2*T/2.D0
        X( 7)=D2/2.D0
        X( 8)=D2/(2.D0*T)
        X( 9)=D2/(2.D0*T2)
        X(10)=D3*T/3.D0
        X(11)=D3/3.D0
        X(12)=D3/(3.D0*T)
        X(13)=D4/4.D0
        X(14)=D5/(5.D0*T)
        X(15)=D5/(5.D0*T2)
        X(16)=D6/(6.D0*T)
        X(17)=D7/(7.D0*T)
        X(18)=D7/(7.D0*T2)
        X(19)=D8/(8.D0*T2)
        X(20)=G1/T2
        X(21)=G1/T3
        X(22)=G2/T2
        X(23)=G2/T4
        X(24)=G3/T2
        X(25)=G3/T3
        X(26)=G4/T2
        X(27)=G4/T4
        X(28)=G5/T2
        X(29)=G5/T3
        X(30)=G6/T2
        X(31)=G6/T3
        X(32)=G6/T4
        P=0
        DO 501 I=1,M
  501   P=P+G(I)*X(I)
        PP=P
        RETURN
        END
      
