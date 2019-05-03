SUBROUTINE RPCALCUL(T,RP,P)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G(32)
      COMMON /CONSTA/G,R,TAO
      DATA TAO/0.0033033259D0/
      DATA R/0.00831434D0/
      DATA G/               .4558980227431D-04,  .1260692007853D-02,
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

       F(RP)=PP-(RP*A1+RP*RP*A2+RP**3*A3+RP**4*A4+RP**5*A5+RP**6*A6
     &        +RP**7*A7+RP**8*A8+RP**9*A9
     &        +(RP**3*B3+RP**5*B5+RP**7*B7+RP**9*B9+RP**11*B11
     &        +RP**13*B13)*DEXP(-TAO*RP*RP))
       FP(RP)=-(A1+2.0D0*RP*A2+3.0D0*RP**2*A3+4.0D0*RP**3*A4
     &         +5.0D0*RP**4*A5+6.0D0*RP**5*A6+7.0D0*RP**6*A7
     &         +8.0D0*RP**7*A8+9.0D0*RP**8*A9+(3.0D0*RP**2*B3
     &         +5.0D0*RP**4*B5+7.0D0*RP**6*B7+9.0D0*RP**8*B9
     &         +11.0D0*RP**10*B11+13.0D0*RP**12*B13)*DEXP(-TAO*RP*RP)
     &        +(RP**3*B3+RP**5*B5+RP**7*B7+RP**9*B9+RP**11*B11
     &        +RP**13*B13)*DEXP(-TAO*RP*RP)*(-2.0D0*TAO*RP))

        PP=P/1.0D6
        A1=R*T
        A2=G(1)*T+G(2)*DSQRT(T)+G(3)+G(4)/T+G(5)/T/T
        A3=G(6)*T+G(7)+G(8)/T+G(9)/T/T
        A4=G(10)*T+G(11)+G(12)/T
        A5=G(13)
        A6=G(14)/T+G(15)/T/T
        A7=G(16)/T
        A8=G(17)/T+G(18)/T/T
        A9=G(19)/T/T
        
        B3=G(20)/T/T+G(21)/(T**3)
        B5=G(22)/T/T+G(23)/(T**4)
        B7=G(24)/T/T+G(25)/(T**3)
        B9=G(26)/T/T+G(27)/(T**4)
        B11=G(28)/T/T+G(29)/(T**3)
        B13=G(30)/T/T+G(31)/(T**3)+G(32)/(T**4)

        EPS=1.0D-6
        IF(T.GT.8.0D0) THEN
        RP0=RP/4.0026D0
        ELSE
        RP0=160.0D0/4.0026D0
        ENDIF
5       F0=F(RP0)
10      FP0=FP(RP0)
        IF(DABS(FP0).LE.EPS) THEN
           RP0=RP0+0.00001D0*RP0
        GOTO 5
        ELSE
           RP1=RP0-F0/FP0
        IF (RP1.GT.0.0D0) GOTO 20
            RP1=RP0+0.1D0*RP0
20         F1=F(RP1)
           IF(DABS(F1).LT.EPS) THEN
                RP=RP1
           ELSE
                 IF(DABS((RP1-RP0)/RP0).LT.1.0D-6) THEN
                   RP0=RP0+0.000001D0*RP0
                   GOTO 5
               ELSE
                   RP0=RP1
                   F0=F1
                   GOTO 10
                ENDIF
             ENDIF
        ENDIF
        RP=RP*4.0026D0

        RETURN
        END
