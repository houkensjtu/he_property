C     This program is used to test enthalpy calculation and density 
      program enthalpyChk
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CONSTA/G,R,TAO
      
      P  = 0.5d6
      T  = 4.0d0
      RP = 0.0
      RP0 = 0.0
      h = 0.0      
      do i=1,270
        T = T + 1.0
        call cpu_time(t1)
c        do j=1,2
          RP = RP0
          call RPCALCUL(T,RP,P)
          h = ENTHM(P,RP,T)
c        end do
        call cpu_time(t2)
        
        print *,t,", density =",rp,", entalpy =",h,", time =",t2-t1

        
      end do
      end program enthalpyChk
