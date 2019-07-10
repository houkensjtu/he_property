      program helium
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CONSTA/G,R,TAO

      P  = 10
      T  = 5.0
      RP = 10.0

      do j=1,5
      P = P * 1.5
      call cpu_time(t1)
      do i=1,3
CC    P  = 1500000.0
      call RPCALCUL(T,RP,P)
CC    print *,T,RP,P
      end do
      call cpu_time(t2)

      print *," "
      print *,"Collapsed CPU time: ", t2-t1
      print *,"P  = ", p
      print *,"RP = ", rp
      print *," "
      end do

      end program helium
