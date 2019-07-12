      program densityOutput
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CONSTA/G,R,TAO
      
      open(1, file='he.dat', status='replace')
      
      P  = 0.3d6
      T  = 5.0d0
      RP = 10.0
      do i=1,2
      T = 5.0d0
      P=P+0.01d6
        do j=1,100000
        T = T + 0.01
        PP = P/1.0d6 
CC      call cpu_time(t1)
CC      P  = 1500000.0
        call RPCALCUL(T,RP,P)
CC      print *,T,RP,P
CC      call cpu_time(t2)
CC      print *,"Collapsed CPU time: ", t2-t1
        
CC      print '("P  = ",f6.3 ,", T  = ", f8.3, ", RP = ", f9.4)',PP,T,RP
        write (1,'(f6.3, ",", f8.3, ",", f9.4)') PP,T,RP
        end do
        print *,"Current i = ",i
      end do
      end program densityOutput
