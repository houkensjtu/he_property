      program rootFindCompare
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CONSTA/G,R,TAO
      
      P  = 0.3d6
      T  = 3.0d0
      RP = 0.0
      RP0 = 0.0
      
      do i=1,600
        T = T * 1.1
        call cpu_time(t1)
        do j=1,100000
CC      P  = 1500000.0
          RP = RP0
          call RPCALCUL(T,RP,P)
CC      print *,T,RP,P
CC      print *,"Collapsed CPU time: ", t2-t1
        
CC      print '("P  = ",f6.3 ,", T  = ", f8.3, ", RP = ", f9.4)',PP,T,RP
CC      write (1,'(f6.3, ",", f8.3, ",", f9.4)') PP,T,RP
        end do
        call cpu_time(t2)
        
        print *,t,",",rp,",",t2-t1
        
      end do
      end program rootFindCompare
