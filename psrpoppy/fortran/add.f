      real function addtwo(a, b)
        integer, intent(in) :: a, b
        addtwo = a+b
      end function


      real function ykr(seed)
c==============================================================================
c
c     uses ykarea to draw a random deviate from Yusifov & Kukuk's 
c     radial distribution. a=1.64,b=4.01,r1=0.55
c
      integer, intent(in) :: seed
      real ykr0
c      ykr=ykr0(seed,1.64,4.01,0.55)
      ykr= real(seed)
      end function   
