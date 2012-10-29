c==============================================================================
      real function drawgamma(pmin,pmax,a,m,seed)
c==============================================================================
c
c            utility to draw periods from a gamma function of the form
c            f(P) = (P/m)^(a-1) exp(-P/m) 
c            where a and m are free parameters
c            drawgamma uses the rejection method (see Press et al) to
c            find the appropriate deviate given a choice of a and m
c
      real p,gammafunc,pmin,pmax,a,m,gmax,pgmax,x,y,psrran
      integer seed
c
c     figure out the size of a box which encompasses the gamma distribution
c
      pgmax=(a-1.0)*m        ! this is the peak of the PDF
      gmax=gammafunc(pgmax,a,m)  ! this is the maximum y value of our box
c
c     draw uniform deviates from box until one lies under the gamma curve
c
      do while (.true.)
         x=pmin+(pmax-pmin)*psrran(seed)
         y=gmax*psrran(seed)
         if (y.le.gammafunc(x,a,m)) then
            drawgamma=x
            return
         endif
      enddo
      end
