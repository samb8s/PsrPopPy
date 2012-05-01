c==============================================================================
        real function draw1d(xmin,xmax,fx,nbins,seed)
c==============================================================================
c
c     Function to draw a random number from the binned probability
c     distribution fx(nbins) which is a function of the variable
c     x(nbins) where seed is for the random number generator
c
c     How it works is to work out the total area under all the bins
c     - assuming they are equally spaced this is proportional to the
c       sum of fx. Then it generates a random number and finds the 
c       nearest x bin into which this would fall. To avoid quantization
c       the final x value is drawn from anywhere inside the xbin
c

      implicit none
      real xmin, xmax, fx(*)
      integer nbins, seed

      integer i
      real psum, ptot, ranno, psrran

      ptot=0.0
      do i=1,nbins
         ptot=ptot+fx(i)
      enddo
c
c     NB this assumes equal bin widths!
c
      ranno=psrran(seed)
      psum=0.0
      do i=1,nbins
        psum=psum+fx(i)
        if (ranno.le.psum/ptot) then
          ranno=psrran(seed)
          draw1d=xmin+(float(i-1)+ranno)*(xmax-xmin)/float(nbins)
          return
        endif
      enddo
      end

