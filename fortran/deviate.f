ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function deviate(seed,binned,maxbin,xmin,xmax,nbins)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     uses the accept/reject technique to return a deviate from a
c     binned distribution passed down as binned(nbins) over the
c     ranges xmin to xmax, where x = xmin + (xmax-xmin)*i/nbins
c     the maximum bin value is passed down as maxbin to save 
c     repeated computations of it within this routine.
c     
      implicit none
      integer seed,nbins,binned(*),maxbin,i
      real x,y,psrran,xmin,xmax,delta
      logical accepted

      accepted=.false.
      delta=xmax-xmin
      do while (.not.accepted)
c
c        generate a deviate within the reference
c        function, a top-hat with the same height as maxbin
c
         x=xmin+delta*psrran(seed)
         y=psrran(seed)*real(maxbin)
c
c        test to see whether this falls under the binned
c        distribution and accept it if it does
c
         i=(x-xmin)*real(nbins)/delta
         accepted = (y.le.real(binned(i)))
      enddo
c
c     return the deviate drawn from the chosen distribution
c
      deviate=x
      end
