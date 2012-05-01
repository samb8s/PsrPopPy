c=============================================================================
      real function lfunc(lmin,m,lmax,seed)
      implicit none
      real lmin,lmax,m,c,loglmin,loglmax,logl,logn,test,n,nmax
      integer seed
      real psrran

      loglmin=log10(lmin)
      loglmax=log10(lmax)
      c=-1.0*loglmax*m
      nmax=10.0**(m*loglmin+c)

 1    logl=loglmin+(loglmax-loglmin)*psrran(seed)
      n=10.0**(m*logl+c)
      test=nmax*psrran(seed)
      if (test.gt.n) goto 1
      lfunc=10.0**logl
      end


c      real function lfunc(lmin,m,c,seed)
c      implicit none
c      real lmin,m,c,loglmin,loglmax,logl,logn,test,n,nmax
c      integer seed
c      real psrran
c
c      loglmin=log10(lmin)
c      loglmax=-1.0*c/m
c      nmax=10.0**(m*loglmin+c)
c
c 1    logl=loglmin+(loglmax-loglmin)*psrran(seed)
c      n=10.0**(m*logl+c)
c      test=nmax*psrran(seed)
c      if (test.gt.n) goto 1
c      lfunc=10.0**logl
c      end

