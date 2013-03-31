C common for the numerical recipes routine to store the path (note this should
C not be changed without changing ODEINT)
      integer maxpts
      parameter (maxpts=2000)
      integer kmax, kount
      real dxsav, xplot(maxpts),yplot(10,maxpts)
      common /path/ kmax,kount,dxsav,xplot,yplot
