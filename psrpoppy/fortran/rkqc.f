c-----------------------------------------------------------------
      subroutine rkqc(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
c-----------------------------------------------------------------
      parameter (nmax=10,fcor=.0666666667,
     *    one=1.,safety=0.9,errcon=6.e-4)
      external derivs
      dimension y(n),dydx(n),yscal(n),ytemp(nmax),ysav(nmax),dysav(nmax)
      pgrow=-0.20
      pshrnk=-0.25
      xsav=x
      do 11 i=1,n
        ysav(i)=y(i)
        dysav(i)=dydx(i)
11    continue
      h=htry
1     hh=0.5*h
      call rk4(ysav,dysav,n,xsav,hh,ytemp,derivs)
      x=xsav+hh
      call derivs(x,ytemp,dydx)
      call rk4(ytemp,dydx,n,x,hh,y,derivs)
      x=xsav+h
c      if(x.eq.xsav)pause 'stepsize not significant in rkqc.'
      call rk4(ysav,dysav,n,xsav,h,ytemp,derivs)
      errmax=0.
      do 12 i=1,n
        ytemp(i)=y(i)-ytemp(i)
        errmax=max(errmax,abs(ytemp(i)/yscal(i)))
12    continue
      errmax=errmax/eps
      if(errmax.gt.one) then
        h=safety*h*(errmax**pshrnk)
        goto 1
      else
        hdid=h
        if(errmax.gt.errcon)then
          hnext=safety*h*(errmax**pgrow)
        else
          hnext=4.*h
        endif
      endif
      do 13 i=1,n
        y(i)=y(i)+ytemp(i)*fcor
13    continue
      return
      end


