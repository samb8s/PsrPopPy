c=============================================================================
      integer function dmidx(nsurv,dm)
c=============================================================================
      include 'survey.inc'
      integer nsurv,i,idx
      real dm,dmin,delta
      dmin=1.0e32
      idx=0
      do i=1,ndms(nsurv)
         delta=abs(dm-dmdd(nsurv,i))
         if (delta.lt.dmin) then
            dmin=delta
            idx=i
         endif
      enddo
      dmidx=idx
      end
