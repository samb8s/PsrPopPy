c==============================================================================
      real function pksoff(ldeg,bdeg)
c==============================================================================
      implicit none
      integer nbeams
      parameter(nbeams=19583)
      real l(nbeams),b(nbeams),ldeg,bdeg,diff,dmin
      integer il,is(200),i,j,k,n,imin,imax,indx,lun,len
      logical first
      character*80 path
      data first/.true./
      save

      if (first) then
         j=1
         call glun(lun)
	 call getpath(path,len)
         len=index(path,' ')-1
         open(lun,file=path(1:len)//"/lookuptables/lb.pksmb")
         is(1)=1
         il=-140
         do i=1,nbeams
            read(lun,*) l(i),b(i)
            k=int(l(i))
            if (k.ne.il) then
               j=j+1
               il=k
               is(j)=i
            endif
         enddo
         close(lun)
         n=j
         first=.false.
      endif

c      goto  1

      imin=max(1,int(ldeg)+137)
      imax=min(n,imin+8)
      dmin=1.0e32
      do i=is(imin),is(imax)
         diff=sqrt( (l(i)-ldeg)**2 + (b(i)-bdeg)**2 )
         if (diff.lt.dmin) then
            dmin=diff
            indx=i
         endif
      enddo

      pksoff=dmin*60.0
      return

 1    continue

      dmin=1.0e32
      do i=1,nbeams
         diff=sqrt( (l(i)-ldeg)**2 + (b(i)-bdeg)**2 )
         if (diff.lt.dmin) then
            dmin=diff
            indx=i
         endif
      enddo
      pksoff=dmin*60.0

      end
