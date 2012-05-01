c==============================================================================
      subroutine get_pkspos(ldeg,bdeg,offset,ibeam)
c==============================================================================
      implicit none
      integer nbeams,ibeam
      parameter(nbeams=40077)
      real l(nbeams),b(nbeams),ldeg,bdeg,diff,dmin,offset
      integer il,is(200),i,j,k,n,imin,imax,indx,lun,len,beam(nbeams)
      logical first
      character*80 path,junk
      data first/.true./
      save

      if (first) then
         j=1
         call glun(lun)
	 call getpath(path,len)
         len=index(path,' ')-1
         open(lun,file=path(1:len)//"/lookuptables/beam.list")
         is(1)=1
         il=1
         do i=1,nbeams
            beam(i)=0
            read(lun,*) l(i),b(i),junk
            if (junk.eq.'1') beam(i)=1
            if (junk.eq.'2') beam(i)=2
            if (junk.eq.'3') beam(i)=3
            if (junk.eq.'4') beam(i)=4
            if (junk.eq.'5') beam(i)=5
            if (junk.eq.'6') beam(i)=6
            if (junk.eq.'7') beam(i)=7
            if (junk.eq.'8') beam(i)=8
            if (junk.eq.'9') beam(i)=9
            if (junk.eq.'A') beam(i)=10
            if (junk.eq.'B') beam(i)=11
            if (junk.eq.'C') beam(i)=12
            if (junk.eq.'D') beam(i)=13
            if (beam(i).eq.0) write(*,*) l(i),b(i),junk
            if (beam(i).eq.0) stop
            k=int(l(i))+141
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

      imin=max(1,int(ldeg)+137)
      imax=min(n,imin+8)
      dmin=1.0e32
      do i=is(imin),is(imax)
         diff=sqrt( (l(i)-ldeg)**2 + (b(i)-bdeg)**2 )
         if (diff.lt.dmin) then
            ibeam=beam(i)
            dmin=diff
            indx=i
         endif
      enddo

c      dmin=1.0e32
c      do i=1,nbeams
c         diff=sqrt( (l(i)-ldeg)**2 + (b(i)-bdeg)**2 )
c         if (diff.lt.dmin) then
c            ibeam=beam(i)
c            dmin=diff
c            indx=i
c         endif
c      enddo

      offset=dmin*60.0
      end
