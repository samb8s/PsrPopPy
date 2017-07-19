c==============================================================================
      program popjoin
c==============================================================================
      implicit none
      integer lin,lout,iargc,narg,i,j,maxpops,maxpsrs
      parameter(maxpops=99,maxpsrs=100000)
      character*240 fname(maxpops),cline(maxpops)
      real freq(maxpops),ldeg(maxpsrs),bdeg(maxpsrs),ppsr(maxpsrs)
      real width(maxpsrs),dmp(maxpsrs),tau(maxpsrs),tsky(maxpsrs)
      real dtrue(maxpsrs),lpsr(maxpsrs),s(maxpsrs),si(maxpsrs)
      real dderi(maxpsrs)
      integer lun(maxpops),lcline(maxpops),ngen(maxpops),ntotal
c==============================================================================
      narg=iargc()
      if (narg.lt.3) then
         write(*,*)
         write(*,*) 'POPJOIN: A utility to join individual populations'
         write(*,*)
         write(*,*) 'usage: popjoin pop1 pop2 .... popn popout'
         write(*,*)
         stop
      endif
      if (narg.gt.maxpops) stop 'Too many populations specified!'
c==============================================================================
      ntotal=0
      do i=1,narg
         call getarg(i,fname(i))
         call glun(lun(i))
         open(lun(i),file=fname(i),form="unformatted",status="unknown")
         if (i.lt.narg) then
            read(lun(i)) lcline(i),cline(i)
            read(lun(i)) ngen(i),freq(i)
            if (freq(i).ne.freq(1)) stop 'Frequency inconsistency'
            ntotal=ntotal+ngen(i)
         else
            write(lun(i)) 23,"file created by popjoin"
            write(lun(i)) ntotal,freq(1)
         endif
      enddo
      do i=1,narg-1
	 write(*,'(a)') fname(i)
	 write(*,*) ngen(i)
         do j=1,ngen(i)
            read(lun(i)) ldeg(j),bdeg(j),ppsr(j),width(j),dmp(j),tau(j),
     &                tsky(j),dtrue(j),dderi(j),lpsr(j),s(j),si(j)
	 enddo
	 do j=1,ngen(i)
            write(lun(narg)) ldeg(j),bdeg(j),ppsr(j),width(j),dmp(j),
     &           tau(j),tsky(j),dtrue(j),dderi(j),lpsr(j),s(j),si(j)
         enddo
         close(lun(i))
      enddo
      close(lun(narg))
      end
c==============================================================================
