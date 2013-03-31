c==============================================================================
      subroutine dmdist(dm,l,b,mod,distance)
c==============================================================================
c
c     Routine to integrate numerically from the sun to a pulsar given its
c     dispersion measure dm (pc/cc), galactic longitude l and latitude b 
c     (degrees). Uses a model for the free electron density distribution in 
c     the galaxy chosen by the integer mod. The distance is then returned 
c     in kpc once the dm is reached.  See psr_ne for list of distance models
c
      implicit none
c
c     passed down variables
c
      real dm, l, b,distance
      integer mod,lout
c
c     local variables...
c
      real dstep, dmint, maxd, locne, x, y, z, dist, sm, 
     &     smtau, smtheta, smiso, lr, br, dmpsr, dkpc 
      character*1 limit
c
c     functions
c
      real psr_ne
c
c     define parameters & initial values
c
      parameter(dstep=0.01)        ! integrating step in kpc
      parameter(maxd=30.0)         ! max distance to go for kpc
      real s1,s2,s3,s4,dr
      parameter(dr=0.0174532925199)
      logical first
      common /outputtxt/ lout
      data first/.true./
      save
      dmint = 0.0                  ! integrated dm
      if (mod.eq.4) then
          if (first) write(lout,*) 'NE2001 distance model'
          if (first) write(lout,*) 
	  first = .false.
         call dmdsm(l*dr,b*dr,1,dm,distance,limit,s1,s2,s3,s4)
         return
      endif
c
c     Main loop, integrate out to given dm
c
      do dist = dstep, maxd, dstep
c
c     calculate local position
c
         call calc_xyz(l,b,dist,x,y,z)
c
c     local electron density
c
         locne = psr_ne(x, y, z, mod)           
c
c     integrated dm
c
         dmint = dmint + locne * dstep * 1000.0 
c
c     integrated far enough?
c
         if (dmint.eq.dm) then
            distance = dist
            return
         else if (dmint.gt.dm) then
            distance = dist - dstep
            return
         end if
      end do

      end
