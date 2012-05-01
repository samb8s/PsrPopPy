ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Quick program to average n_e as a function of Galactocentric radius
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program ner
      implicit none
      real x,y,z,r0,r,ne,psr_ne,t,twopi,psrran
      integer i,seed,nav,dmod
      character*80 cline

      call getarg(1,cline)
      if (cline.eq.'') stop 'usage: ner model'
      read(cline,*) dmod

      nav=10000
      z=0
      seed=0
      twopi=8.0*atan(1.0)
      do r0=0.25,14.5,0.25
         ne=0.0
         do i=1,nav
            r=r0+psrran(seed)*0.25
            t=twopi*psrran(seed)
            x=r*cos(t)
            y=r*sin(t)
            ne=ne+psr_ne(x,y,z,dmod)
         enddo
         write(*,*) r0+0.125,ne/real(nav)
      enddo
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
