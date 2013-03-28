ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Quick program to average n_e as a function of Galactocentric radius
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program ner
      implicit none
      real x,y,z,a,r,ne,psr_ne,t,twopi,psrran
      integer i,seed,nav,dmod
      character*80 cline

      call getarg(1,cline)
      if (cline.eq.'') stop 'usage: ner model'
      read(cline,*) dmod

      nav=10000
      seed=0
      twopi=8.0*atan(1.0)
      do a=-2,2,0.1
         ne=0.0
         do i=1,nav
            r=15.0*sqrt(psrran(seed))
            t=twopi*psrran(seed)
            x=r*cos(t)
            y=r*sin(t)
            z=a-0.05+psrran(seed)*0.1
            ne=ne+psr_ne(x,y,z,dmod)
         enddo
         write(*,*) a,ne/real(nav)
      enddo
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
