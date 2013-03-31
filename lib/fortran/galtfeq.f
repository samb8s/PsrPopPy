c==============================================================================
      subroutine galtfeq(ll, bb, ra, dec, conv)
c==============================================================================
c
c     Converts from/to Galactic coords (ll,bb) to/from equatorial coords (ra,dec)
c     with conv set to +1/-1  respectively. All angles are passed down
c     in degrees using the standard sla utilities to do the hard work.
c     The sla routine converts in the J2000 frame.
c     93/05/26 DRL @ JB Switched to sla_routine
c     Last change 93/07/26 DRL @ JB added back transform.
c
      implicit none
      real ll, bb, ra, dec
      integer conv
c
c     Local variables & pi
c
      double precision ld, bd, rd, dd
      real pi
      parameter(pi=3.1415927)
      if (conv.eq.1) then
c
c       Convert from galactic to equatorial coords
c       Firstly, get l,b into radians and dprecs
c
        ld = dble(pi*ll/180.0)
        bd = dble(pi*bb/180.0)
c
c       Now call the sla routine
c
        call sla_galeq(ld,bd,rd,dd)
c
c       Convert the ra & decs to reals in degrees
c
        ra = rd*180.0/pi
        dec= dd*180.0/pi
      else if (conv.eq.-1) then
c
c       Convert from equatorial to galactic coords
c       Firstly, get ra,dec in radians and dprecs
c
        rd = dble(pi*ra/180.0)
        dd = dble(pi*dec/180.0)
c
c       Now call the sla routine
c
        call sla_eqgal(rd,dd,ld,bd)
c
c       Convert the galactic coords to reals in degrees
c
        ll = ld*180.0/pi
        bb = bd*180.0/pi
      else
c
c       Warn user of incorrect call
c
        write(*,*) 'Unspecified conversion option passed to galtfeq:'
        write(*,*) '+1 gives Galactic   --> Equatorial'
        write(*,*) '-1 gives Equatorial --> Galactic'
        write(*,*) 'No action taken in this case'
      endif
      end

