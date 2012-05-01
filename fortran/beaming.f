ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function beaming(alpha,rho,seed)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     returns true if pulsar is beaming at us, given inclination angle
c     alpha and beamwidth rho (both in degrees), or false otherwise
c
      implicit none
      integer seed
      real alpha,rho,psrran,beamfrac
      beaming = (psrran(seed).lt.beamfrac(alpha,rho))
      end
