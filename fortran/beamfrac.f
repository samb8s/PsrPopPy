ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function beamfrac(alpha,rho)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     returns the fraction of 4pi sr beamed to the sky given a dipolar beam 
c     with inclination alpha and beamwidth rho (both in degrees) using the 
c     formula given in Emmering & Chevalier 89
c
      implicit none
      real alpha,rho,deg2rad,thetal,thetau
      parameter (deg2rad=.01745329278)

      thetal=deg2rad*max(0.0,alpha-rho)
      thetau=deg2rad*min(90.0,alpha+rho)

      beamfrac=cos(thetal)-cos(thetau)
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
