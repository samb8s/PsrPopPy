c==============================================================================
      real function psr_exp(seed, origin, scale)
c==============================================================================
c
c     Produces an exponential distribution centred about origin with
c     a 1/e length given by scale. Seed is the seed for the random
c     number generator.
c
c     Created DRL 92/06/21 @ JB
c
      implicit none
      integer seed
      real origin, scale
c
c     Local variables...
c
      real psrran, rn1, rn2, sgn
      
      rn1=psrran(seed)
      rn2=psrran(seed)
      sgn=1.0
      if (rn1.lt.0.5) sgn=-1.0

      psr_exp = origin + sgn * scale * log(rn2)

      end
