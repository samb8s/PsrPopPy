c---------------------------------------------
      real function psr_normal(seed, mu, sigma)
c---------------------------------------------
c
c     Gaussian distribution function
c
      real mu, sigma, rnd1, rnd2
      integer*4 seed
      real psrran
      rnd1 = psrran(seed)
      rnd2 = psrran(seed)
      if (rnd1 .eq. 0.0) rnd1 = psrran(seed)
      psr_normal = ((sigma * ((- (2.0 * log(rnd1))) ** 0.5)) 
     &               * cos((2.0 * 3.1415926) * rnd2)) + mu
      end
