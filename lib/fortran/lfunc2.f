ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function lfunc2(seed,lmin,lmax,slope)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     returns a luminosity sampled from a power law in log l between 
c     lmin and lmax (n.b. these must passed down in mJy kpc^2)
c     for a given slope such that dlogN/dlogL=slope. 
c
      implicit none
      integer seed
      real lmin,lmax,slope,psrran,p,lmx,lmi,l
      p=psrran(seed)
      if (slope.eq.0.0) then
         l=log10(lmin)+p*(log10(lmax)-log10(lmin))
         l=10.0**l
      else
         lmx=lmax**slope
         lmi=lmin**slope
         l=(p*(lmx-lmi)+lmi)**(1.0/slope)
      endif
      lfunc2=l
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
