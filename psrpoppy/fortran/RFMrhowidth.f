ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine RFMrhowidth(period, alpha, freq, rho, width, iseed)

c     CALCULATE EMISSION HEIGHT AND CORRESPONDING RHO VALUE
C     EQUATIONS FROM KARASTERGIOU & JOHNSTON 2007
      implicit none
      real pi,rsun,rad
      parameter(pi=3.1415927,rsun=8.5,rad=57.29577951)

c     -- input: 
      real period, freq, rho, width, alpha  ! period in ms, alpha in deg
      integer iseed

c     -- local:
      real percut, drho, drho0, dd, psrran, rholaw, fghz
      real logrho, log10, beta, sind, asind, height_km

      external sind, asind
      
      data percut/30./
      data drho/0.3/  ! in log10
      data drho0/.3/  ! in log10
      
c     USE KARASTERGIOU/MITRA&RANKIN TO COMPUTE EMISSION HEIGHT
c     need freq in GHz

      if (period.gt.percut) then
         fghz = freq/1000.
         height_km = 100. + (fghz/300.)**(-.666)
c        Need to convert height to metres -> tho in radians
         rho = 3. * sqrt(pi * height_km * 1000./ (2. * 3.e8 * period))
         rho = rho * 180./pi ! convert to degrees
         dd=dd*drho
      else
         rho=rholaw(percut)
         dd=dd*drho0
      endif
      logrho = log10(rho) + dd
      rho = 10.**logrho

      beta = (2.0*psrran(iseed)-1.0)*rho
      width = sind(0.5*rho)**2-sind(0.5*beta)**2
      width = width/(sind(alpha)*sind(alpha+beta))
      
      if (width.lt.0.0) then
         width=0.0
         rho=0.0
         return
      endif
         
      width = sqrt(width)

      if (width.gt.1.0) then
         width=0.0
         rho=0.0
         return
      endif
         
      width=asind(width)*4.0

      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
