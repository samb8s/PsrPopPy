ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getrhowidth(per,alpha,iseed,rho,width)

      implicit none

c     -- input: 
      real per, alpha  ! period in ms, alpha in deg
      integer iseed

c     -- output: 
      real rho, width

c     -- local:
      real percut, rho0, drho, drho0, dd, psrran, rholaw
      real logrho, log10, beta, sind, asind

      external sind, asind
      
      data percut/30./
      data drho/0.3/  ! in log10
      data drho0/.3/  ! in log10
      
 1    rho0=rholaw(percut)

      dd=(psrran(iseed)-0.5)
      if (per.gt.percut) then
         rho=5.4/sqrt(0.001*per)
         dd=dd*drho
      else
         rho=rho0
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
