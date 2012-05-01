ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function rholaw(per)

      implicit none

      real per ! ms

      rholaw = 5.4/sqrt(0.001*per)

      end
