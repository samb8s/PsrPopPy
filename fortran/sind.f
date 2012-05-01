ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function sind(x)

      implicit none

      real x, deg2rad ! ms
      parameter (deg2rad=3.1415927/180.0)

      sind = sin(x*deg2rad)

      end
