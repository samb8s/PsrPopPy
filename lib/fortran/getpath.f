      subroutine getpath(path,lpth)
      implicit none

      character*150 path
      integer lpth

      path="/Users/sbates/Documents/Physics/Pulsars/code/python/"//
     +"pypop/lib/fortran"
      lpth=index(path,' ')-1

      end
