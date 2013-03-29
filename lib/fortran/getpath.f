      subroutine getpath(path,lpth)
      implicit none

      character*150 path
      integer lpth

      path="PATHTOCODE"
      lpth=index(path,' ')-1

      end
