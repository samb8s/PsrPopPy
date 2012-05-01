c--------------------------------------------
      subroutine whittle(string,start,finish)
c--------------------------------------------
c
c     routine to whittle the input string
c     of blank characters from start to finish
c
      implicit none
      character*(*) string
      integer start,finish
      start=0
      do while(.true.)
        start=start+1
        if (string(start:start).ne.' ') goto 10
      enddo
 10   finish=start
      do while(.true.)
        finish=finish+1
        if (string(finish:finish).eq.' ') goto 20
      enddo
 20   finish=finish-1
      end
