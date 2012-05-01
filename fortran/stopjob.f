      subroutine stopjob(job)
      character*80 path,job
      logical filex,first
      integer lpth
      data first/.true./
      save
      if (first) then 
         call getpath(path,lpth)
	 first=.false.
      endif
      inquire(file=path(1:lpth)//job,exist=filex)
      if (filex) stop
      end
