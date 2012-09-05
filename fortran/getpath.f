	subroutine getpath(path,lpth)
	implicit none

	character*80 path
	integer lpth

    path='/Users/sbates/soft/psrpop-GPS'
	lpth=index(path,' ')-1

	end
