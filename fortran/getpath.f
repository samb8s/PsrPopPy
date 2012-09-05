	subroutine getpath(path,lpth)
	implicit none

	character*80 path
	integer lpth

    path='/Users/sbates/Documents/Physics/Pulsars/code/Python/pypop/fortran'
	lpth=index(path,' ')-1

	end
