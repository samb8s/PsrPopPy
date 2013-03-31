      program underlying
      implicit none
      character*80 filename,cline*240,param
      integer i,lcline,ngen,lpop,nread
      real ldeg,bdeg,ppsr,width,dmp,tau,x,y,z,
     &                tsky,dtrue,dderi,lpsr,s,si,freq
      logical filex

      call getarg(1,filename)
      inquire(file=filename,exist=filex)
      if (.not.filex) then
         write(*,*) 'usage: underlying population_filename num par'
         stop
      endif
      
      call getarg(2,cline)
      read(cline,*) nread
      call getarg(3,param)

      call glun(lpop)
      open(lpop,file=filename,form="unformatted")
      read(lpop) lcline,cline
      read(lpop) ngen,freq
      if (nread.gt.ngen.or.nread.eq.0) nread=ngen
      do i=1,nread
         read(lpop) ldeg,bdeg,ppsr,width,dmp,tau,
     &                tsky,dtrue,dderi,lpsr,s,si
	 if (param.eq.'true_xyz') then
		call calc_xyz(ldeg,bdeg,dtrue,x,y,z)
		write(*,*) x,y,z
	 endif
	 if (param.eq.'derived_xyz') then
		call calc_xyz(ldeg,bdeg,dderi,x,y,z)
		write(*,*) x,y,z
	 endif
         if (param.eq.'period') write(*,*) ppsr
         if (param.eq.'derived_dist') write(*,*) dderi
         if (param.eq.'true_dist') write(*,*) dtrue
         if (param.eq.'true_radial') then
            call calc_xyz(ldeg,bdeg,dtrue,x,y,z)
            write(*,*) sqrt(x*x+y*y)
         else if (param.eq.'derived_radial') then
            call calc_xyz(ldeg,bdeg,dderi,x,y,z)
            write(*,*) sqrt(x*x+y*y)
         else if (param.eq.'luminosity') then
            write(*,*) lpsr
         endif
      enddo
      close(lpop)
      
      end
