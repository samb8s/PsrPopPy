ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c common set of survey parameters for storing up to "maxsurv" surveys
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer maxsurv
      parameter(maxsurv=16)
      character*80 surveyfile(maxsurv),dmlist(maxsurv)
      common /surveynames/ surveyfile,dmlist
      real beta(maxsurv),gain(maxsurv),tint(maxsurv),
     &     tsamp(maxsurv),tsys(maxsurv),nu(maxsurv),bw(maxsurv),
     &     dnu(maxsurv),fwhm(maxsurv),ramin(maxsurv),ramax(maxsurv),
     &     decmin(maxsurv),decmax(maxsurv),
     &     lomin(maxsurv),lomax(maxsurv),
     &     abmin(maxsurv),abmax(maxsurv),sfrac(maxsurv),
     &     npol(maxsurv),snmin(maxsurv)
      common /surveyparms/ beta,gain,tint,
     &     tsamp,tsys,nu,bw,
     &     dnu,fwhm,ramin,ramax,
     &     decmin,decmax,lomin,lomax,
     &     abmin,abmax,sfrac,npol,snmin	
c
c     optional DM list and tsamp information
c
      integer maxdms,ndms(maxsurv)
      parameter(maxdms=1024)
      common /numberofdms/ ndms
      real dmdd(maxsurv,maxdms),acts(maxsurv,maxdms)
      common /dmandtsamp/ dmdd,acts
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
