c==============================================================================
      real function surveysn(seed,pms,wms,dm,tsc,flux,l,b,tsky,fmod,id,
     &                     scindx)
c==============================================================================
c
c     returns S/N for an arbitrary pulsar survey with:
c
c     period (ms) : pms
c     wms (ms)    : wms
c     DM (pc/cc)  : dm
c     scat time   : tsc (ms)
c     flux density: flux (mJy)
c     Gal. long   : l (deg)
c     Gal. lat    : b (deg)
c     Sky back T  : tsky (K)
c     fmod: model pulsar reference frequency (MHz)
c     scindx: spectral index for scattering
c
c     modification history:
c
c     2005/05/16 - drl@jb.man.ac.uk - added RRAT option (signified by P=0)
c                  so that integration time matches effective pulse width.
c
c     2006/01/11 - drl@jb.man.ac.uk - added number of polarizations (was 2)
c
c     2006/02/06 - drl@jb.man.ac.uk - added option to use pksoff for PMSURV
c     2006/03/27 - drl@jb.man.ac.uk - use get_pkspos rather than pksoff
c     2006/06/15 - drl@wvu.edu      - added gain and beamwidth for PMSURV
c     2006/11/06 - drl@wvu.edu      - fixed bug on degfac line - see below
c     2011/10/12 - sam.d.bates@gmail.com - adding scindx parameter
c==============================================================================
      implicit none
      integer seed,i,id
      real pms,wms,dm,flux,l,b,tsky,tsc,eta,ttot,etsc,wint, scindx
      real tband,fmod,ra,dec,ae,psr_normal,wfac,fftsn,pmfft
      real long,psrran,pbeam,offset,weff,tdm,delta,snfac,degfac
      integer lun,lout,beam_number,idx,dmidx
      logical first,inlong,inlat,inra,indec
      include 'survey.inc'
      common /outputtxt/ lout
      data first/.true./
      save
c
c     return baseline sensitivity first time around [sqrt(w/p=0.04)=0.2]
c
c      if (first) then
c         do i=1,maxsurv
c            if (pms.eq.0.0) tint(i)=0.001 ! optimally detect a 1-ms pulse
c            wfac=1.0                      ! no factor for SP searching
c            if (pms.gt.0.0) wfac=0.2      ! sqrt(dc=4%) for periodicity search
c            if (surveyfile(i).ne.'') then
c            write(lout,*) 'Baseline sensitivity for ',
c     &      surveyfile(i)(1:index(surveyfile(i),' ')-1),
c     &      beta(i)*1000.0*tsys(i)*wfac/gain(i)/
c     &      sqrt(npol(i)*bw(i)*tint(i)),' uJy/sigma'
c            endif
c         enddo
c         first=.false.
c      endif
c
c     set S/N to be zero initially!
c
      surveysn=0.0
c
c     is the pulsar in the search area?
c
      long=l
      if (l.gt.180.0) long=l-360.0
c      if (l.lt.0) long=long+360.0
      call galtfeq(l,b,ra,dec,1)
      if (ra.lt.ramin(id).or.ra.gt.ramax(id)) return
      if (dec.lt.decmin(id).or.dec.gt.decmax(id)) return
      if (long.lt.lomin(id).or.long.gt.lomax(id)) return
      if (abs(b).lt.abmin(id).or.abs(b).gt.abmax(id)) return
      if (psrran(seed).gt.sfrac(id)) return
      if (surveyfile(id).eq.'PMSURV') then
c
c        Use the pksoff routine which looks up the nearest beam position
c        and returns the offset (in arcmin) and beam number
c                 
         call get_pkspos(long,b,offset,beam_number)
c
c        Use info from Table 1 of Manchester et al to set gain and fwhm
c        depending on the beam number used in this detection
c
         if (beam_number.eq.1) then
            gain(id)=0.735
            fwhm(id)=14.0
         else if (beam_number.lt.8) then
            gain(id)=0.690
            fwhm(id)=14.1
         else
            gain(id)=0.581
            fwhm(id)=14.5
         endif
      else
c
c        Calculate apparent flux density due to being away from beam centre
c        offset is a randomly chosen offset within FWHM/2
c
         offset = (fwhm(id)*sqrt(psrran(seed))/2.0)
      endif
c
c     -2.7726 = 4 * log(2) ... fixed bug in this line Nov 6, 2006
c     where degfac was getting erroneously multiplied by frequency :(
c     this has no effect on stuff done at 1.4 GHz, but overpredicted
c     numbers for higher frequency surveys and underpredicted lower frq
c
      degfac = exp(-2.7726*offset*offset/fwhm(id)/fwhm(id))
c
c     calculate total noise temperature
c
      ttot = tsys(id) + tsky*(nu(id)/fmod)**(-2.6)
c
c     calculate width-independent S/N factor
c
      snfac= flux*degfac*gain(id)*sqrt(npol(id)*bw(id)*tint(id))
     &        /beta(id)/ttot
      if (ndms(id).eq.0) then
c
c        in absence of a DM table, assume we dedisperse at the exact DM
c
         tband = 0.0 
      else
c
c        use DM table information to find closest trial DM and sample time
c
         idx=dmidx(id,dm)
         tband = 8.3e6*abs(dm-dmdd(id,idx))*bw(id)/nu(id)/nu(id)/nu(id)
         tsamp(id)=acts(id,idx)
      endif
      tdm = 8.3e6*dm*dnu(id)/nu(id)/nu(id)/nu(id) ! dispersion smearing in ms
c
c     calculate effective pulse width
c
      etsc = tsc*(nu(id)/fmod)**(scindx) ! effective scattering time 
      weff = sqrt(wms**2+tsamp(id)**2+tdm**2+etsc**2+tband**2)
      wint=wms
      wms = weff                ! pass back width as the observed value
      if (pms.eq.0.0) then
c
c        this is a RRAT - change the integration time to be equal
c        to the effective pulse width. i.e. assume a perfectly 
c        matched filter
c
         if (wms.gt.1000.0) then
            surveysn=-1.0
            return
         endif
         tint(id)=wms/1000.0 ! note ms->s conversion
         delta=0.5           ! so that multiplication below has no effect
      else
c
c        this is a pulsar, calculate the duty cycle in the standard way
c
         delta = weff/pms
      endif
c
c     Has the pulse been smeared out?
c
      if (delta.gt.1.0) then
         surveysn=-1.0
         return
      endif
c
c     All tests passed. Calculate and return S/N
c
      surveysn = snfac * sqrt((1.0-delta)/delta) 

c         write(66,*) pms,dm,flux,wint,surveysn,delta
c         write(66,*) degfac,gain(id),npol(id),bw(id),tint(id),
c     &        beta(id),ttot
c         stop

c
c     call PMFFTSN
c
c      if (surveyfile(id).eq.'PMSURV'.and.pms.gt.0) then
c         fftsn=2*pmfft(pms,dm,flux*degfac,wms,etsc,beam_number,tsky)
c         if (fftsn.lt.6.0) surveysn=0
c      endif
      end
