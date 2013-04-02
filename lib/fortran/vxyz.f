      block data odeint_data
C common for the numerical recipes routine to store the path (note this should
C not be changed without changing ODEINT)
      include 'vxyz.inc'
      data kmax,dxsav/maxpts,0.3/
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine vxyz(eps, x0, y0, z0, nmillion, vx0, vy0, vz0,
     &                    xt, yt, zt, vxt, vyt, vzt, bound)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   eps - accuracy of numerical integration (recommend setting to 0.005) 
c   x0,y0,z0 - initial x,y,z positions (kpc)
c   nmillion - age to integrate (Myr)
c   vx0,vy0,vz0 - initial vx,vy,vz (km/s)
c   xt,yt,zt - final x,y,z positions (kpc)
c   vxt,vyt,vzt - final vx,vy,vz (km/s)
c   bound - logical = .true. if pulsar stays bound to potential
c                   = .false. if pulsar escapes the potential
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      DOUBLE PRECISION DPI,D2PI,DR2DEG,DR2AS,DDEG2R,DAS2R
      REAL SPI,S2PI
*  PI, 2PI
      PARAMETER (DPI=3.141592653589793238462643D0)
      PARAMETER (D2PI=6.283185307179586476925287D0)
      PARAMETER (SPI=3.1415926535897,S2PI=6.28318530717958)
*  RADIANS TO DEGREES
      PARAMETER (DR2DEG=180.D0/DPI)

*  RADIANS TO ARC SECONDS
      PARAMETER (DR2AS=0.2062648062470963560D+06)

*  DEGREES TO RADIANS
      PARAMETER (DDEG2R=1.D0/DR2DEG)

*  ARC SECONDS TO RADIANS
      PARAMETER (DAS2R=0.4848136811095359949D-05)

*  Parsecs to metres
      DOUBLE PRECISION DPC2M
      PARAMETER (DPC2M=3.085678D16)

*  Seconds to years
      DOUBLE PRECISION DSEC2YR
      PARAMETER (DSEC2YR=1.0D0/(60*60*24*365.25))

* kpc/Myr to km/s
      double precision dkpcmy2kms
      parameter (dkpcmy2kms=dpc2m*dsec2yr*1.0d-6)

* Angles to velocities used with angles in mas/yr and distances in kpc
* produces velocities in kms-1
      double precision dang2vel
      PARAMETER(DANG2VEL=4.7405d0)

c      include 'psrpop_parameters.inc'

C starting positions (kpc) and velocities (kms-1)
      real x0, y0, z0, vx0, vy0, vz0

C flight time in millions of years
      real nmillion

C final positions (kpc) and velocities (kms-1)
      real xt, yt, zt, vxt, vyt, vzt
      logical bound

C the initial angular momentum
      REAL LZ0
      COMMON/GALMOD/LZ0

      real y(5)
      real phi, vphi, r, dpdz, dpdr, vrot, cphi, sphi
      real etot

      integer nvar, nok, nbad
      real eps, h1, hmin, rsun
C
      real kggalpot
      external galmodderivs, rkqc
      data nvar,h1,hmin/5,0.1,2.e-5/

      include 'vxyz.inc'

      external odeint_data
      rsun = 8.5
      r = sqrt(x0**2+y0**2)
      phi = atan2(y0,x0)
      call kggalmod(r, 0.0, dpdr, dpdz)
      vrot=sqrt(r*dpdr)
      vphi=(-vx0*sin(phi) + vy0*cos(phi))/dkpcmy2kms - vrot
      vrot = vrot * dkpcmy2kms
c      write(*,*) vphi * dkpcmy2kms

c calculate the total energy to see whether pulsar is bound
      etot = ((vx0+vrot*sin(phi))**2 + (vy0-vrot*cos(phi))**2
     &        +vz0*vz0)*0.5e6+
     &       kggalpot(r,z0)
      bound=(etot.lt.0.0)

C set up initial conditions
      y(1)=r
      y(2)=z0
      y(3)=phi
      y(4)=(vx0*cos(phi)+vy0*sin(phi))/dkpcmy2kms
      y(5)=vz0/dkpcmy2kms
      lz0=r*vphi

C n.b. the integration package has the following units
C distance kpc
C time Myr
C mass in terms of the solar mass
      call odeint(y,nvar,0.0,nmillion, eps, h1, hmin, nok, nbad,
     &              galmodderivs, rkqc)

c$$$      write(*,'(x,2(a,x,i3,x))') 'good steps=',nok,'bad steps=',nbad
      cphi=cos(y(3))
      sphi=sin(y(3))
      
      xt=y(1)*cphi
      yt=y(1)*sphi
      zt=y(2)
      call kggalmod(rsun, 0.0, dpdr, dpdz)
      vrot = sqrt(rsun*dpdr)

C just take off rotation of the earth about gc
C i.e this is in the LSR of the SUN
      vphi = (lz0/y(1))*dkpcmy2kms
      y(4) = y(4)*dkpcmy2kms
      
      vxt = (y(4)*cphi-vphi*sphi)
      vyt = (y(4)*sphi+vphi*cphi) - vrot * dkpcmy2kms
      vzt = y(5)*dkpcmy2kms

      end

      SUBROUTINE ODEINT(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,DERIVS,
     * RKQC)
C adapted from Numerical recipes
      PARAMETER (MAXSTP=100000,NMAX=10,TWO=2.0,ZERO=0.0,TINY=1.E-10)
      include 'vxyz.inc'
      DIMENSION YSTART(NVAR),YSCAL(NMAX),Y(NMAX),DYDX(NMAX)
      X=X1
      H=SIGN(H1,X2-X1)
      NOK=0
      NBAD=0
      KOUNT=0
      DO 11 I=1,NVAR
        Y(I)=YSTART(I)
11    CONTINUE
      XSAV=X-DXSAV*TWO
      DO 16 NSTP=1,MAXSTP
        CALL DERIVS(X,Y,DYDX)
        DO 12 I=1,NVAR
          YSCAL(I)=max(ABS(Y(I))+ABS(H*DYDX(I)),TINY)
12      CONTINUE
        IF(KMAX.GT.0)THEN
          IF(ABS(X-XSAV).GT.ABS(DXSAV)) THEN
            IF(KOUNT.LT.KMAX-1)THEN
              KOUNT=KOUNT+1
              XPLOT(KOUNT)=X
              DO 13 I=1,NVAR
                YPLOT(I,KOUNT)=Y(I)
13            CONTINUE
              XSAV=X
            ENDIF
          ENDIF
        ENDIF
        IF((X+H-X2)*(X+H-X1).GT.ZERO) H=X2-X
        CALL RKQC(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS)
        IF(HDID.EQ.H)THEN
          NOK=NOK+1
        ELSE
          NBAD=NBAD+1
        ENDIF
        IF((X-X2)*(X2-X1).GE.ZERO)THEN
          DO 14 I=1,NVAR
            YSTART(I)=Y(I)
14        CONTINUE
          IF(KMAX.NE.0)THEN
            KOUNT=KOUNT+1
            XPLOT(KOUNT)=X
            DO 15 I=1,NVAR
              YPLOT(I,KOUNT)=Y(I)
15          CONTINUE
          ENDIF
          RETURN
        ENDIF
        IF(ABS(HNEXT).LT.HMIN) then
           write(*,'(x,a,x,e10.4)')
     :       'Stepsize smaller than minimum.',hnext
           hnext=sign(hmin,hnext)
        endif
        H=HNEXT
16    CONTINUE
      PAUSE 'Too many steps.'
      RETURN
      END

C--------------------------------------------------------------------------
      subroutine galmodderivs(X, Y, DYDX)
C real function to supply derivatives for numerical recipes for galactic model from
C Paul Harrison

      REAL X, Y(*), DYDX(*)
      REAL LZ0
      COMMON/GALMOD/LZ0

C Y(1) is R
C Y(2) is z
C Y(3) is phi
C Y(4) is dR/dt
C Y(5) is dz/dt
C and the DYDX array contains their time derivatives
      real dpdz, dpdr

      call kggalmod(y(1),y(2),dpdr,dpdz)
      dydx(1)=y(4)
      dydx(2)=y(5)
      dydx(3)=lz0/(y(1)*y(1))
      dydx(4)=dydx(3)*dydx(3)*y(1)-dpdr
      dydx(5)=-dpdz
      end

C--------------------------------------------------------------------------
      subroutine kggalmod(r,z,dpdr, dpdz)

C Kuijken & Gilmore MNRAS 239 571-603 (89)
C returns the derivatives of the galactic potential
C Paul Harrison

      real r, z, dpdr, dpdz

C  parameter     disk/halo        nucleus        bulge
C   mass(m0)      1.45e11          9.3e9          1e10
C     beta1


C masses in solar masses
      real mdisk, mnuc, mbulge
      data mdisk, mnuc, mbulge/1.45e11,9.3e9,1e10/
      real beta(3)
      data beta/0.4,0.5,0.1/
C distances in Kpc
      real h(3)
      data h/0.325,0.090,0.125/
      real a
      data a/2.4/
      real b_disk, b_nuc, b_bulge
      data b_disk, b_nuc, b_bulge/5.5,0.25,1.5/

C gravitational constant times Msun in units of kpc**3 Myr**-2
      double precision gmsun
      parameter(gmsun=4.498502167e-12)

C TEMP VARIABLES FOR SPEED
      real b1shz,b2shz,b3shz, shz1, shz2, shz3,rbz

      shz1=sqrt(h(1)*h(1)+z*z)
      shz2=sqrt(h(2)*h(2)+z*z)
      shz3=sqrt(h(3)*h(3)+z*z)
      b1shz=beta(1)*shz1
      b2shz=beta(2)*shz2
      b3shz=beta(3)*shz3

C disk/halo
      dpdr=(mdisk*r/(b_disk*b_disk + r*r + 
     -         (a + b1shz + b2shz + b3shz)**2)**1.5)

      dpdz=(mdisk*(beta(1)*z/shz1 + beta(2)*z/shz2 +beta(3)*z/shz3)*
     -     (a + b1shz + b2shz +b3shz)/
     -     (b_disk*b_disk + r*r + (a + b1shz + b2shz + b3shz)**2)**1.5)

C NUCLEUS

      rbz = (z*z+b_nuc*b_nuc+r*r)**(-1.5)
      dpdz=dpdz+mnuc*z*rbz
      dpdr=dpdr+mnuc*r*rbz

C bulge

      rbz = (z*z+b_bulge*b_bulge+r*r)**(-1.5)
      dpdz=(dpdz+mbulge*z*rbz)*gmsun
      dpdr=(dpdr+mbulge*r*rbz)*gmsun
      

      end
C-----------------------------------------------------------------------
      real function kggalpot(r,z)
C Kuijken & Gilmore MNRAS 239 571-603 (89)
C returns the galactic potential (in
C r & z in Kpc
      real r, z

C  parameter     disk/halo        nucleus        bulge
C   mass(m0)      1.45e11          9.3e9          1e10
C     beta1


C masses in solar masses
      real mdisk, mnuc, mbulge
      data mdisk, mnuc, mbulge/1.45e11,9.3e9,1e10/
      real beta(3)
      data beta/0.4,0.5,0.1/
C distances in Kpc
      real h(3)
      data h/0.325,0.090,0.125/
      real a
      data a/2.4/
      real b_disk, b_nuc, b_bulge
      data b_disk, b_nuc, b_bulge/5.5,0.25,1.5/

C gravitational constant times Msun in units of m**3 s**-2
      double precision gmsun
      parameter(gmsun=1.32712497d20)
      double precision dpc2m
      parameter (dpc2m=3.08567802d16)

      
C TEMP VARIABLES FOR SPEED
      real b1shz, b2shz, b3shz, shz1, shz2, shz3, z2, r2

      r2=r*r
      z2=z*z
      shz1=sqrt(h(1)*h(1)+z2)
      shz2=sqrt(h(2)*h(2)+z2)
      shz3=sqrt(h(3)*h(3)+z2)
      b1shz=beta(1)*shz1
      b2shz=beta(2)*shz2
      b3shz=beta(3)*shz3

C disk/halo
      kggalpot=-mdisk/sqrt((a+b1shz+b2shz+b3shz)**2+b_disk*b_disk+r2)

C NUCLEUS
      kggalpot=kggalpot-mnuc/sqrt(z2+r2+b_nuc*b_nuc)
C bulge
      kggalpot=(kggalpot-mbulge/sqrt(z2+r2+b_bulge*b_bulge))*gmsun/
     & (dpc2m*1000.0)
      end
