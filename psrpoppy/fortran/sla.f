      SUBROUTINE sla_GALEQ (DL, DB, DR, DD)
*+
*     - - - - - -
*      G A L E Q
*     - - - - - -
*
*  Transformation from IAU 1958 galactic coordinates to
*  J2000.0 equatorial coordinates (double precision)
*
*  Given:
*     DL,DB       dp       galactic longitude and latitude L2,B2
*
*  Returned:
*     DR,DD       dp       J2000.0 RA,Dec
*
*  (all arguments are radians)
*
*  Called:
*     sla_DCS2C, sla_DIMXV, sla_DCC2S, sla_DRANRM, sla_DRANGE
*
*  Note:
*     The equatorial coordinates are J2000.0.  Use the routine
*     sla_GE50 if conversion to B1950.0 'FK4' coordinates is
*     required.
*
*  Reference:
*     Blaauw et al, Mon.Not.R.Astron.Soc.,121,123 (1960)
*
*  P.T.Wallace   Starlink   November 1988
*-

      IMPLICIT NONE

      DOUBLE PRECISION DL,DB,DR,DD

      DOUBLE PRECISION sla_DRANRM,sla_DRANGE

      DOUBLE PRECISION V1(3),V2(3)

*
*  L2,B2 system of galactic coordinates
*
*  P = 192.25       RA of galactic north pole (mean B1950.0)
*  Q =  62.6        inclination of galactic to mean B1950.0 equator
*  R =  33          longitude of ascending node
*
*  P,Q,R are degrees
*
*  Equatorial to galactic rotation matrix (J2000.0), obtained by
*  applying the standard FK4 to FK5 transformation, for inertially
*  zero proper motion, to the columns of the B1950 equatorial to
*  galactic rotation matrix:
*
      DOUBLE PRECISION RMAT(3,3)
      DATA RMAT(1,1),RMAT(1,2),RMAT(1,3),
     :     RMAT(2,1),RMAT(2,2),RMAT(2,3),
     :     RMAT(3,1),RMAT(3,2),RMAT(3,3)/
     : -0.054875539726D0,-0.873437108010D0,-0.483834985808D0,
     : +0.494109453312D0,-0.444829589425D0,+0.746982251810D0,
     : -0.867666135858D0,-0.198076386122D0,+0.455983795705D0/



*  Spherical to Cartesian
      CALL sla_DCS2C(DL,DB,V1)

*  Galactic to equatorial
      CALL sla_DIMXV(RMAT,V1,V2)

*  Cartesian to spherical
      CALL sla_DCC2S(V2,DR,DD)

*  Express in conventional ranges
      DR=sla_DRANRM(DR)
      DD=sla_DRANGE(DD)

      END
      SUBROUTINE sla_DMXV (DM, VA, VB)
*+
*     - - - - - 
*      D M X V
*     - - - - -
*
*  Performs the 3-D forward unitary transformation:
*
*     vector VB = matrix DM * vector VA
*
*  (double precision)
*
*  Given:
*     DM       dp(3,3)    matrix
*     VA       dp(3)      vector
*
*  Returned:
*     VB       dp(3)      result vector
*
*  P.T.Wallace   Starlink   March 1986
*-

      IMPLICIT NONE

      DOUBLE PRECISION DM(3,3),VA(3),VB(3)

      INTEGER I,J
      DOUBLE PRECISION W,VW(3)


*  Matrix DM * vector VA -> vector VW
      DO J=1,3
         W=0D0
         DO I=1,3
            W=W+DM(J,I)*VA(I)
         END DO
         VW(J)=W
      END DO

*  Vector VW -> vector VB
      DO J=1,3
         VB(J)=VW(J)
      END DO

      END
      DOUBLE PRECISION FUNCTION sla_DRANGE (ANGLE)
*+
*     - - - - - - -
*      D R A N G E
*     - - - - - - -
*
*  Normalise angle into range +/- pi  (double precision)
*
*  Given:
*     ANGLE     dp      the angle in radians
*
*  The result (double precision) is ANGLE expressed in the range +/- pi.
*
*  P.T.Wallace   Starlink   6 April 1990
*-

      IMPLICIT NONE

      DOUBLE PRECISION ANGLE

      DOUBLE PRECISION DPI,D2PI
      PARAMETER (DPI=3.141592653589793238462643D0)
      PARAMETER (D2PI=6.283185307179586476925287D0)


      sla_DRANGE=MOD(ANGLE,D2PI)
      IF (ABS(sla_DRANGE).GE.DPI)
     :          sla_DRANGE=sla_DRANGE-SIGN(D2PI,ANGLE)

      END
      SUBROUTINE sla_DCC2S (V, A, B)
*+
*     - - - - - -
*      D C C 2 S
*     - - - - - -
*
*  Direction cosines to spherical coordinates (double precision)
*
*  Given:
*     V     d(3)   x,y,z vector
*
*  Returned:
*     A,B   d      spherical coordinates in radians
*
*  The spherical coordinates are longitude (+ve anticlockwise
*  looking from the +ve latitude pole) and latitude.  The
*  Cartesian coordinates are right handed, with the x axis
*  at zero longitude and latitude, and the z axis at the
*  +ve latitude pole.
*
*  If V is null, zero A and B are returned.
*  At either pole, zero A is returned.
*
*  P.T.Wallace   Starlink   July 1989
*-

      IMPLICIT NONE

      DOUBLE PRECISION V(3),A,B

      DOUBLE PRECISION X,Y,Z,R


      X = V(1)
      Y = V(2)
      Z = V(3)
      R = SQRT(X*X+Y*Y)

      IF (R.EQ.0D0) THEN
         A = 0D0
      ELSE
         A = ATAN2(Y,X)
      END IF

      IF (Z.EQ.0D0) THEN
         B = 0D0
      ELSE
         B = ATAN2(Z,R)
      END IF

      END
      SUBROUTINE sla_DIMXV (DM, VA, VB)
*+
*     - - - - - - 
*      D I M X V
*     - - - - - -
*
*  Performs the 3-D backward unitary transformation:
*
*     vector VB = (inverse of matrix DM) * vector VA
*
*  (double precision)
*
*  (n.b.  the matrix must be unitary, as this routine assumes that
*   the inverse and transpose are identical)
*
*  Given:
*     DM       dp(3,3)    matrix
*     VA       dp(3)      vector
*
*  Returned:
*     VB       dp(3)      result vector
*
*  P.T.Wallace   Starlink   March 1986
*-

      IMPLICIT NONE

      DOUBLE PRECISION DM(3,3),VA(3),VB(3)

      INTEGER I,J
      DOUBLE PRECISION W,VW(3)



*  Inverse of matrix DM * vector VA -> vector VW
      DO J=1,3
         W=0D0
         DO I=1,3
            W=W+DM(I,J)*VA(I)
         END DO
         VW(J)=W
      END DO

*  Vector VW -> vector VB
      DO J=1,3
         VB(J)=VW(J)
      END DO

      END
      SUBROUTINE sla_EQGAL (DR, DD, DL, DB)
*+
*     - - - - - -
*      E Q G A L
*     - - - - - -
*
*  Transformation from J2000.0 equatorial coordinates to
*  IAU 1958 galactic coordinates (double precision)
*
*  Given:
*     DR,DD       dp       J2000.0 RA,Dec
*
*  Returned:
*     DL,DB       dp       galactic longitude and latitude L2,B2
*
*  (all arguments are radians)
*
*  Called:
*     sla_DCS2C, sla_DMXV, sla_DCC2S, sla_DRANRM, sla_DRANGE
*
*  Note:
*     The equatorial coordinates are J2000.0.  Use the routine
*     sla_EG50 if conversion from B1950.0 'FK4' coordinates is
*     required.
*
*  Reference:
*     Blaauw et al, Mon.Not.R.Astron.Soc.,121,123 (1960)
*
*  P.T.Wallace   Starlink   November 1988
*-

      IMPLICIT NONE

      DOUBLE PRECISION DR,DD,DL,DB

      DOUBLE PRECISION sla_DRANRM,sla_DRANGE

      DOUBLE PRECISION V1(3),V2(3)

*
*  L2,B2 system of galactic coordinates
*
*  P = 192.25       RA of galactic north pole (mean B1950.0)
*  Q =  62.6        inclination of galactic to mean B1950.0 equator
*  R =  33          longitude of ascending node
*
*  P,Q,R are degrees
*
*  Equatorial to galactic rotation matrix (J2000.0), obtained by
*  applying the standard FK4 to FK5 transformation, for inertially
*  zero proper motion, to the columns of the B1950 equatorial to
*  galactic rotation matrix:
*
      DOUBLE PRECISION RMAT(3,3)
      DATA RMAT(1,1),RMAT(1,2),RMAT(1,3),
     :     RMAT(2,1),RMAT(2,2),RMAT(2,3),
     :     RMAT(3,1),RMAT(3,2),RMAT(3,3)/
     : -0.054875539726D0,-0.873437108010D0,-0.483834985808D0,
     : +0.494109453312D0,-0.444829589425D0,+0.746982251810D0,
     : -0.867666135858D0,-0.198076386122D0,+0.455983795705D0/



*  Spherical to Cartesian
      CALL sla_DCS2C(DR,DD,V1)

*  Equatorial to galactic
      CALL sla_DMXV(RMAT,V1,V2)

*  Cartesian to spherical
      CALL sla_DCC2S(V2,DL,DB)

*  Express in conventional ranges
      DL=sla_DRANRM(DL)
      DB=sla_DRANGE(DB)

      END
      SUBROUTINE sla_DCS2C (A, B, V)
*+
*     - - - - - -
*      D C S 2 C
*     - - - - - -
*
*  Spherical coordinates to direction cosines (double precision)
*
*  Given:
*     A,B       dp      spherical coordinates in radians
*                        (RA,Dec), (Long,Lat) etc
*
*  Returned:
*     V         dp(3)   x,y,z unit vector
*
*  The spherical coordinates are longitude (+ve anticlockwise
*  looking from the +ve latitude pole) and latitude.  The
*  Cartesian coordinates are right handed, with the x axis
*  at zero longitude and latitude, and the z axis at the
*  +ve latitude pole.
*
*  P.T.Wallace   Starlink   October 1984
*-

      IMPLICIT NONE

      DOUBLE PRECISION A,B,V(3)

      DOUBLE PRECISION COSB



      COSB=COS(B)

      V(1)=COS(A)*COSB
      V(2)=SIN(A)*COSB
      V(3)=SIN(B)

      END
      DOUBLE PRECISION FUNCTION sla_DRANRM (ANGLE)
*+
*     - - - - - - -
*      D R A N R M
*     - - - - - - -
*
*  Normalise angle into range 0-2 pi  (double precision)
*
*  Given:
*     ANGLE     dp      the angle in radians
*
*  The result is ANGLE expressed in the range 0-2 pi (double
*  precision).
*
*  P.T.Wallace   Starlink   December 1984
*-

      IMPLICIT NONE

      DOUBLE PRECISION ANGLE

      DOUBLE PRECISION D2PI
      PARAMETER (D2PI=6.283185307179586476925287D0)


      sla_DRANRM=MOD(ANGLE,D2PI)
      IF (sla_DRANRM.LT.0D0) sla_DRANRM=sla_DRANRM+D2PI

      END
