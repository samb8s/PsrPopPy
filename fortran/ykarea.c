/* ykarea.f -- translated by f2c (version 20050501).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static real c_b3 = (float)1.64;
static real c_b4 = (float)4.01;
static real c_b5 = (float).55;
static real c_b7 = (float)3.51;
static real c_b8 = (float)7.89;
static real c_b9 = (float)0.;
static real c_b11 = (float)500.;

/* ============================================================================== */
real ykarea_(r__, amax, a, b, r1)
real *r__, *amax, *a, *b, *r1;
{
    /* System generated locals */
    real ret_val, r__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double exp(), pow_dd();

    /* Local variables */
    static real x, r0;

/* ============================================================================== */

/*     if amax=0.0, this function returns the integral  x*rho(x) dx */
/*     for the limits 0<x<r, where rho(x) is the function derived by */
/*     Yusifov & Kukuk (2004) */

/*     if amax>0.0, this function returns the value of r when the */
/*     integral is equal to amax. */

    r0 = *r1 + (float)8.5;
    ret_val = (float)0.;
    r__1 = *r__;
    for (x = (float)0.; x <= r__1; x += (float).01) {
	d__1 = (doublereal) ((x + *r1) / r0);
	d__2 = (doublereal) (*a);
	ret_val += pow_dd(&d__1, &d__2) * exp(*b * (float)-1. * (x - (float)
		8.5) / r0) * x * (float).01;
	if (*amax > (float)0. && ret_val > *amax) {
	    ret_val = x;
	    return ret_val;
	}
    }
    return ret_val;
} /* ykarea_ */

/* ============================================================================== */
real ykr_(seed)
integer *seed;
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    extern real ykr0_();

/* ============================================================================== */

/*     uses ykarea to draw a random deviate from Yusifov & Kukuk's */
/*     radial distribution. a=1.64,b=4.01,r1=0.55 */

    ret_val = ykr0_(seed, &c_b3, &c_b4, &c_b5);
    return ret_val;
} /* ykr_ */

/* ============================================================================== */
real llfr_(seed)
integer *seed;
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    extern real ykr0_();

/* ============================================================================== */

/*     uses ykarea to draw a random deviate from Yusifov & Kukuk's */
/*     radial distribution. a=1.64,b=4.01,r1=0.0 */

    ret_val = ykr0_(seed, &c_b7, &c_b8, &c_b9);
    return ret_val;
} /* llfr_ */

/* ============================================================================== */
real ykr0_(seed, a, b, r1)
integer *seed;
real *a, *b, *r1;
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    real ret_val;

    /* Local variables */
    static real area, amax;
    extern real ykarea_(), psrran_();

/* ============================================================================== */

/*     uses ykarea to draw a random deviate from Yusifov & Kukuk's */
/*     radial distribution. a=1.64,b=4.01,r1=0.55 */

    if (first) {
	amax = ykarea_(&c_b11, &c_b9, a, b, r1);
	first = FALSE_;
    }
    area = psrran_(seed) * amax;
    ret_val = ykarea_(&c_b11, &area, a, b, r1);
    return ret_val;
} /* ykr0_ */

