/* pxcone.f -- translated by f2c and hacked by Leif Lönnblad to avoid
   linking with libf2c.
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
//#include "Rivet/Projections/pxcone.h"
using namespace std;

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Table of constant values, which are actually non const to be able
   to be used as fortran arguments. */
static int MAXV = 20000;
static int VDIM = 3;
static double twopi = 6.283185307;

// The standard fortran SIGN function for doubles.
double d_sign(double *a, double * b) {
  return *b < 0.0? -fabs(*a): fabs(*a);
}

// The standard fortran MOD function for doubles.
double d_mod(double * a, double * p) {
  return (*a) - int((*a)/(*p))*(*p);
}

// The main PXCONE function.
void pxcone_(int *mode, int *ntrak, int *itkdm, 
	double *ptrak, double *coner, double *epslon, double *
	ovlim, int *mxjet, int *njet, double *pjet, int *
	ipass, int *ijmul, int *ierr)
{
    /* Initialized data */

    static int ncall = 0;
    static int nprint = 0;
    static double rold = 0.;
    static double epsold = 0.;
    static double ovold = 0.;

    /* System generated locals */
    int ptrak_dim1, ptrak_offset, i__1, i__2;
    double d__1, d__2, d__3;

    /* Local variables */
    static double cosr, rsep, ppsq, ptsq, cos2r;
    static int i__, j, n;
    static double vseed[3];
    static int iterr;
    extern /* Subroutine */ int pxord_(double *, int *, int *, 
	    int *, double *);
    static int n1, n2;
    static double pj[20000]	/* was [4][5000] */, pp[20000]	/* was [4][
	    5000] */;
    static int mu;
    static double pu[15000]	/* was [3][5000] */, cosval;
    extern /* Subroutine */ int pxaddv_(int *, double *, double *,
	     double *, int *);
    static int jetlis[25000000]	/* was [5000][5000] */;
    extern double pxmdpi_(double *);
    extern /* Subroutine */ int pxsear_(int *, double *, int *, 
	    double *, double *, double *, int *, int *, 
	    double *, int *, int *), pxolap_(int *, int *,
	     int *, int *, double *, double *, double *);
    static int unstbl;
    extern /* Subroutine */ int pxuvec_(int *, double *, double *,
	     int *), pxzeri_(int *, int *), pxnorv_(int *, 
	    double *, double *, int *), pxzerv_(int *, 
	    double *);
    static double vec1[3], vec2[3];

/* .********************************************************* */
/* . ------ */
/* . PXCONE */
/* . ------ */
/* . */
/* . Code downloaded from the following web page */
/* . */
/* .   http://aliceinfo.cern.ch/alicvs/viewvc/JETAN/pxcone.F?view=markup&pathrev=v4-05-04 */
/* . */
/* . on 17/10/2006 by G. Salam. Permission subsequently granted by Michael */
/* . H. Seymour (on behalf of the PxCone authors) for this code to be */
/* . distributed together with FastJet under the terms of the GNU Public */
/* . License v2 (see the file COPYING in the main FastJet directory). */
/* . */
/* .********** Pre Release Version 26.2.93 */
/* . */
/* . Driver for the Cone  Jet finding algorithm of L.A. del Pozo. */
/* . Based on algorithm from D.E. Soper. */
/* . Finds jets inside cone of half angle CONER with energy > EPSLON. */
/* . Jets which receive more than a fraction OVLIM of their energy from */
/* . overlaps with other jets are excluded. */
/* . Output jets are ordered in energy. */
/* . If MODE.EQ.2 momenta are stored as (eta,phi,<empty>,pt) */
/* . Usage     : */
/* . */
/* .      INTEGER  ITKDM,MXTRK */
/* .      PARAMETER  (ITKDM=4.or.more,MXTRK=1.or.more) */
/* .      INTEGER  MXJET, MXTRAK, MXPROT */
/* .      PARAMETER  (MXJET=10,MXTRAK=500,MXPROT=500) */
/* .      INTEGER  IPASS (MXTRAK),IJMUL (MXJET) */
/* .      INTEGER  NTRAK,NJET,IERR,MODE */
/* .      DOUBLE PRECISION  PTRAK (ITKDM,MXTRK),PJET (5,MXJET) */
/* .      DOUBLE PRECISION  CONER, EPSLON, OVLIM */
/* .      NTRAK = 1.to.MXTRAK */
/* .      CONER   = ... */
/* .      EPSLON  = ... */
/* .      OVLIM   = ... */
/* .      CALL PXCONE (MODE,NTRAK,ITKDM,PTRAK,CONER,EPSLON,OVLIM,MXJET, */
/* .     +             NJET,PJET,IPASS,IJMUL,IERR) */
/* . */
/* . INPUT     :  MODE      1=>e+e-, 2=>hadron-hadron */
/* . INPUT     :  NTRAK     Number of particles */
/* . INPUT     :  ITKDM     First dimension of PTRAK array */
/* . INPUT     :  PTRAK     Array of particle 4-momenta (Px,Py,Pz,E) */
/* . INPUT     :  CONER     Cone size (half angle) in radians */
/* . INPUT     :  EPSLON    Minimum Jet energy (GeV) */
/* . INPUT     :  OVLIM     Maximum fraction of overlap energy in a jet */
/* . INPUT     :  MXJET     Maximum possible number of jets */
/* . OUTPUT    :  NJET      Number of jets found */
/* . OUTPUT    :  PJET      5-vectors of jets */
/* . OUTPUT    :  IPASS(k)  Particle k belongs to jet number IPASS(k) */
/* .                        IPASS = -1 if not assosciated to a jet */
/* . OUTPUT    :  IJMUL(i)  Jet i contains IJMUL(i) particles */
/* . OUTPUT    :  IERR      = 0 if all is OK ;   = -1 otherwise */
/* . */
/* . CALLS     : PXSEAR, PXSAME, PXNEW, PXTRY, PXORD, PXUVEC, PXOLAP */
/* . CALLED    : User */
/* . */
/* . AUTHOR    :  L.A. del Pozo */
/* . CREATED   :  26-Feb-93 */
/* . LAST MOD  :   2-Mar-93 */
/* . */
/* . Modification Log. */
/* . 25-Feb-07: G P Salam   - fix bugs concerning 2pi periodicity in eta phi mode */
/* .                        - added commented code to get consistent behaviour */
/* .                          regardless of particle order (replaces n-way */
/* .                          midpoints with 2-way midpoints however...) */
/* . 2-Jan-97: M Wobisch    - fix bug concerning COS2R in eta phi mode */
/* . 4-Apr-93: M H Seymour  - Change 2d arrays to 1d in PXTRY & PXNEW */
/* . 2-Apr-93: M H Seymour  - Major changes to add boost-invariant mode */
/* . 1-Apr-93: M H Seymour  - Increase all array sizes */
/* . 30-Mar-93: M H Seymour - Change all REAL variables to DOUBLE PRECISION */
/* . 30-Mar-93: M H Seymour - Change OVLIM into an input parameter */
/* . 2-Mar-93: L A del Pozo - Fix Bugs in PXOLAP */
/* . 1-Mar-93: L A del Pozo - Remove Cern library routine calls */
/* . 1-Mar-93: L A del Pozo - Add Print out of welcome and R and Epsilon */
/* . */
/* .********************************************************* */
/* +SEQ,DECLARE. */
/* ** External Arrays */
/* ** Internal Arrays */
/* ** Used in the routine. */
/* MWobisch */
/* MWobisch */
    /* Parameter adjustments */
    --ipass;
    ptrak_dim1 = *itkdm;
    ptrak_offset = 1 + ptrak_dim1 * 1;
    ptrak -= ptrak_offset;
    --ijmul;
    pjet -= 6;

    /* Function Body */
/* MWobisch */
/* *************************************** */
    rsep = 2.;
/* *************************************** */
/* MWobisch */
    *ierr = 0;

/* ** INITIALIZE */
    if (ncall <= 0) {
	rold = (float)0.;
	epsold = (float)0.;
	ovold = (float)0.;
    }
    ++ncall;

/* ** Print welcome and Jetfinder parameters */
    if ((*coner != rold || *epslon != epsold || *ovlim != ovold) && nprint <= 
	    10) {
      printf("%s\n", " *********** PXCONE: Cone Jet-finder ***********");
      printf("%s\n", "    Written by Luis Del Pozo of OPAL");
      printf("%s\n", "    Modified for eta-phi by Mike Seymour");
      printf("%s\n", "    Includes bug fixes by Wobisch, Salam");
      printf("%s\n", "    Translated to c(++) by Leif Lonnblad");
      printf("%s%5.2f%s\n", "    Cone Size R = ",*coner," Radians");
      printf("%s%5.2f%s\n", "    Min Jet energy Epsilon = ",*epslon," GeV");
      printf("%s%5.2f\n", "   Overlap fraction parameter = ",*ovlim);
      printf("%s\n", "    PXCONE is not a supported product and is");
      printf("%s\n", "    is provided for comparative purposes only");
      printf("%s\n", " ***********************************************");
/*         WRITE (6,*) */
/*         WRITE (6,*) ' *********** PXCONE: Cone Jet-finder ***********' */
/*         WRITE (6,*) '    Written by Luis Del Pozo of OPAL' */
/*         WRITE (6,*) '    Modified for eta-phi by Mike Seymour' */
/*         WRITE (6,*) '    Includes bug fixes by Wobisch, Salam' */
/*         WRITE(6,1000)'   Cone Size R = ',CONER,' Radians' */
/*         WRITE(6,1001)'   Min Jet energy Epsilon = ',EPSLON,' GeV' */
/*         WRITE(6,1002)'   Overlap fraction parameter = ',OVLIM */
/*         WRITE (6,*) '    PXCONE is not a supported product and is' */
/*         WRITE (6,*) '    is provided for comparative purposes only' */
/*         WRITE (6,*) ' ***********************************************' */
/* MWobisch */
	if (rsep < (float)1.999) {
          printf("%s\n", " ******************************************");
          printf("%s\n", " ******************************************");
          printf("%s\n", " M Wobisch: private change !!!!!!!!!!!! ");
          printf("%s%5.2f\n", "      Rsep is set to ",rsep);
          printf("%s\n", " this is ONLY meaningful in a NLO calculation");
          printf("%s\n", "      ------------------------  ");
          printf("%s\n", "  please check what you're doing!!");
          printf("%s\n", " ******************************************");
          printf("%s\n", " ******************************************");
          printf("%s\n", " ******************************************");
          printf("%s\n", "");
          printf("%s\n", "");
          printf("%s\n", "");
/*            WRITE(6,*) ' ' */
/*            WRITE (6,*) ' ******************************************' */
/*            WRITE (6,*) ' ******************************************' */
/*            WRITE(6,*) ' M Wobisch: private change !!!!!!!!!!!! ' */
/*            WRITE(6,*) '      Rsep is set to ',RSEP */
/*            WRITE(6,*) ' this is ONLY meaningful in a NLO calculation' */
/*            WRITE(6,*) '      ------------------------  ' */
/*            WRITE(6,*) '  please check what you''re doing!!' */
/*            WRITE(6,*) '   or ask:  Markus.Wobisch@desy.de --' */
/*            WRITE (6,*) ' ******************************************' */
/*            WRITE (6,*) ' ******************************************' */
/*            WRITE (6,*) ' ******************************************' */
/*            WRITE(6,*) ' ' */
/*            WRITE(6,*) ' ' */
	}
/* MWobisch */
/*          WRITE (6,*) */
/* 1000     FORMAT(A18,F5.2,A10) */
/* 1001     FORMAT(A29,F5.2,A5) */
/* 1002     FORMAT(A33,F5.2) */
	++nprint;
	rold = *coner;
	epsold = *epslon;
	ovold = *ovlim;
    }

/* ** Copy calling array PTRAK  to internal array PP(4,NTRAK) */

    if (*ntrak > 5000) {
/*         WRITE (6,*) ' PXCONE: Ntrak too large: ',NTRAK */
      printf("%s%d\n", " PXCONE: Ntrak too large: ", *ntrak);
	*ierr = -1;
	return;
    }
    if (*mode != 2) {
	i__1 = *ntrak;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 4; ++j) {
		pp[j + (i__ << 2) - 5] = ptrak[j + i__ * ptrak_dim1];
/* L101: */
	    }
/* L100: */
	}
    } else {
/* ** Converting to eta,phi,pt if necessary */
	i__1 = *ntrak;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	    d__1 = ptrak[i__ * ptrak_dim1 + 1];
/* Computing 2nd power */
	    d__2 = ptrak[i__ * ptrak_dim1 + 2];
	    ptsq = d__1 * d__1 + d__2 * d__2;
/* Computing 2nd power */
	    d__3 = ptrak[i__ * ptrak_dim1 + 3];
/* Computing 2nd power */
	    d__2 = sqrt(ptsq + d__3 * d__3) + (d__1 = ptrak[i__ * ptrak_dim1 
		    + 3], abs(d__1));
	    ppsq = d__2 * d__2;
	    if (ptsq <= ppsq * (float)4.25e-18) {
		pp[(i__ << 2) - 4] = 20.;
	    } else {
		pp[(i__ << 2) - 4] = log(ppsq / ptsq) * (float).5;
	    }
	    pp[(i__ << 2) - 4] = d_sign(&pp[(i__ << 2) - 4], &ptrak[i__ * 
		    ptrak_dim1 + 3]);
	    if (ptsq == 0.) {
		pp[(i__ << 2) - 3] = 0.;
	    } else {
		pp[(i__ << 2) - 3] = atan2(ptrak[i__ * ptrak_dim1 + 2], ptrak[
			i__ * ptrak_dim1 + 1]);
	    }
	    pp[(i__ << 2) - 2] = 0.;
	    pp[(i__ << 2) - 1] = sqrt(ptsq);
	    pu[i__ * 3 - 3] = pp[(i__ << 2) - 4];
	    pu[i__ * 3 - 2] = pp[(i__ << 2) - 3];
	    pu[i__ * 3 - 1] = pp[(i__ << 2) - 2];
/* L104: */
	}
    }

/* ** Zero output variables */

    *njet = 0;
    i__1 = *ntrak;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 5000; ++j) {
	    jetlis[j + i__ * 5000 - 5001] = false;
/* L103: */
	}
/* L102: */
    }
    pxzerv_(&MAXV, pj);
    pxzeri_(mxjet, &ijmul[1]);

    if (*mode != 2) {
	cosr = cos(*coner);
	cos2r = cos(*coner);
    } else {
/* ** Purely for convenience, work in terms of 1-R**2 */
/* Computing 2nd power */
	d__1 = *coner;
	cosr = 1 - d__1 * d__1;
/* MW -- select Rsep: 1-(Rsep*CONER)**2 */
/* Computing 2nd power */
	d__1 = rsep * *coner;
	cos2r = 1 - d__1 * d__1;
/* ORIGINAL         COS2R =  1-(2*CONER)**2 */
    }
    unstbl = false;
    if (*mode != 2) {
	pxuvec_(ntrak, pp, pu, ierr);
	if (*ierr != 0) {
	    return;
	}
    }
/* ** Look for jets using particle diretions as seed axes */

    i__1 = *ntrak;
    for (n = 1; n <= i__1; ++n) {
	for (mu = 1; mu <= 3; ++mu) {
	    vseed[mu - 1] = pu[mu + n * 3 - 4];
/* L120: */
	}
	pxsear_(mode, &cosr, ntrak, pu, pp, vseed, njet, jetlis, pj, &unstbl, 
		ierr);
	if (*ierr != 0) {
	    return;
	}
/* L110: */
    }
/* MW - for Rsep=1 goto 145 */
/*      GOTO 145 */
/* ** Now look between all pairs of jets as seed axes. */
/*      NJTORG = NJET           ! GPS -- to get consistent behaviour (2-way midpnts) */
/*      DO 140 N1 = 1,NJTORG-1  ! GPS -- to get consistent behaviour (2-way midpnts) */
    i__1 = *njet - 1;
    for (n1 = 1; n1 <= i__1; ++n1) {
	vec1[0] = pj[(n1 << 2) - 4];
	vec1[1] = pj[(n1 << 2) - 3];
	vec1[2] = pj[(n1 << 2) - 2];
	if (*mode != 2) {
	    pxnorv_(&VDIM, vec1, vec1, &iterr);
	}
/*         DO 150 N2 = N1+1,NJTORG ! GPS -- to get consistent behaviour */
	i__2 = *njet;
	for (n2 = n1 + 1; n2 <= i__2; ++n2) {
	    vec2[0] = pj[(n2 << 2) - 4];
	    vec2[1] = pj[(n2 << 2) - 3];
	    vec2[2] = pj[(n2 << 2) - 2];
	    if (*mode != 2) {
		pxnorv_(&VDIM, vec2, vec2, &iterr);
	    }
	    pxaddv_(&VDIM, vec1, vec2, vseed, &iterr);
	    if (*mode != 2) {
		pxnorv_(&VDIM, vseed, vseed, &iterr);
	    } else {
		vseed[0] /= 2;
/* VSEED(2)=VSEED(2)/2 */
/* GPS 25/02/07 */
		d__2 = vec2[1] - vec1[1];
		d__1 = vec1[1] + pxmdpi_(&d__2) * .5;
		vseed[1] = pxmdpi_(&d__1);
	    }
/* ---ONLY BOTHER IF THEY ARE BETWEEN 1 AND 2 CONE RADII APART */
	    if (*mode != 2) {
		cosval = vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * 
			vec2[2];
	    } else {
		if (abs(vec1[0]) >= 20. || abs(vec2[0]) >= 20.) {
		    cosval = -1e3;
		} else {
/* Computing 2nd power */
		    d__1 = vec1[0] - vec2[0];
		    d__3 = vec1[1] - vec2[1];
/* Computing 2nd power */
		    d__2 = pxmdpi_(&d__3);
		    cosval = 1 - (d__1 * d__1 + d__2 * d__2);
		}
	    }
	    if (cosval <= cosr && cosval >= cos2r) {
		pxsear_(mode, &cosr, ntrak, pu, pp, vseed, njet, jetlis, pj, &
			unstbl, ierr);
	    }
/*            CALL PXSEAR(MODE,COSR,NTRAK,PU,PP,VSEED,NJET, */
/*     +           JETLIS,PJ,UNSTBL,IERR) */
	    if (*ierr != 0) {
		return;
	    }
/* L150: */
	}
/* L140: */
    }
    if (unstbl) {
	*ierr = -1;
/*        WRITE (6,*) ' PXCONE: Too many iterations to find a proto-jet' */
        printf(" PXCONE: Too many iterations to find a proto-jet\n");
	return;
    }
/* L145: */
/* ** Now put the jet list into order by jet energy, eliminating jets */
/* ** with energy less than EPSLON. */
    pxord_(epslon, njet, ntrak, jetlis, pj);

/* ** Take care of jet overlaps */
    pxolap_(mode, njet, ntrak, jetlis, pj, pp, ovlim);

/* ** Order jets again as some have been eliminated, or lost energy. */
    pxord_(epslon, njet, ntrak, jetlis, pj);

/* ** All done!, Copy output into output arrays */
    if (*njet > *mxjet) {
/*         WRITE (6,*) ' PXCONE:  Found more than MXJET jets' */
      printf(" PXCONE:  Found more than MXJET jets\n");
	*ierr = -1;
	goto L99;
    }
    if (*mode != 2) {
	i__1 = *njet;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 4; ++j) {
		pjet[j + i__ * 5] = pj[j + (i__ << 2) - 5];
/* L310: */
	    }
/* L300: */
	}
    } else {
	i__1 = *njet;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    pjet[i__ * 5 + 1] = pj[(i__ << 2) - 1] * cos(pj[(i__ << 2) - 3]);
	    pjet[i__ * 5 + 2] = pj[(i__ << 2) - 1] * sin(pj[(i__ << 2) - 3]);
	    pjet[i__ * 5 + 3] = pj[(i__ << 2) - 1] * sinh(pj[(i__ << 2) - 4]);
	    pjet[i__ * 5 + 4] = pj[(i__ << 2) - 1] * cosh(pj[(i__ << 2) - 4]);
/* L315: */
	}
    }
    i__1 = *ntrak;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ipass[i__] = -1;
	i__2 = *njet;
	for (j = 1; j <= i__2; ++j) {
	    if (jetlis[j + i__ * 5000 - 5001]) {
		++ijmul[j];
		ipass[i__] = j;
	    }
/* L330: */
	}
/* L320: */
    }
L99:
    return;
} /* pxcone_ */

/* CMZ :  1.06/00 28/02/94  15.44.44  by  P. Schleper */
/* -- Author : */
/* ----------------------------------------------------------------------- */
/* Subroutine */ int pxnorv_(int *n, double *a, double *b, 
	int *iterr)
{
    /* System generated locals */
    int i__1;
    double d__1;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    static double c__;
    static int i__;

    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
    c__ = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = a[i__];
	c__ += d__1 * d__1;
/* L10: */
    }
    if (c__ <= 0.) {
	return 0;
    }
    c__ = 1 / sqrt(c__);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b[i__] = a[i__] * c__;
/* L20: */
    }
    return 0;
} /* pxnorv_ */

/* CMZ :  2.00/00 10/01/95  10.17.57  by  P. Schleper */
/* CMZ :  1.06/00 15/03/94  12.17.46  by  P. Schleper */
/* -- Author : */

/* +DECK,PXOLAP. */
/* Subroutine */ int pxolap_(int *mode, int *njet, int *ntrak, 
	int *jetlis, double *pj, double *pp, double *ovlim)
{
    /* Initialized data */

    static int ijmin = 0;

    /* System generated locals */
    int i__1, i__2, i__3;
    double d__1, d__2, d__3;

    /* Local variables */
    static int ijet[5000];
    static double thet, cost;
    static int i__, j, n;
    static double eover, thmin;
    static int iterr;
    extern /* Subroutine */ int pxang3_(double *, double *, 
	    double *, double *, int *);
    static int nj, mu;
    static int ovelap;
    extern double pxmdpi_(double *);
    static double vec1[3], vec2[3];


/* ** Looks for particles assigned to more than 1 jet, and reassigns them */
/* ** If more than a fraction OVLIM of a jet's energy is contained in */
/* ** higher energy jets, that jet is neglected. */
/* ** Particles assigned to the jet closest in angle (a la CDF, Snowmass). */
/* +SEQ,DECLARE. */
    /* Parameter adjustments */
    pp -= 5;
    pj -= 5;
    jetlis -= 5001;

    /* Function Body */

    if (*njet <= 1) {
	return 0;
    }
/* ** Look for jets with large overlaps with higher energy jets. */
    i__1 = *njet;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* ** Find overlap energy between jets I and all higher energy jets. */
	eover = (float)0.;
	i__2 = *ntrak;
	for (n = 1; n <= i__2; ++n) {
	    ovelap = false;
	    i__3 = i__ - 1;
	    for (j = 1; j <= i__3; ++j) {
		if (jetlis[i__ + n * 5000] && jetlis[j + n * 5000]) {
		    ovelap = true;
		}
/* L120: */
	    }
	    if (ovelap) {
		eover += pp[(n << 2) + 4];
	    }
/* L110: */
	}
/* ** Is the fraction of energy shared larger than OVLIM? */
	if (eover > *ovlim * pj[(i__ << 2) + 4]) {
/* ** De-assign all particles from Jet I */
	    i__2 = *ntrak;
	    for (n = 1; n <= i__2; ++n) {
		jetlis[i__ + n * 5000] = false;
/* L130: */
	    }
	}
/* L100: */
    }
/* ** Now there are no big overlaps, assign every particle in */
/* ** more than 1 jet to the closet jet. */
/* ** Any particles now in more than 1 jet are assigned to the CLOSET */
/* ** jet (in angle). */
    i__1 = *ntrak;
    for (i__ = 1; i__ <= i__1; ++i__) {
	nj = 0;
	i__2 = *njet;
	for (j = 1; j <= i__2; ++j) {
	    if (jetlis[j + i__ * 5000]) {
		++nj;
		ijet[nj - 1] = j;
	    }
/* L150: */
	}
	if (nj > 1) {
/* ** Particle in > 1 jet - calc angles... */
	    vec1[0] = pp[(i__ << 2) + 1];
	    vec1[1] = pp[(i__ << 2) + 2];
	    vec1[2] = pp[(i__ << 2) + 3];
	    thmin = (float)0.;
	    i__2 = nj;
	    for (j = 1; j <= i__2; ++j) {
		vec2[0] = pj[(ijet[j - 1] << 2) + 1];
		vec2[1] = pj[(ijet[j - 1] << 2) + 2];
		vec2[2] = pj[(ijet[j - 1] << 2) + 3];
		if (*mode != 2) {
		    pxang3_(vec1, vec2, &cost, &thet, &iterr);
		} else {
/* Computing 2nd power */
		    d__1 = vec1[0] - vec2[0];
		    d__3 = vec1[1] - vec2[1];
/* Computing 2nd power */
		    d__2 = pxmdpi_(&d__3);
		    thet = d__1 * d__1 + d__2 * d__2;
		}
		if (j == 1) {
		    thmin = thet;
		    ijmin = ijet[j - 1];
		} else if (thet < thmin) {
		    thmin = thet;
		    ijmin = ijet[j - 1];
		}
/* L160: */
	    }
/* ** Assign track to IJMIN */
	    i__2 = *njet;
	    for (j = 1; j <= i__2; ++j) {
		jetlis[j + i__ * 5000] = false;
/* L170: */
	    }
	    jetlis[ijmin + i__ * 5000] = true;
	}
/* L140: */
    }
/* ** Recompute PJ */
    i__1 = *njet;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (mu = 1; mu <= 4; ++mu) {
	    pj[mu + (i__ << 2)] = (float)0.;
/* L210: */
	}
	i__2 = *ntrak;
	for (n = 1; n <= i__2; ++n) {
	    if (jetlis[i__ + n * 5000]) {
		if (*mode != 2) {
		    for (mu = 1; mu <= 4; ++mu) {
			pj[mu + (i__ << 2)] += pp[mu + (n << 2)];
/* L230: */
		    }
		} else {
		    pj[(i__ << 2) + 1] += pp[(n << 2) + 4] / (pp[(n << 2) + 4]
			     + pj[(i__ << 2) + 4]) * (pp[(n << 2) + 1] - pj[(
			    i__ << 2) + 1]);
/* GPS 25/02/07 */
		    d__2 = pp[(n << 2) + 2] - pj[(i__ << 2) + 2];
		    d__1 = pj[(i__ << 2) + 2] + pp[(n << 2) + 4] / (pp[(n << 
			    2) + 4] + pj[(i__ << 2) + 4]) * pxmdpi_(&d__2);
		    pj[(i__ << 2) + 2] = pxmdpi_(&d__1);
/*                PJ(2,I)=PJ(2,I) */
/*     +               + PP(4,N)/(PP(4,N)+PJ(4,I))*PXMDPI(PP(2,N)-PJ(2,I)) */
		    pj[(i__ << 2) + 4] += pp[(n << 2) + 4];
		}
	    }
/* L220: */
	}
/* L200: */
    }
    return 0;
} /* pxolap_ */

/* CMZ :  2.00/00 10/01/95  10.17.57  by  P. Schleper */
/* CMZ :  1.06/00 14/03/94  15.37.45  by  P. Schleper */
/* -- Author : */

/* +DECK,PXORD. */
/* Subroutine */ int pxord_(double *epslon, int *njet, int *ntrak,
	 int *jetlis, double *pj)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static int i__, j, index[5000];
    static double elist[5000], ptemp[20000]	/* was [4][5000] */;
    static int logtmp[25000000]	/* was [5000][5000] */;
    extern /* Subroutine */ int pxsorv_(int *, double *, int *, 
	    char);


/* ** Routine to put jets into order and eliminate tose less than EPSLON */
/* +SEQ,DECLARE. */
/* ** Puts jets in order of energy: 1 = highest energy etc. */
/* ** Then Eliminate jets with energy below EPSLON */

/* ** Copy input arrays. */
    /* Parameter adjustments */
    pj -= 5;
    jetlis -= 5001;

    /* Function Body */
    i__1 = *njet;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 4; ++j) {
	    ptemp[j + (i__ << 2) - 5] = pj[j + (i__ << 2)];
/* L110: */
	}
	i__2 = *ntrak;
	for (j = 1; j <= i__2; ++j) {
	    logtmp[i__ + j * 5000 - 5001] = jetlis[i__ + j * 5000];
/* L120: */
	}
/* L100: */
    }
    i__1 = *njet;
    for (i__ = 1; i__ <= i__1; ++i__) {
	elist[i__ - 1] = pj[(i__ << 2) + 4];
/* L150: */
    }
/* ** Sort the energies... */
    pxsorv_(njet, elist, index, 'I');
/* ** Fill PJ and JETLIS according to sort ( sort is in ascending order!!) */
    i__1 = *njet;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 4; ++j) {
	    pj[j + (i__ << 2)] = ptemp[j + (index[*njet + 1 - i__ - 1] << 2) 
		    - 5];
/* L210: */
	}
	i__2 = *ntrak;
	for (j = 1; j <= i__2; ++j) {
	    jetlis[i__ + j * 5000] = logtmp[index[*njet + 1 - i__ - 1] + j * 
		    5000 - 5001];
/* L220: */
	}
/* L200: */
    }
/* * Jets are now in order */
/* ** Now eliminate jets with less than Epsilon energy */
    i__1 = *njet;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (pj[(i__ << 2) + 4] < *epslon) {
	    --(*njet);
	    pj[(i__ << 2) + 4] = (float)0.;
	}
/* L300: */
    }
    return 0;
} /* pxord_ */

/* ******************************************************************* */
/* CMZ :  2.00/00 10/01/95  10.17.57  by  P. Schleper */
/* CMZ :  1.06/00 14/03/94  15.37.44  by  P. Schleper */
/* -- Author : */
/* +DECK,PXSEAR. */
/* Subroutine */ int pxsear_(int *mode, double *cosr, int *ntrak, 
	double *pu, double *pp, double *vseed, int *njet, 
	int *jetlis, double *pj, int *unstbl, int *ierr)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int iter;
    static double pnew[4];
    static int n;
    static double naxis[3], oaxis[3];
    extern int pxnew_(int *, int *, int *, int *);
    extern /* Subroutine */ int pxtry_(int *, double *, int *, 
	    double *, double *, double *, double *, 
	    double *, int *, int *);
    static int ok;
    static int mu;
    static int oldlis[5000];
    extern int pxsame_(int *, int *, int *);
    static int newlis[5000];


/* +SEQ,DECLARE. */
/* ** Using VSEED as a trial axis , look for a stable jet. */
/* ** Check stable jets against those already found and add to PJ. */
/* ** Will try up to MXITER iterations to get a stable set of particles */
/* ** in the cone. */

    /* Parameter adjustments */
    pj -= 5;
    jetlis -= 5001;
    --vseed;
    pp -= 5;
    pu -= 4;

    /* Function Body */
    for (mu = 1; mu <= 3; ++mu) {
	oaxis[mu - 1] = vseed[mu];
/* L100: */
    }
    i__1 = *ntrak;
    for (n = 1; n <= i__1; ++n) {
	oldlis[n - 1] = false;
/* L110: */
    }
    for (iter = 1; iter <= 30; ++iter) {
	pxtry_(mode, cosr, ntrak, &pu[4], &pp[5], oaxis, naxis, pnew, newlis, 
		&ok);
/* ** Return immediately if there were no particles in the cone. */
	if (! ok) {
	    return 0;
	}
	if (pxsame_(newlis, oldlis, ntrak)) {
/* ** We have a stable jet. */
	    if (pxnew_(newlis, &jetlis[5001], ntrak, njet)) {
/* ** And the jet is a new one. So add it to our arrays. */
/* ** Check arrays are big anough... */
		if (*njet == 5000) {
/*             WRITE (6,*) ' PXCONE:  Found more than MXPROT proto-jets' */
                  printf(" PXCONE:  Found more than MXPROT proto-jets\n");
		    *ierr = -1;
		    return 0;
		}
		++(*njet);
		i__1 = *ntrak;
		for (n = 1; n <= i__1; ++n) {
		    jetlis[*njet + n * 5000] = newlis[n - 1];
/* L130: */
		}
		for (mu = 1; mu <= 4; ++mu) {
		    pj[mu + (*njet << 2)] = pnew[mu - 1];
/* L140: */
		}
	    }
	    return 0;
	}
/* ** The jet was not stable, so we iterate again */
	i__1 = *ntrak;
	for (n = 1; n <= i__1; ++n) {
	    oldlis[n - 1] = newlis[n - 1];
/* L150: */
	}
	for (mu = 1; mu <= 3; ++mu) {
	    oaxis[mu - 1] = naxis[mu - 1];
/* L160: */
	}
/* L120: */
    }
    *unstbl = true;
    return 0;
} /* pxsear_ */

/* CMZ :  1.06/00 28/02/94  15.44.44  by  P. Schleper */
/* -- Author : */
/* ----------------------------------------------------------------------- */
/* Subroutine */ int pxsorv_(int *n, double *a, int *k, char opt)
{
    /* System generated locals */
    int i__1;

    /* Builtin functions */
    //    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static double b[5000];
    static int i__, j, il[5000], ir[5000];

/*     Sort A(N) into ascending order */
/*     OPT = 'I' : return index array K only */
/*     OTHERWISE : return sorted A and index array K */
/* ----------------------------------------------------------------------- */

/*      INT N,I,J,K(N),IL(NMAX),IR(NMAX) */
/* LUND */

/*      DOUBLE PRECISION A(N),B(NMAX) */
/* LUND */
    /* Parameter adjustments */
    --k;
    --a;

    /* Function Body */
    if (*n > 5000) {
      // WRITE	s_stop("Sorry, not enough room in Mike's PXSORV", (ftnlen)39);
      printf("Sorry, not enough room in Mike's PXSORV\n");
      abort();
    }
    il[0] = 0;
    ir[0] = 0;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	il[i__ - 1] = 0;
	ir[i__ - 1] = 0;
	j = 1;
L2:
	if (a[i__] > a[j]) {
	    goto L5;
	}
/* L3: */
	if (il[j - 1] == 0) {
	    goto L4;
	}
	j = il[j - 1];
	goto L2;
L4:
	ir[i__ - 1] = -j;
	il[j - 1] = i__;
	goto L10;
L5:
	if (ir[j - 1] <= 0) {
	    goto L6;
	}
	j = ir[j - 1];
	goto L2;
L6:
	ir[i__ - 1] = ir[j - 1];
	ir[j - 1] = i__;
L10:
	;
    }
    i__ = 1;
    j = 1;
    goto L8;
L20:
    j = il[j - 1];
L8:
    if (il[j - 1] > 0) {
	goto L20;
    }
L9:
    k[i__] = j;
    b[i__ - 1] = a[j];
    ++i__;
    if ((i__1 = ir[j - 1]) < 0) {
	goto L12;
    } else if (i__1 == 0) {
	goto L30;
    } else {
	goto L13;
    }
L13:
    j = ir[j - 1];
    goto L8;
L12:
    j = -ir[j - 1];
    goto L9;
L30:
    if ( opt == 'I') {
	return 0;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L31: */
	a[i__] = b[i__ - 1];
    }
/* L999: */
    return 0;
} /* pxsorv_ */

/* ******************************************************************** */
/* CMZ :  2.00/00 10/01/95  10.17.57  by  P. Schleper */
/* CMZ :  1.06/00 14/03/94  15.37.44  by  P. Schleper */
/* -- Author : */

/* +DECK,PXTRY. */
/* Subroutine */ int pxtry_(int *mode, double *cosr, int *ntrak, 
	double *pu, double *pp, double *oaxis, double *naxis, 
	double *pnew, int *newlis, int *ok)
{
    /* System generated locals */
    int i__1;
    double d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    static double norm;
    static int n, mu;
    static double cosval;
    extern double pxmdpi_(double *);
    static double normsq;
    static int npp, npu;


/* +SEQ,DECLARE. */
/* ** Note that although PU and PP are assumed to be 2d arrays, they */
/* ** are used as 1d in this routine for efficiency */
/* ** Finds all particles in cone of size COSR about OAXIS direction. */
/* ** Calculates 4-momentum sum of all particles in cone (PNEW) , and */
/* ** returns this as new jet axis NAXIS (Both unit Vectors) */

    /* Parameter adjustments */
    --newlis;
    --pnew;
    --naxis;
    --oaxis;
    --pp;
    --pu;

    /* Function Body */
    *ok = false;
    for (mu = 1; mu <= 4; ++mu) {
	pnew[mu] = (float)0.;
/* L100: */
    }
    npu = -3;
    npp = -4;
    i__1 = *ntrak;
    for (n = 1; n <= i__1; ++n) {
	npu += 3;
	npp += 4;
	if (*mode != 2) {
	    cosval = (float)0.;
	    for (mu = 1; mu <= 3; ++mu) {
		cosval += oaxis[mu] * pu[mu + npu];
/* L120: */
	    }
	} else {
	    if ((d__1 = pu[npu + 1], abs(d__1)) >= 20. || abs(oaxis[1]) >= 
		    20.) {
		cosval = -1e3;
	    } else {
/* Computing 2nd power */
		d__1 = oaxis[1] - pu[npu + 1];
		d__3 = oaxis[2] - pu[npu + 2];
/* Computing 2nd power */
		d__2 = pxmdpi_(&d__3);
		cosval = 1 - (d__1 * d__1 + d__2 * d__2);
	    }
	}
	if (cosval >= *cosr) {
	    newlis[n] = true;
	    *ok = true;
	    if (*mode != 2) {
		for (mu = 1; mu <= 4; ++mu) {
		    pnew[mu] += pp[mu + npp];
/* L130: */
		}
	    } else {
		pnew[1] += pp[npp + 4] / (pp[npp + 4] + pnew[4]) * (pp[npp + 
			1] - pnew[1]);
/*                PNEW(2)=PNEW(2) */
/*     +              + PP(4+NPP)/(PP(4+NPP)+PNEW(4)) */
/*     +               *PXMDPI(PP(2+NPP)-PNEW(2)) */
/* GPS 25/02/07 */
		d__2 = pp[npp + 2] - pnew[2];
		d__1 = pnew[2] + pp[npp + 4] / (pp[npp + 4] + pnew[4]) * 
			pxmdpi_(&d__2);
		pnew[2] = pxmdpi_(&d__1);
		pnew[4] += pp[npp + 4];
	    }
	} else {
	    newlis[n] = false;
	}
/* L110: */
    }
/* ** If there are particles in the cone, calc new jet axis */
    if (*ok) {
	if (*mode != 2) {
	    normsq = (float)0.;
	    for (mu = 1; mu <= 3; ++mu) {
/* Computing 2nd power */
		d__1 = pnew[mu];
		normsq += d__1 * d__1;
/* L140: */
	    }
	    norm = sqrt(normsq);
	} else {
	    norm = 1.;
	}
	for (mu = 1; mu <= 3; ++mu) {
	    naxis[mu] = pnew[mu] / norm;
/* L150: */
	}
    }
    return 0;
} /* pxtry_ */

/* ******************************************************************** */
/* CMZ :  2.00/00 10/01/95  10.17.57  by  P. Schleper */
/* CMZ :  1.06/00 28/02/94  15.44.44  by  P. Schleper */
/* -- Author : */
/* +DECK,PXUVEC. */

/* Subroutine */ int pxuvec_(int *ntrak, double *pp, double *pu, 
	int *ierr)
{
    /* System generated locals */
    int i__1;
    double d__1;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    static int n, mu;
    static double mag;


/* ** Routine to calculate unit vectors PU of all particles PP */
/* +SEQ,DECLARE. */
    /* Parameter adjustments */
    pu -= 4;
    pp -= 5;

    /* Function Body */
    i__1 = *ntrak;
    for (n = 1; n <= i__1; ++n) {
	mag = (float)0.;
	for (mu = 1; mu <= 3; ++mu) {
/* Computing 2nd power */
	    d__1 = pp[mu + (n << 2)];
	    mag += d__1 * d__1;
/* L110: */
	}
	mag = sqrt(mag);
	if (mag == (float)0.) {
/*             WRITE(6,*)' PXCONE: An input particle has zero mod(p)' */
          printf(" PXCONE: An input particle has zero mod(p)\n");
	    *ierr = -1;
	    return 0;
	}
	for (mu = 1; mu <= 3; ++mu) {
	    pu[mu + n * 3] = pp[mu + (n << 2)] / mag;
/* L120: */
	}
/* L100: */
    }
    return 0;
} /* pxuvec_ */

/* CMZ :  1.06/00 28/02/94  15.44.44  by  P. Schleper */
/* -- Author : */
/* ----------------------------------------------------------------------- */
/* Subroutine */ int pxzeri_(int *n, int *a)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__;

    /* Parameter adjustments */
    --a;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__] = 0;
/* L10: */
    }
    return 0;
} /* pxzeri_ */

/* CMZ :  1.06/00 28/02/94  15.44.44  by  P. Schleper */
/* -- Author : */
/* ----------------------------------------------------------------------- */
/*     This is a set of routines written by Mike Seymour to provide the */
/*     services presumably normally provided by standard OPAL routines */
/*     PXZERV zeroes a vector */
/*     PXZERI zeroes a vector of integers */
/*     PXNORV normalizes a vector */
/*     PXADDV adds two vectors */
/*     PXSORV sorts a vector (copied from HERWIG) */
/*     PXANG3 finds the angle (and its cosine) between two vectors */
/*     PXMDPI moves its argument onto the range [-pi,pi) */
/* ----------------------------------------------------------------------- */
/* Subroutine */ int pxzerv_(int *n, double *a)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__;

    /* Parameter adjustments */
    --a;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__] = 0.;
/* L10: */
    }
    return 0;
} /* pxzerv_ */

/* -- Author : */
/* ----------------------------------------------------------------------- */
/* Subroutine */ int pxaddv_(int *n, double *a, double *b, 
	double *c__, int *iterr)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__;

    /* Parameter adjustments */
    --c__;
    --b;
    --a;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c__[i__] = a[i__] + b[i__];
/* L10: */
    }
    return 0;
} /* pxaddv_ */

/* CMZ :  1.06/00 28/02/94  15.44.44  by  P. Schleper */
/* -- Author : */
/* ----------------------------------------------------------------------- */
/* Subroutine */ int pxang3_(double *a, double *b, double *cost, 
	double *thet, int *iterr)
{
    /* System generated locals */
    double d__1, d__2, d__3, d__4, d__5, d__6;

    /* Builtin functions */
    double sqrt(double), acos(double);

    /* Local variables */
    static double c__;

    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
/* Computing 2nd power */
    d__1 = a[1];
/* Computing 2nd power */
    d__2 = a[2];
/* Computing 2nd power */
    d__3 = a[3];
/* Computing 2nd power */
    d__4 = b[1];
/* Computing 2nd power */
    d__5 = b[2];
/* Computing 2nd power */
    d__6 = b[3];
    c__ = (d__1 * d__1 + d__2 * d__2 + d__3 * d__3) * (d__4 * d__4 + d__5 * 
	    d__5 + d__6 * d__6);
    if (c__ <= 0.) {
	return 0;
    }
    c__ = 1 / sqrt(c__);
    *cost = (a[1] * b[1] + a[2] * b[2] + a[3] * b[3]) * c__;
    *thet = acos(*cost);
    return 0;
} /* pxang3_ */

/* CMZ :  1.06/00 14/03/94  15.41.57  by  P. Schleper */
/* -- Author :    P. Schleper   28/02/94 */
int pxnew_(int *tstlis, int *jetlis, int *ntrak, int *
	njet)
{
    /* System generated locals */
    int i__1, i__2;
    int ret_val;

    /* Local variables */
    static int i__, n;
    static int match;
    static int in;


/* ** Note that although JETLIS is assumed to be a 2d array, it */
/* ** it is used as 1d in this routine for efficiency */
/* ** Checks to see if TSTLIS entries correspond to a jet already found */
/* ** and entered in JETLIS */

    /* Parameter adjustments */
    --jetlis;
    --tstlis;

    /* Function Body */
    ret_val = true;
    i__1 = *njet;
    for (i__ = 1; i__ <= i__1; ++i__) {
	match = true;
	in = i__ - 5000;
	i__2 = *ntrak;
	for (n = 1; n <= i__2; ++n) {
	    in += 5000;
	    if (tstlis[n] != jetlis[in]) {
		match = false;
		goto L100;
	    }
/* L110: */
	}
	if (match) {
	    ret_val = false;
	    return ret_val;
	}
L100:
	;
    }
    return ret_val;
} /* pxnew_ */

/* CMZ :  1.06/00 14/03/94  15.41.57  by  P. Schleper */
/* -- Author :    P. Schleper   28/02/94 */
int pxsame_(int *list1, int *list2, int *n)
{
    /* System generated locals */
    int i__1;
    int ret_val;

    /* Local variables */
    static int i__;


/* ** Returns T if the first N elements of LIST1 are the same as the */
/* ** first N elements of LIST2. */

    /* Parameter adjustments */
    --list2;
    --list1;

    /* Function Body */
    ret_val = true;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (list1[i__] != list2[i__]) {
	    ret_val = false;
	    return ret_val;
	}
/* L100: */
    }
    return ret_val;
} /* pxsame_ */

/* CMZ :  1.06/00 28/02/94  15.44.44  by  P. Schleper */
/* -- Author : */
/* ----------------------------------------------------------------------- */
double pxmdpi_(double *phi)
{
    /* System generated locals */
    double ret_val, d__1;


/* ---RETURNS PHI, MOVED ONTO THE RANGE [-PI,PI) */
    ret_val = *phi;
    if (ret_val <= 3.141592654) {
	if (ret_val > -3.141592654) {
	    goto L100;
	} else if (ret_val > -9.424777961) {
	    ret_val += 6.283185307;
	} else {
	    d__1 = 3.141592654 - ret_val;
	    ret_val = -d_mod(&d__1, &twopi) + 3.141592654;
	}
    } else if (ret_val <= 9.424777961) {
	ret_val += -6.283185307;
    } else {
	d__1 = ret_val + 3.141592654;
	ret_val = d_mod(&d__1, &twopi) - 3.141592654;
    }
L100:
    if (abs(ret_val) < 1e-15) {
	ret_val = 0.;
    }
    return ret_val;
} /* pxmdpi_ */

#ifdef __cplusplus
//	}
#endif

int main() {
  int mode = 0;
  int ntrak = 0;
  int itkdm = 0;
  int mxjet = 0;
  int njet = 0;
  int ipass = 0;
  int ijmul = 0;
  int ierr = 0;
  double ptrak = 0.0;
  double coner = 0.0;
  double epslon = 0.0;
  double ovlim = 0.0;
  double pjet = 0.0;

  pxcone_(&mode, &ntrak, &itkdm, &ptrak, &coner, &epslon, &ovlim,
          &mxjet, &njet, &pjet, &ipass, &ijmul, &ierr);
}
#ifdef __cplusplus
	}
#endif
