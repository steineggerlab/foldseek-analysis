/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2003 Washington University School of Medicine
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* histogram.c
 * SRE, Sat Jan 20 16:16:17 1996
 *
 * Accumulation, printing, and fitting of score histograms
 * from database searches.
 *
 * CVS $Id: histogram.c,v 1.18 2003/04/14 16:00:16 eddy Exp $
 ************************************************************
 * Basic API:
 *
 * struct histogram_s *h;
 *
 * h = AllocHistogram(min_hint, max_hint, lumpsize);
 *
 * while (getting scores x) AddToHistogram(h, x);
 *
 * ExtremeValueFitHistogram(h, high_hint);
 * PrintASCIIHistogram(fp, h);
 * FreeHistogram(h);
 */


#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <algorithm>

#ifndef MIN
#define MIN(a,b)         (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b)         (((a)>(b))?(a):(b))
#endif

/* Function: AllocHistogram()
 *
 * Purpose:  Allocate and return a histogram structure.
 *           min and max are your best guess. They need
 *           not be absolutely correct; the histogram
 *           will expand dynamically to accomodate scores
 *           that exceed these suggested bounds. The amount
 *           that the histogram grows by is set by "lumpsize".
 *
 * Args:     min:      minimum score (integer)
 *           max:      maximum score (integer)
 *           lumpsize: when reallocating histogram, pad the reallocation
 *                     by this much (saves excessive reallocation)
 */

struct histogram_s {
    int *histogram;		/* counts of hits                     */
    int  min;			/* elem 0 of histogram == min         */
    int  max;                     /* last elem of histogram == max      */
    int  highscore;		/* highest active elem has this score */
    int  lowscore;		/* lowest active elem has this score  */
    int  lumpsize;		/* when resizing, overalloc by this   */
    int  total;			/* total # of hits counted            */

    float *expect;		/* expected counts of hits            */
    int    fit_type;		/* flag indicating distribution type  */
    float  param[3];		/* parameters used for fits           */
    float  chisq;			/* chi-squared val for goodness of fit*/
    float  chip;			/* P value for chisquared             */
};
#define HISTFIT_NONE     0	/* no fit done yet               */
#define HISTFIT_EVD      1	/* fit type = extreme value dist */
#define HISTFIT_GAUSSIAN 2	/* fit type = Gaussian           */
#define EVD_MU		 0	/* EVD fit parameter mu          */
#define EVD_LAMBDA       1	/* EVD fit parameter lambda      */
#define EVD_WONKA        2      /* EVD fit fudge factor          */
#define GAUSS_MEAN       0	/* Gaussian parameter mean       */
#define GAUSS_SD         1	/* Gaussian parameter std. dev.  */





/* Function: EVDDensity()
 * Date:     SRE, Sat Nov 15 19:37:52 1997 [St. Louis]
 *
 * Purpose:  Return the extreme value density P(S=x) at
 *           a given point x, for an EVD controlled by
 *           parameters mu and lambda.
 */
double
EVDDensity(float x, float mu, float lambda)
{
    return (lambda * exp(-1. * lambda * (x - mu)
                         - exp(-1. * lambda * (x - mu))));
}

/* Function: EVDDistribution()
 * Date:     SRE, Tue Nov 18 08:02:22 1997 [St. Louis]
 *
 * Purpose:  Returns the extreme value distribution P(S < x)
 *           evaluated at x, for an EVD controlled by parameters
 *           mu and lambda.
 */
double
EVDDistribution(float x, float mu, float lambda)
{
    return (exp(-1. * exp(-1. * lambda * (x - mu))));
}



/* Function: ExtremeValueP()
 *
 * Purpose:  Calculate P(S>x) according to an extreme
 *           value distribution, given x and the parameters
 *           of the distribution (characteristic
 *           value mu, decay constant lambda).
 *
 *           This function is exquisitely prone to
 *           floating point exceptions if it isn't coded
 *           carefully.
 *
 * Args:     x      = score
 *           mu     = characteristic value of extreme value distribution
 *           lambda = decay constant of extreme value distribution
 *
 * Return:   P(S>x)
 */
double
ExtremeValueP(float x, float mu, float lambda)
{
    double y;
    /* avoid exceptions near P=1.0 */
    /* typical 32-bit sys: if () < -3.6, return 1.0 */
    if ((lambda * (x - mu)) <= -1. * log(-1. * log(DBL_EPSILON))) return 1.0;
    /* avoid underflow fp exceptions near P=0.0*/
    if ((lambda * (x - mu)) >= 2.3 * (double) DBL_MAX_10_EXP)     return 0.0;
    /* a roundoff issue arises; use 1 - e^-x --> x for small x */
    y = exp(-1. * lambda * (x - mu));
    if       (y < 1e-7) return y;
    else     return (1.0 - exp(-1. * y));
}



/* Function: ExtremeValueP2()
 *
 * Purpose:  Calculate P(S>x) in a database of size N,
 *           using P(S>x) for a single sequence, according
 *           to a Poisson distribution.
 *
 * Args:     x      = score
 *           mu     = characteristic value of extreme value distribution
 *           lambda = decay constant of extreme value distribution
 *           N      = number of trials (number of sequences)
 *
 * Return:   P(S>x) for database of size N
 */
double
ExtremeValueP2(float x, float mu, float lambda, int N)
{
    double y;
    y = N * ExtremeValueP(x,mu,lambda);
    if (y < 1e-7) return y;
    else          return (1.0 - exp(-1. * y));
}

/* Function: ExtremeValueE()
 *
 * Purpose:  Calculate E(S>x) in a database of size N,
 *           using P(S>x) for a single sequence: simply np.
 *
 * Args:     x      = score
 *           mu     = characteristic value of extreme value distribution
 *           lambda = decay constant of extreme value distribution
 *           N      = number of trials (number of sequences)
 *
 * Return:   E(S>x) for database of size N
 */
double
ExtremeValueE(float x, float mu, float lambda, int N)
{
    return (double)N * ExtremeValueP(x,mu,lambda);
}


static int sre_randseed = 42;	/* default seed for sre_random()   */

/* Function: sre_random()
 *
 * Purpose:  Return a uniform deviate x, 0.0 <= x < 1.0.
 *
 *           sre_randseed is a static variable, set
 *           by sre_srandom(). When it is non-zero,
 *           we re-seed.
 *
 *           Implements L'Ecuyer's algorithm for combining output
 *           of two linear congruential generators, plus a Bays-Durham
 *           shuffle. This is essentially ran2() from Numerical Recipes,
 *           sans their nonhelpful Rand/McNally-esque code obfuscation.
 *
 *           Overflow errors are avoided by Schrage's algorithm:
 *               az % m = a(z%q) - r(z/q) (+m if <0)
 *           where q=m/a, r=m%a
 *
 *           Requires that long int's have at least 32 bits.
 *           This function uses statics and is NOT THREADSAFE.
 *
 * Reference: Press et al. Numerical Recipes in C, 1992.
 *
 * Reliable and portable, but slow. Benchmarks on wrasse,
 * using Linux gcc and Linux glibc rand() (see randspeed, in Testsuite):
 *     sre_random():    0.5 usec/call
 *     rand():          0.2 usec/call
 */
double
sre_random(void)
{
    static long  rnd1;		/* random number from LCG1 */
    static long  rnd2;            /* random number from LCG2 */
    static long  rnd;             /* random number we return */
    static long  tbl[64];		/* table for Bays/Durham shuffle */
    long x,y;
    int i;

    /* Magic numbers a1,m1, a2,m2 from L'Ecuyer, for 2 LCGs.
     * q,r derive from them (q=m/a, r=m%a) and are needed for Schrage's algorithm.
     */
    long a1 = 40014;
    long m1 = 2147483563;
    long q1 = 53668;
    long r1 = 12211;

    long a2 = 40692;
    long m2 = 2147483399;
    long q2 = 52774;
    long r2 = 3791;

    if (sre_randseed > 0)
    {
        rnd1 = sre_randseed;
        rnd2 = sre_randseed;
        /* Fill the table for Bays/Durham */
        for (i = 0; i < 64; i++) {
            x    = a1*(rnd1%q1);   /* LCG1 in action... */
            y    = r1*(rnd1/q1);
            rnd1 = x-y;
            if (rnd1 < 0) rnd1 += m1;

            x    = a2*(rnd2%q2);   /* LCG2 in action... */
            y    = r2*(rnd2/q2);
            rnd2 = x-y;
            if (rnd2 < 0) rnd2 += m2;

            tbl[i] = rnd1-rnd2;
            if (tbl[i] < 0) tbl[i] += m1;
        }
        sre_randseed = 0;		/* drop the flag. */
    }/* end of initialization*/


    x    = a1*(rnd1%q1);   /* LCG1 in action... */
    y    = r1*(rnd1/q1);
    rnd1 = x-y;
    if (rnd1 < 0) rnd1 += m1;

    x    = a2*(rnd2%q2);   /* LCG2 in action... */
    y    = r2*(rnd2/q2);
    rnd2 = x-y;
    if (rnd2 < 0) rnd2 += m2;

    /* Choose our random number from the table... */
    i   = (int) (((double) rnd / (double) m1) * 64.);
    rnd = tbl[i];
    /* and replace with a new number by L'Ecuyer. */
    tbl[i] = rnd1-rnd2;
    if (tbl[i] < 0) tbl[i] += m1;

    return ((double) rnd / (double) m1);
}



/* Function: EVDrandom()
 *
 * Purpose:  Randomly sample an x from an EVD.
 *           Trivially done by the transformation method, since
 *           the distribution is analytical:
 *              x = \mu - \frac{\log \left[ -\log P(S<x) \right]}{\lambda}
 *           where P(S<x) is sampled uniformly on 0 < P(S<x) < 1.
 */
float
EVDrandom(float mu, float lambda)
{
    float p = 0.0;

    /* Very unlikely, but possible,
     * that sre_random() would give us exactly 0 or 1
     */
    while (p == 0. || p == 1.) p = sre_random();
    return mu - log(-1. * log(p)) / lambda;
}





struct histogram_s *
AllocHistogram(int min, int max, int lumpsize)
{
    struct histogram_s *h;
    int            newsize;
    int            i;

    newsize = max - min + 1;

    h = (struct histogram_s *) malloc(sizeof(struct histogram_s));
    h->min       = min;
    h->max       = max;
    h->total     = 0;
    h->lowscore  = INT_MAX;
    h->highscore = INT_MIN;
    h->lumpsize  = lumpsize;
    h->histogram = (int *) malloc (sizeof(int) * newsize);
    for (i = 0; i < newsize; i++) h->histogram[i] = 0;

    h->expect    = NULL;
    h->fit_type  = HISTFIT_NONE;

    return h;
}


/* Function: FreeHistogram()
 *
 * Purpose:  free a histogram structure.
 */
void
FreeHistogram(struct histogram_s *h)
{
    free(h->histogram);
    if (h->expect != NULL) free(h->expect);
    free(h);
}

/* Function: UnfitHistogram()
 *
 * Purpose:  Free only the theoretical fit part of a histogram.
 */
void
UnfitHistogram(struct histogram_s *h)
{
    if (h->expect != NULL) free(h->expect);
    h->expect   = NULL;
    h->fit_type = HISTFIT_NONE;
}



/* Function: Lawless416()
 * Date:     SRE, Thu Nov 13 11:48:50 1997 [St. Louis]
 *
 * Purpose:  Equation 4.1.6 from [Lawless82], pg. 143, and
 *           its first derivative with respect to lambda,
 *           for finding the ML fit to EVD lambda parameter.
 *           This equation gives a result of zero for the maximum
 *           likelihood lambda.
 *
 *           Can either deal with a histogram or an array.
 *
 *           Warning: beware overflow/underflow issues! not bulletproof.
 *
 * Args:     x      - array of sample values (or x-axis of a histogram)
 *           y      - NULL (or y-axis of a histogram)
 *           n      - number of samples (or number of histogram bins)
 *           lambda - a lambda to test
 *           ret_f  - RETURN: 4.1.6 evaluated at lambda
 *           ret_df - RETURN: first derivative of 4.1.6 evaluated at lambda
 *
 * Return:   (void)
 */
void
Lawless416(float *x, int *y, int n, float lambda, float *ret_f, float *ret_df)
{

    double esum;			/* \sum e^(-lambda xi)      */
    double xesum;			/* \sum xi e^(-lambda xi)   */
    double xxesum;		/* \sum xi^2 e^(-lambda xi) */
    double xsum;			/* \sum xi                  */
    double mult;			/* histogram count multiplier */
    double total;			/* total samples            */
    int i;


    esum = xesum = xsum  = xxesum = total = 0.;
    for (i = 0; i < n; i++)
    {
        mult = (y == NULL) ? 1. : (double) y[i];
        xsum   += mult * x[i];
        xesum  += mult * x[i] * exp(-1. * lambda * x[i]);
        xxesum += mult * x[i] * x[i] * exp(-1. * lambda * x[i]);
        esum   += mult * exp(-1. * lambda * x[i]);
        total  += mult;
    }
    *ret_f  = 1./lambda - xsum / total + xesum / esum;
    *ret_df = ((xesum / esum) * (xesum / esum))
              - (xxesum / esum)
              - (1. / (lambda * lambda));

    return;
}


/* Function: Lawless422()
 * Date:     SRE, Mon Nov 17 09:42:48 1997 [St. Louis]
 *
 * Purpose:  Equation 4.2.2 from [Lawless82], pg. 169, and
 *           its first derivative with respect to lambda,
 *           for finding the ML fit to EVD lambda parameter
 *           for Type I censored data.
 *           This equation gives a result of zero for the maximum
 *           likelihood lambda.
 *
 *           Can either deal with a histogram or an array.
 *
 *           Warning: beware overflow/underflow issues! not bulletproof.
 *
 * Args:     x      - array of sample values (or x-axis of a histogram)
 *           y      - NULL (or y-axis of a histogram)
 *           n      - number of observed samples (or number of histogram bins)
 *           z      - number of censored samples
 *           c      - censoring value; all observed x_i >= c
 *           lambda - a lambda to test
 *           ret_f  - RETURN: 4.2.2 evaluated at lambda
 *           ret_df - RETURN: first derivative of 4.2.2 evaluated at lambda
 *
 * Return:   (void)
 */
void
Lawless422(float *x, int *y, int n, int z, float c,
           float lambda, float *ret_f, float *ret_df)
{
    double esum;			/* \sum e^(-lambda xi)      + z term    */
    double xesum;			/* \sum xi e^(-lambda xi)   + z term    */
    double xxesum;		/* \sum xi^2 e^(-lambda xi) + z term    */
    double xsum;			/* \sum xi                  (no z term) */
    double mult;			/* histogram count multiplier */
    double total;			/* total samples            */
    int i;

    esum = xesum = xsum  = xxesum = total = 0.;
    for (i = 0; i < n; i++)
    {
        mult = (y == NULL) ? 1. : (double) y[i];
        xsum   += mult * x[i];
        esum   += mult *               exp(-1. * lambda * x[i]);
        xesum  += mult * x[i] *        exp(-1. * lambda * x[i]);
        xxesum += mult * x[i] * x[i] * exp(-1. * lambda * x[i]);
        total  += mult;
    }

    /* Add z terms for censored data
     */
    esum   += (double) z *         exp(-1. * lambda * c);
    xesum  += (double) z * c *     exp(-1. * lambda * c);
    xxesum += (double) z * c * c * exp(-1. * lambda * c);

    *ret_f  = 1./lambda - xsum / total + xesum / esum;
    *ret_df = ((xesum / esum) * (xesum / esum))
              - (xxesum / esum)
              - (1. / (lambda * lambda));

    return;
}





/* Function: AddToHistogram()
 *
 * Purpose:  Bump the appropriate counter in a histogram
 *           structure, given a score. The score is
 *           rounded off from float precision to the
 *           next lower integer.
 */
void
AddToHistogram(struct histogram_s *h, float sc)
{
    int score;
    int moveby;
    int prevsize;
    int newsize;
    int i;

    /* Adding to a histogram conflicts with existing fit:
     * prohibit this.
     */
    if (h->fit_type != HISTFIT_NONE)
        exit(-1);


    /* histogram bins are defined as:  score >= bin value, < bin+1
     * -1.9 -> -2    -0.4 -> -1    1.9 -> 1
     * -2.1 -> -3     0.4 -> 0     2.1 -> 2
     */
    score = (int) floor(sc);

    /* Check to see if we must reallocate the histogram.
     */
    if (score < h->min)
    {
        prevsize = h->max - h->min + 1;
        moveby   = (h->min - score) + h->lumpsize;
        newsize  = prevsize + moveby;
        h->min  -= moveby;

        h->histogram = (int *) realloc(h->histogram, sizeof(int) * newsize);
        memmove(h->histogram+moveby, h->histogram, sizeof(int) * prevsize);
        for (i = 0; i < moveby; i++)
            h->histogram[i] = 0;
    }
    else if (score > h->max)
    {
        prevsize = h->max - h->min + 1;
        h->max   = h->lumpsize + score;
        newsize  = h->max - h->min + 1;

        h->histogram = (int *) realloc(h->histogram, sizeof(int) * newsize);
        for (i = prevsize; i < newsize; i++)
            h->histogram[i] = 0;
    }

    /* Bump the correct bin.
     * The bin number is score - h->min
     */
    h->histogram[score - h->min]++;
    h->total++;
    if (score < h->lowscore) h->lowscore   = score;
    if (score > h->highscore) h->highscore = score;

    printf("AddToHistogram(): added %.1f; rounded to %d; in bin %d (%d-%d)\n",
            sc, score, score-h->min, h->min, h->max);
    return;
}


/* Function: EVDMaxLikelyFit()
 * Date:     SRE, Fri Nov 14 07:56:29 1997 [St. Louis]
 *
 * Purpose:  Given a list or a histogram of EVD-distributed samples,
 *           find maximum likelihood parameters lambda and
 *           mu.
 *
 * Algorithm: Uses approach described in [Lawless82]. Solves
 *           for lambda using Newton/Raphson iterations;
 *           then substitutes lambda into Lawless' equation 4.1.5
 *           to get mu.
 *
 *           Newton/Raphson algorithm developed from description in
 *           Numerical Recipes in C [Press88].
 *
 * Args:     x          - list of EVD distributed samples or x-axis of histogram
 *           c          - NULL, or y-axis of histogram
 *           n          - number of samples, or number of histogram bins
 *           ret_mu     : RETURN: ML estimate of mu
 *           ret_lambda : RETURN: ML estimate of lambda
 *
 * Return:   1 on success; 0 on any failure
 */
int
EVDMaxLikelyFit(float *x, int *c, int n, float *ret_mu, float *ret_lambda)
{
    float  lambda, mu;
    float  fx;			/* f(x)  */
    float  dfx;			/* f'(x) */
    double esum;                  /* \sum e^(-lambda xi) */
    double mult;
    double total;
    float  tol = 1e-5;
    int    i;

    /* 1. Find an initial guess at lambda: linear regression here?
     */
    lambda = 0.2;

    /* 2. Use Newton/Raphson to solve Lawless 4.1.6 and find ML lambda
     */
    for (i = 0; i < 100; i++)
    {
        Lawless416(x, c, n, lambda, &fx, &dfx);
        if (fabs(fx) < tol) break;             /* success */
        lambda = lambda - fx / dfx;	     /* Newton/Raphson is simple */
        if (lambda <= 0.) lambda = 0.001;      /* but be a little careful  */
    }

    /* 2.5: If we did 100 iterations but didn't converge, Newton/Raphson failed.
     *      Resort to a bisection search. Worse convergence speed
     *      but guaranteed to converge (unlike Newton/Raphson).
     *      We assume (!?) that fx is a monotonically decreasing function of x;
     *      i.e. fx > 0 if we are left of the root, fx < 0 if we
     *      are right of the root.
     */
    if (i == 100)
    {
        float left, right, mid;
        printf(("EVDMaxLikelyFit(): Newton/Raphson failed; switchover to bisection"));

        /* First we need to bracket the root */
        lambda = right = left = 0.2;
        Lawless416(x, c, n, lambda, &fx, &dfx);
        if (fx < 0.)
        {			/* fix right; search left. */
            do {
                left -= 0.1;
                if (left < 0.) {
                    printf(("EVDMaxLikelyFit(): failed to bracket root"));
                    return 0;
                }
                Lawless416(x, c, n, left, &fx, &dfx);
            } while (fx < 0.);
        }
        else
        {			/* fix left; search right. */
            do {
                right += 0.1;
                Lawless416(x, c, n, right, &fx, &dfx);
                if (right > 100.) {
                    printf(("EVDMaxLikelyFit(): failed to bracket root"));
                    return 0;
                }
            } while (fx > 0.);
        }
        /* now we bisection search in left/right interval */
        for (i = 0; i < 100; i++)
        {
            mid = (left + right) / 2.;
            Lawless416(x, c, n, mid, &fx, &dfx);
            if (fabs(fx) < tol) break;             /* success */
            if (fx > 0.)	left = mid;
            else          right = mid;
        }
        if (i == 100) {
            printf(("EVDMaxLikelyFit(): even the bisection search failed"));
            return 0;
        }
        lambda = mid;
    }

    /* 3. Substitute into Lawless 4.1.5 to find mu
     */
    esum = 0.;
    total = 0.;
    for (i = 0; i < n; i++)
    {
        mult   = (c == NULL) ? 1. : (double) c[i];
        esum  += mult * exp(-1 * lambda * x[i]);
        total += mult;
    }
    mu = -1. * log(esum / total) / lambda;

    *ret_lambda = lambda;
    *ret_mu     = mu;
    return 1;
}


/* Function: EVDCensoredFit()
 * Date:     SRE, Mon Nov 17 10:01:05 1997 [St. Louis]
 * 
 * Purpose:  Given a /left-censored/ list or histogram of EVD-distributed 
 *           samples, as well as the number of censored samples z and the
 *           censoring value c,              
 *           find maximum likelihood parameters lambda and
 *           mu. 
 *           
 * Algorithm: Uses approach described in [Lawless82]. Solves
 *           for lambda using Newton/Raphson iterations;
 *           then substitutes lambda into Lawless' equation 4.2.3
 *           to get mu. 
 *           
 *           Newton/Raphson algorithm developed from description in
 *           Numerical Recipes in C [Press88]. 
 *           
 * Args:     x          - list of EVD distributed samples or x-axis of histogram
 *           y          - NULL, or y-axis of histogram
 *           n          - number of observed samples,or number of histogram bins 
 *           z          - number of censored samples
 *           c          - censoring value (all x_i >= c)
 *           ret_mu     : RETURN: ML estimate of mu
 *           ret_lambda : RETURN: ML estimate of lambda
 *           
 * Return:   (void)
 */
int
EVDCensoredFit(float *x, int *y, int n, int z, float c,
               float *ret_mu, float *ret_lambda)
{
    float  lambda, mu;
    float  fx;			/* f(x)  */
    float  dfx;			/* f'(x) */
    double esum;                  /* \sum e^(-lambda xi) */
    double mult;
    double total;
    float  tol = 1e-5;
    int    i;

    /* 1. Find an initial guess at lambda: linear regression here?
     */
    lambda = 0.2;

    /* 2. Use Newton/Raphson to solve Lawless 4.2.2 and find ML lambda
     */
    for (i = 0; i < 100; i++)
    {
        Lawless422(x, y, n, z, c, lambda, &fx, &dfx);
        if (fabs(fx) < tol) break;             /* success */
        lambda = lambda - fx / dfx;	     /* Newton/Raphson is simple */
        if (lambda <= 0.) lambda = 0.001;      /* but be a little careful  */
    }

    /* 2.5: If we did 100 iterations but didn't converge, Newton/Raphson failed.
      *      Resort to a bisection search. Worse convergence speed
      *      but guaranteed to converge (unlike Newton/Raphson).
      *      We assume (!?) that fx is a monotonically decreasing function of x;
      *      i.e. fx > 0 if we are left of the root, fx < 0 if we
      *      are right of the root.
      */
    if (i == 100)
    {
        float left, right, mid;
        /* First we need to bracket the root */
        printf(("EVDCensoredFit(): Newton/Raphson failed; switched to bisection"));
        lambda = right = left = 0.2;
        Lawless422(x, y, n, z, c, lambda, &fx, &dfx);
        if (fx < 0.)
        {			/* fix right; search left. */
            do {
                left -= 0.03;
                if (left < 0.) {
                    printf(("EVDCensoredFit(): failed to bracket root"));
                    return 0;
                }
                Lawless422(x, y, n, z, c, left, &fx, &dfx);
            } while (fx < 0.);
        }
        else
        {			/* fix left; search right. */
            do {
                right += 0.1;
                Lawless422(x, y, n, z, c, left, &fx, &dfx);
                if (right > 100.) {
                    printf(("EVDCensoredFit(): failed to bracket root"));
                    return 0;
                }
            } while (fx > 0.);
        }
        /* now we bisection search in left/right interval */
        for (i = 0; i < 100; i++)
        {
            mid = (left + right) / 2.;
            Lawless422(x, y, n, z, c, left, &fx, &dfx);
            if (fabs(fx) < tol) break;             /* success */
            if (fx > 0.)	left = mid;
            else          right = mid;
        }
        if (i == 100) {
            printf(("EVDCensoredFit(): even the bisection search failed"));
            return 0;
        }
        lambda = mid;
    }

    /* 3. Substitute into Lawless 4.2.3 to find mu
     */
    esum =  total = 0.;
    for (i = 0; i < n; i++)
    {
        mult   = (y == NULL) ? 1. : (double) y[i];
        esum  += mult * exp(-1. * lambda * x[i]);
        total += mult;
    }
    esum += (double) z * exp(-1. * lambda * c);    /* term from censored data */
    mu = -1. * log(esum / total) / lambda;

    *ret_lambda = lambda;
    *ret_mu     = mu;
    return 1;
}

/* Function: Gammln()
 *
 * Returns the natural log of the gamma function of x.
 * x is > 0.0.
 *
 * Adapted from a public domain implementation in the
 * NCBI core math library. Thanks to John Spouge and
 * the NCBI. (According to the NCBI, that's Dr. John
 * "Gammas Galore" Spouge to you, pal.)
 */
double
Gammln(double x)
{
    int i;
    double xx, tx;
    double tmp, value;
    static double cof[11] = {
            4.694580336184385e+04,
            -1.560605207784446e+05,
            2.065049568014106e+05,
            -1.388934775095388e+05,
            5.031796415085709e+04,
            -9.601592329182778e+03,
            8.785855930895250e+02,
            -3.155153906098611e+01,
            2.908143421162229e-01,
            -2.319827630494973e-04,
            1.251639670050933e-10
    };

    /* Protect against x=0. We see this in Dirichlet code,
     * for terms alpha = 0. This is a severe hack but it is effective
     * and (we think?) safe. (due to GJM)
     */
    if (x <= 0.0) return 999999.;

    xx       = x - 1.0;
    tx = tmp = xx + 11.0;
    value    = 1.0;
    for (i = 10; i >= 0; i--)	/* sum least significant terms first */
    {
        value += cof[i] / tmp;
        tmp   -= 1.0;
    }
    value  = log(value);
    tx    += 0.5;
    value += 0.918938533 + (xx+0.5)*log(tx) - tx;
    return value;
}

/* Function: IncompleteGamma()
 *
 * Purpose:  Returns 1 - P(a,x) where:
 *           P(a,x) = \frac{1}{\Gamma(a)} \int_{0}^{x} t^{a-1} e^{-t} dt
 *                  = \frac{\gamma(a,x)}{\Gamma(a)}
 *                  = 1 - \frac{\Gamma(a,x)}{\Gamma(a)}
 *
 *           Used in a chi-squared test: for a X^2 statistic x
 *           with v degrees of freedom, call:
 *                  p = IncompleteGamma(v/2., x/2.)
 *           to get the probability p that a chi-squared value
 *           greater than x could be obtained by chance even for
 *           a correct model. (i.e. p should be large, say
 *           0.95 or more).
 *
 * Method:   Based on ideas from Numerical Recipes in C, Press et al.,
 *           Cambridge University Press, 1988.
 *
 * Args:     a  - for instance, degrees of freedom / 2     [a > 0]
 *           x  - for instance, chi-squared statistic / 2  [x >= 0]
 *
 * Return:   1 - P(a,x).
 */
double
IncompleteGamma(double a, double x)
{
    int iter;			/* iteration counter */

    if (a <= 0.) exit(-1);
    if (x <  0.) exit(-1);

    /* For x > a + 1 the following gives rapid convergence;
     * calculate 1 - P(a,x) = \frac{\Gamma(a,x)}{\Gamma(a)}:
     *     use a continued fraction development for \Gamma(a,x).
     */
    if (x > a+1)
    {
        double oldp;		/* previous value of p    */
        double nu0, nu1;		/* numerators for continued fraction calc   */
        double de0, de1;		/* denominators for continued fraction calc */

        nu0 = 0.;			/* A_0 = 0       */
        de0 = 1.;			/* B_0 = 1       */
        nu1 = 1.;			/* A_1 = 1       */
        de1 = x;			/* B_1 = x       */

        oldp = nu1;
        for (iter = 1; iter < 100; iter++)
        {
            /* Continued fraction development:
             * set A_j = b_j A_j-1 + a_j A_j-2
             *     B_j = b_j B_j-1 + a_j B_j-2
                 * We start with A_2, B_2.
             */
            /* j = even: a_j = iter-a, b_j = 1 */
            /* A,B_j-2 are in nu0, de0; A,B_j-1 are in nu1,de1 */
            nu0 = nu1 + ((double)iter - a) * nu0;
            de0 = de1 + ((double)iter - a) * de0;

            /* j = odd: a_j = iter, b_j = x */
            /* A,B_j-2 are in nu1, de1; A,B_j-1 in nu0,de0 */
            nu1 = x * nu0 + (double) iter * nu1;
            de1 = x * de0 + (double) iter * de1;

            /* rescale */
            if (de1 != 0.)
            {
                nu0 /= de1;
                de0 /= de1;
                nu1 /= de1;
                de1 =  1.;
            }
            /* check for convergence */
            if (fabs((nu1-oldp)/nu1) < 1.e-7)
                return nu1 * exp(a * log(x) - x - Gammln(a));

            oldp = nu1;
        }
        exit(-1);
    }
    else /* x <= a+1 */
    {
        double p;			/* current sum               */
        double val;		/* current value used in sum */

        /* For x <= a+1 we use a convergent series instead:
         *   P(a,x) = \frac{\gamma(a,x)}{\Gamma(a)},
         * where
         *   \gamma(a,x) = e^{-x}x^a \sum_{n=0}{\infty} \frac{\Gamma{a}}{\Gamma{a+1+n}} x^n
         * which looks appalling but the sum is in fact rearrangeable to
         * a simple series without the \Gamma functions:
         *   = \frac{1}{a} + \frac{x}{a(a+1)} + \frac{x^2}{a(a+1)(a+2)} ...
         * and it's obvious that this should converge nicely for x <= a+1.
         */

        p = val = 1. / a;
        for (iter = 1; iter < 10000; iter++)
        {
            val *= x / (a+(double)iter);
            p   += val;

            if (fabs(val/p) < 1.e-7)
                return 1. - p * exp(a * log(x) - x - Gammln(a));
        }
        exit(-1);
    }
    /*NOTREACHED*/
    return 0.;
}



/* Function: ExtremeValueSetHistogram()
 *
 * Purpose:  Instead of fitting the histogram to an EVD,
 *           simply set the EVD parameters from an external source.
 *
 * Args:     h        - the histogram to set
 *           mu       - mu location parameter
 *           lambda   - lambda scale parameter
 *           lowbound - low bound of the histogram that was fit
 *           highbound- high bound of histogram that was fit
 *           ndegrees - extra degrees of freedom to subtract in X^2 test:
 *                        typically 0 if mu, lambda are parametric,
 *                        else 2 if mu, lambda are estimated from data
 */
void
ExtremeValueSetHistogram(struct histogram_s *h, float mu, float lambda,
                         float lowbound, float highbound, int ndegrees)
{
    int   sc;
    int   hsize, idx;
    int   nbins;
    float delta;

    UnfitHistogram(h);
    h->fit_type          = HISTFIT_EVD;
    h->param[EVD_LAMBDA] = lambda;
    h->param[EVD_MU]     = mu;

    hsize     = h->max - h->min + 1;
    h->expect = (float *) malloc(sizeof(float) * hsize);
    for (idx = 0; idx < hsize; idx++)
        h->expect[idx] = 0.;

    /* Calculate the expected values for the histogram.
     */
    for (sc = h->min; sc <= h->max; sc++)
        h->expect[sc - h->min] =
                ExtremeValueE((float)(sc), h->param[EVD_MU], h->param[EVD_LAMBDA],
                              h->total) -
                ExtremeValueE((float)(sc+1), h->param[EVD_MU], h->param[EVD_LAMBDA],
                              h->total);

    /* Calculate the goodness-of-fit (within whole region)
     */
    h->chisq = 0.;
    nbins    = 0;
    for (sc = lowbound; sc <= highbound; sc++)
        if (h->expect[sc-h->min] >= 5. && h->histogram[sc-h->min] >= 5)
        {
            delta = (float) h->histogram[sc-h->min] - h->expect[sc-h->min];
            h->chisq += delta * delta / h->expect[sc-h->min];
            nbins++;
        }

    /* Since we fit the whole histogram, there is at least
     * one constraint on chi-square: the normalization to h->total.
     */
    if (nbins > 1 + ndegrees)
        h->chip = (float) IncompleteGamma((double)(nbins-1-ndegrees)/2.,
                                          (double) h->chisq/2.);
    else
        h->chip = 0.;
}

/* Function: ExtremeValueFitHistogram()
 * Date:     SRE, Sat Nov 15 17:16:15 1997 [St. Louis]
 *
 * Purpose:  Fit a score histogram to the extreme value
 *           distribution. Set the parameters lambda
 *           and mu in the histogram structure. Calculate
 *           a chi-square test as a measure of goodness of fit.
 *
 * Methods:  Uses a maximum likelihood method [Lawless82].
 *           Lower outliers are removed by censoring the data below the peak.
 *           Upper outliers are removed iteratively using method
 *           described by [Mott92].
 *
 * Args:     h         - histogram to fit
 *           censor    - TRUE to censor data left of the peak
 *           high_hint - score cutoff; above this are `real' hits that aren't fit
 *
 * Return:   1 if fit is judged to be valid.
 *           else 0 if fit is invalid (too few seqs.)
 */
int
ExtremeValueFitHistogram(struct histogram_s *h, int censor, float high_hint)
{
    float *x;                     /* array of EVD samples to fit */
    int   *y;                     /* histogram counts            */
    int    n;			/* number of observed samples  */
    int    z;			/* number of censored samples  */
    int    hsize;			/* size of histogram           */
    float  lambda, mu;		/* new estimates of lambda, mu */
    int    sc;		        /* loop index for score        */
    int    lowbound;		/* lower bound of fitted region*/
    int    highbound;		/* upper bound of fitted region*/
    int    new_highbound;
    int    iteration;

    /* Determine lower bound on fitted region;
     * if we're censoring the data, choose the peak of the histogram.
     * if we're not, then we take the whole histogram.
     */
    lowbound = h->lowscore;
    if (censor)
    {
        int max = -1;
        for (sc = h->lowscore; sc <= h->highscore; sc++)
            if (h->histogram[sc - h->min] > max)
            {
                max      = h->histogram[sc - h->min];
                lowbound = sc;
            }
    }

    /* Determine initial upper bound on fitted region.
     */
    highbound = MIN(high_hint, h->highscore);

    /* Now, iteratively converge on our lambda, mu:
     */
    for (iteration = 0; iteration < 100; iteration++)
    {
        /* Construct x, y vectors.
         */
        x = NULL;
        y = NULL;
        hsize = highbound - lowbound + 1;
        if (hsize < 5) goto FITFAILED; /* require at least 5 bins or we don't fit */

        x = (float*) malloc(sizeof(float) * hsize);
        y = (int*)   malloc(sizeof(int)   * hsize);
        n = 0;
        for (sc = lowbound; sc <= highbound; sc++)
        {
            x[sc-lowbound] = (float) sc + 0.5; /* crude, but tests OK */
            y[sc-lowbound] = h->histogram[sc - h->min];
            n             += h->histogram[sc - h->min];
        }

        if (n < 100) goto FITFAILED;  /* require fitting to at least 100 points */

        /* If we're censoring, estimate z, the number of censored guys
         * left of the bound. Our initial estimate is crudely that we're
         * missing e^-1 of the total distribution (which would be exact
         * if we censored exactly at mu; but we censored at the observed peak).
         * Subsequent estimates are more exact based on our current estimate of mu.
         */
        if (censor)
        {
            if (iteration == 0)
                z = MIN(h->total-n, (int) (0.58198 * (float) n));
            else
            {
                double psx;
                psx = EVDDistribution((float) lowbound, mu, lambda);
                z = MIN(h->total-n, (int) ((double) n * psx / (1. - psx)));
            }
        }

        /* Do an ML fit
         */
        if (censor) {
            if (! EVDCensoredFit(x, y, hsize, z, (float) lowbound, &mu, &lambda))
                goto FITFAILED;
        } else
        if (! EVDMaxLikelyFit(x, y, hsize, &mu, &lambda))
            goto FITFAILED;

        /* Find the Eval = 1 point as a new highbound;
         * the total number of samples estimated to "belong" to the EVD is n+z
         */
        new_highbound = (int)
                (mu - (log (-1. * log((double) (n+z-1) / (double)(n+z))) / lambda));

        free(x);
        free(y);
        if (new_highbound >= highbound) break;
        highbound = new_highbound;
    }

    /* Set the histogram parameters;
     * - we fit from lowbound to highbound; thus we lose 2 degrees of freedom
     *   for fitting mu, lambda, but we get 1 back because we're unnormalized
     *   in this interval, hence we pass 2-1 = 1 as ndegrees.
     */
    ExtremeValueSetHistogram(h, mu, lambda, lowbound, highbound, 1);
    return 1;

    FITFAILED:
    UnfitHistogram(h);
    if (x != NULL) free(x);
    if (y != NULL) free(y);
    return 0;
}



/* Function: PrintASCIIHistogram()
 *
 * Purpose:  Print a "prettified" histogram to a file pointer.
 *           Deliberately a look-and-feel clone of Bill Pearson's
 *           excellent FASTA output.
 *
 * Args:     fp     - open file to print to (stdout works)
 *           h      - histogram to print
 */
void
PrintASCIIHistogram(FILE *fp, struct histogram_s *h)
{
    int units;
    int maxbar;
    int num;
    int i, idx;
    char buffer[81];		/* output line buffer */
    int  pos;			/* position in output line buffer */
    int  lowbound, lowcount;	/* cutoffs on the low side  */
    int  highbound, highcount;	/* cutoffs on the high side */
    int  emptybins = 3;

    /* Find out how we'll scale the histogram.
     * We have 59 characters to play with on a
     * standard 80-column terminal display:
     * leading "%5d %6d %6d|" occupies 20 chars.
     * Save the peak position, we'll use it later.
     */
    maxbar = 0;
    for (i = h->lowscore - h->min; i <= h->highscore - h->min; i++)
        if (h->histogram[i] > maxbar)
        {
            maxbar   = h->histogram[i];     /* max height    */
            lowbound = i + h->min;     	/* peak position */
        }

    /* Truncate histogram display on both sides, ad hoc fashion.
     * Start from the peak; then move out until we see <emptybins> empty bins,
     * and stop.
     */
    highbound = lowbound;		/* start at peak position */
    for (num = 0; lowbound > h->lowscore; lowbound--)
    {
        i = lowbound - h->min;
        if (h->histogram[i] > 0) { num = 0;               continue; } /* reset */
        if (++num == emptybins)  { lowbound += emptybins; break;    } /* stop  */
    }
    for (num = 0; highbound < h->highscore; highbound++)
    {
        i = highbound - h->min;
        if (h->histogram[i] > 0) { num = 0;                continue; } /* reset */
        if (++num == emptybins)  { highbound -= emptybins; break;    } /* stop  */
    }
    /* collect counts outside of bounds */
    for (lowcount = 0, i = h->lowscore - h->min; i <= lowbound - h->min; i++)
        lowcount += h->histogram[i];
    for (highcount = 0, i = h->highscore - h->min; i >= highbound - h->min; i--)
        highcount += h->histogram[i];

    /* maxbar might need raised now; then set our units  */
    if (lowcount  > maxbar) maxbar = lowcount;
    if (highcount > maxbar) maxbar = highcount;
    units = ((maxbar-1)/ 59) + 1;


    /* Print the histogram
     */
    fprintf(fp, "%5s %6s %6s  (one = represents %d sequences)\n",
            "score", "obs", "exp", units);
    fprintf(fp, "%5s %6s %6s\n", "-----", "---", "---");
    buffer[80] = '\0';
    buffer[79] = '\n';
    for (i = h->lowscore; i <= h->highscore; i++)
    {
        memset(buffer, ' ', 79 * sizeof(char));
        idx = i - h->min;

        /* Deal with special cases at edges
         */
        if      (i < lowbound)  continue;
        else if (i > highbound) continue;
        else if (i == lowbound && i != h->lowscore)
        {
            sprintf(buffer, "<%4d %6d %6s|", i+1, lowcount, "-");
            if (lowcount > 0) {
                num = 1+(lowcount-1) / units;
                if (num > 60) exit(-1);
                for (pos = 20; num > 0; num--)  buffer[pos++] = '=';
            }
            fputs(buffer, fp);
            continue;
        }
        else if (i == highbound && i != h->highscore)
        {
            sprintf(buffer, ">%4d %6d %6s|", i, highcount, "-");
            if (highcount > 0) {
                num = 1+(highcount-1) / units;
                for (pos = 20; num > 0; num--)  buffer[pos++] = '=';
            }
            fputs(buffer, fp);
            continue;
        }

        /* Deal with most cases
         */
        if (h->fit_type != HISTFIT_NONE)
            sprintf(buffer, "%5d %6d %6d|",
                    i, h->histogram[idx], (int) h->expect[idx]);
        else
            sprintf(buffer, "%5d %6d %6s|", i, h->histogram[idx], "-");
        buffer[20] = ' ';		/* sprintf writes a null char */

        /* Mark the histogram bar for observed hits
         */
        if (h->histogram[idx] > 0) {
            num = 1 + (h->histogram[idx]-1) / units;
            for (pos = 20; num > 0; num--)  buffer[pos++] = '=';
        }

        /* Mark the theoretically expected value
         */
        if (h->fit_type != HISTFIT_NONE && (int) h->expect[idx] > 0)
        {
            pos = 20 + (int)(h->expect[idx]-1) / units;
            if (pos >= 78) pos = 78; /* be careful of buffer bounds */
            buffer[pos] = '*';
        }

        /* Print the line
         */
        fputs(buffer, fp);
    }

    /* Print details about the statistics
     */
    switch (h->fit_type) {
        case HISTFIT_NONE:
            fprintf(fp, "\n\n%% No statistical fit available\n");
            break;

        case HISTFIT_EVD:
            fprintf(fp, "\n\n%% Statistical details of theoretical EVD fit:\n");
            fprintf(fp, "              mu = %10.4f\n", h->param[EVD_MU]);
            fprintf(fp, "          lambda = %10.4f\n", h->param[EVD_LAMBDA]);
            fprintf(fp, "chi-sq statistic = %10.4f\n", h->chisq);
            fprintf(fp, "  P(chi-square)  = %10.4g\n", h->chip);
            break;

        case HISTFIT_GAUSSIAN:
            fprintf(fp, "\n\n%% Statistical details of theoretical Gaussian fit:\n");
            fprintf(fp, "            mean = %10.4f\n", h->param[GAUSS_MEAN]);
            fprintf(fp, "              sd = %10.4f\n", h->param[GAUSS_SD]);
            fprintf(fp, "chi-sq statistic = %10.4f\n", h->chisq);
            fprintf(fp, "  P(chi-square)  = %10.4g\n", h->chip);
            break;
    }
    return;
}



/* Function: PrintXMGRHistogram()
 * Date:     SRE, Wed Nov 12 11:02:00 1997 [St. Louis]
 *
 * Purpose:  Print an XMGR data file that contains two data sets:
 *               - xy data for the observed histogram
 *               - xy data for the theoretical histogram
 */
void
PrintXMGRHistogram(FILE *fp, struct histogram_s *h)
{
    int sc;			/* integer score in histogram structure */
    double val;

    /* First data set is the observed histogram
     */
    for (sc = h->lowscore; sc <= h->highscore; sc++)
        if (h->histogram[sc - h->min] > 0)
            fprintf(fp, "%-6d %f\n", sc,
                    (float) h->histogram[sc - h->min]/ (float) h->total);
    fprintf(fp, "&\n");

    /* Second data set is the theoretical histogram
     */
    if (h->fit_type != HISTFIT_NONE)
    {
        for (sc = h->lowscore; sc <= h->highscore; sc++)
        {
            val =
                    (1. - ExtremeValueP((float)sc+1, h->param[EVD_MU], h->param[EVD_LAMBDA]))-
                    (1. - ExtremeValueP((float)sc, h->param[EVD_MU], h->param[EVD_LAMBDA]));
            fprintf(fp, "%-6d %f\n", sc, val);
        }
        fprintf(fp, "&\n");
    }
}

/* Function: PrintXMGRDistribution()
 * Date:     SRE, Wed Nov 12 11:02:09 1997 [St. Louis]
 *
 * Purpose:  Print an XMGR data file that contains two data sets:
 *               - xy data for the observed distribution P(S<x)
 *               - xy data for the theoretical distribution P(S<x)
 */
void
PrintXMGRDistribution(FILE *fp, struct histogram_s *h)
{
    int sc;			/* integer score in histogram structure */
    int cum;			/* cumulative count */
    double val;

    /* First data set is the observed distribution;
     * histogram bin x contains # of scores between x and x+1,
     * hence the sc+1 offset.
     */
    for (cum = 0, sc = h->lowscore; sc <= h->highscore; sc++)
    {
        cum += h->histogram[sc - h->min];
        fprintf(fp, "%-6d %f\n", sc + 1, (float) cum / (float) h->total);
    }
    fprintf(fp, "&\n");

    /* Second data set is the theoretical histogram
     */
    if (h->fit_type != HISTFIT_NONE)
    {
        for (sc = h->lowscore; sc <= h->highscore; sc++)
        {
            val = (1. - ExtremeValueP((float) sc, h->param[EVD_MU],
                                      h->param[EVD_LAMBDA]));
            fprintf(fp, "%-6d %f\n", sc, val);
        }
        fprintf(fp, "&\n");
    }
}

/* Function: PrintXMGRRegressionLine()
 * Date:     SRE, Wed Nov 12 11:02:19 1997 [St. Louis]
 *
 * Purpose:  Print an XMGR data file that contains two data sets:
 *               - xy data for log log transform of observed distribution P(S<x)
 *               - xy data for log log transform of theoretical distribution P(S<x)
 */
void
PrintXMGRRegressionLine(FILE *fp, struct histogram_s *h)
{
    int sc;			/* integer score in histogram structure */
    int cum;
    double val;			/* log log transform */

    /* First data set is the observed distribution;
     * histogram bin x contains # of scores between x and x+1,
     * hence the sc+1 offset.
     */
    for (cum = 0, sc = h->lowscore; sc <= h->highscore; sc++)
    {
        cum += h->histogram[sc - h->min];
        val = log (-1. * log((double) cum /  (double) h->total));
        if (cum < h->total)
            fprintf(fp, "%-6d %f\n", sc + 1, val);
    }
    fprintf(fp, "&\n");

    /* Second data set is the theoretical histogram
     */
    if (h->fit_type != HISTFIT_NONE)
    {
        for (sc = h->lowscore; sc <= h->highscore; sc++)
        {
            val = log(-1. * log(1. - ExtremeValueP((float) sc, h->param[EVD_MU],
                                                   h->param[EVD_LAMBDA])));
            fprintf(fp, "%-6d %f\n", sc, val);
        }
        fprintf(fp, "&\n");
    }
}


/* Function: Linefit()
 *
 * Purpose:  Given points x[0..N-1] and y[0..N-1], fit to
 *           a straight line y = a + bx.
 *           a, b, and the linear correlation coefficient r
 *           are filled in for return.
 *
 * Args:     x     - x values of data
 *           y     - y values of data
 *           N     - number of data points
 *           ret_a - RETURN: intercept
 *           ret_b - RETURN: slope
 *           ret_r - RETURN: correlation coefficient
 *
 * Return:   1 on success, 0 on failure.
 */
int
Linefit(float *x, float *y, int N, float *ret_a, float *ret_b, float *ret_r)
{
    float xavg, yavg;
    float sxx, syy, sxy;
    int   i;

    /* Calculate averages, xavg and yavg
     */
    xavg = yavg = 0.0;
    for (i = 0; i < N; i++)
    {
        xavg += x[i];
        yavg += y[i];
    }
    xavg /= (float) N;
    yavg /= (float) N;

    sxx = syy = sxy = 0.0;
    for (i = 0; i < N; i++)
    {
        sxx    += (x[i] - xavg) * (x[i] - xavg);
        syy    += (y[i] - yavg) * (y[i] - xavg);
        sxy    += (x[i] - xavg) * (y[i] - yavg);
    }
    *ret_b = sxy / sxx;
    *ret_a = yavg - xavg*(*ret_b);
    *ret_r = sxy / (sqrt(sxx) * sqrt(syy));
    return 1;
}

/* Function: EVDBasicFit()
 * Date:     SRE, Wed Nov 12 11:02:27 1997 [St. Louis]
 *
 * Purpose:  Fit a score histogram to the extreme value
 *           distribution. Set the parameters lambda
 *           and mu in the histogram structure. Fill in the
 *           expected values in the histogram. Calculate
 *           a chi-square test as a measure of goodness of fit.
 *
 *           This is the basic version of ExtremeValueFitHistogram(),
 *           in a nonrobust form: simple linear regression with no
 *           outlier pruning.
 *
 * Methods:  Uses a linear regression fitting method [Collins88,Lawless82]
 *
 * Args:     h         - histogram to fit
 *
 * Return:   (void)
 */
void
EVDBasicFit(struct histogram_s *h)
{
    float *d;            /* distribution P(S < x)          */
    float *x;            /* x-axis of P(S<x) for Linefit() */
    int    hsize;
    int    sum;
    int    sc, idx;		/* loop indices for score or score-h->min   */
    float  slope, intercept;	/* m,b fit from Linefit()                   */
    float  corr;			/* correlation coeff of line fit, not used  */
    float  lambda, mu;		/* slope, intercept converted to EVD params */

    /* Allocations for x, y axes
     * distribution d runs from min..max with indices 0..max-min
     *     i.e. score - min = index into d, x, histogram, and expect
     */
    hsize = h->highscore - h->lowscore + 1;
    d         = (float *) malloc(sizeof(float) * hsize);
    x         = (float *) malloc(sizeof(float) * hsize);
    for (idx = 0; idx < hsize; idx++)
        d[idx] = x[idx] = 0.;

    /* Calculate P(S < x) distribution from histogram.
     * note off-by-one of sc, because histogram bin contains scores between
     * x and x+1.
     */
    sum = 0;
    for (sc = h->lowscore; sc <= h->highscore; sc++)
    {
        sum += h->histogram[sc - h->min];
        d[sc - h->lowscore] = (float) sum / (float) h->total;
        x[sc - h->lowscore] = (float) (sc + 1);
    }

    /* Do a linear regression fit to the log[-log(P(S<x))] "line".
     * we have log[-log(1-P(S>x))]  = -lambda * x + lambda * mu
     * so lambda = -m  and mu = b/lambda
     */
    /* convert y axis to log[-log(P(S<x))]  */
    for (sc = h->lowscore; sc < h->highscore; sc++)
        d[sc - h->lowscore] = log(-1. * log(d[sc - h->lowscore]));

    /* do the linear regression */
    Linefit(x, d, hsize-1, &intercept, &slope, &corr);
    /* calc mu, lambda */
    lambda = -1. * slope;
    mu     = intercept / lambda;

    /* Set the EVD parameters in the histogram;
     * pass 2 for additional lost degrees of freedom because we fit mu, lambda.
     */
    ExtremeValueSetHistogram(h, mu, lambda, h->lowscore, h->highscore, 2);

    free(x);
    free(d);
    return;
}


/* Function: GaussianFitHistogram()
 *
 * Purpose:  Fit a score histogram to a Gaussian distribution.
 *           Set the parameters mean and sd in the histogram
 *           structure, as well as a chi-squared test for
 *           goodness of fit.
 *
 * Args:     h         - histogram to fit
 *           high_hint - score cutoff; above this are `real' hits that aren't fit
 *
 * Return:   1 if fit is judged to be valid.
 *           else 0 if fit is invalid (too few seqs.)
 */
int
GaussianFitHistogram(struct histogram_s *h, float high_hint)
{
    float sum;
    float sqsum;
    float delta;
    int   sc;
    int   nbins;
    int   hsize, idx;

    /* Clear any previous fitting from the histogram.
     */
    UnfitHistogram(h);

    /* Determine if we have enough hits to fit the histogram;
     * arbitrarily require 1000.
     */
    if (h->total < 1000) { h->fit_type = HISTFIT_NONE; return 0; }

    /* Simplest algorithm for mean and sd;
     * no outlier detection yet (not even using high_hint)
     *
     * Magic 0.5 correction is because our histogram is for
     * scores between x and x+1; we estimate the expectation
     * (roughly) as x + 0.5.
     */
    sum = sqsum = 0.;
    for (sc = h->lowscore; sc <= h->highscore; sc++)
    {
        delta  = (float) sc + 0.5;
        sum   += (float) h->histogram[sc-h->min] * delta;
        sqsum += (float) h->histogram[sc-h->min] * delta * delta;
    }
    h->fit_type          = HISTFIT_GAUSSIAN;
    h->param[GAUSS_MEAN] = sum / (float) h->total;
    h->param[GAUSS_SD]   = sqrt((sqsum - (sum*sum/(float)h->total)) /
                                (float)(h->total-1));

    /* Calculate the expected values for the histogram.
     * Note that the magic 0.5 correction appears again.
     * Calculating difference between distribution functions for Gaussian
     * would be correct but hard.
     */
    hsize     = h->max - h->min + 1;
    h->expect = (float *) malloc(sizeof(float) * hsize);
    for (idx = 0; idx < hsize; idx++)
        h->expect[idx] = 0.;

    for (sc = h->min; sc <= h->max; sc++)
    {
        delta = (float) sc + 0.5 - h->param[GAUSS_MEAN];
        h->expect[sc - h->min] =
                (float) h->total * ((1. / (h->param[GAUSS_SD] * sqrt(2.*3.14159))) *
                                    (exp(-1.* delta*delta / (2. * h->param[GAUSS_SD] * h->param[GAUSS_SD]))));
    }

    /* Calculate the goodness-of-fit (within region that was fitted)
     */
    h->chisq = 0.;
    nbins    = 0;
    for (sc = h->lowscore; sc <= h->highscore; sc++)
        if (h->expect[sc-h->min] >= 5. && h->histogram[sc-h->min] >= 5)
        {
            delta = (float) h->histogram[sc-h->min] - h->expect[sc-h->min];
            h->chisq += delta * delta / h->expect[sc-h->min];
            nbins++;
        }
    /* -1 d.f. for normalization; -2 d.f. for two free parameters */
    if (nbins > 3)
        h->chip = (float) IncompleteGamma((double)(nbins-3)/2.,
                                          (double) h->chisq/2.);
    else
        h->chip = 0.;

    return 1;
}


/* Function: GaussianSetHistogram()
 *
 * Purpose:  Instead of fitting the histogram to a Gaussian,
 *           simply set the Gaussian parameters from an external source.
 */
void
GaussianSetHistogram(struct histogram_s *h, float mean, float sd)
{
    int   sc;
    int   hsize, idx;
    int   nbins;
    float delta;

    UnfitHistogram(h);
    h->fit_type          = HISTFIT_GAUSSIAN;
    h->param[GAUSS_MEAN] = mean;
    h->param[GAUSS_SD]   = sd;

    /* Calculate the expected values for the histogram.
     */
    hsize     = h->max - h->min + 1;
    h->expect = (float *) malloc(sizeof(float) * hsize);
    for (idx = 0; idx < hsize; idx++)
        h->expect[idx] = 0.;

    /* Note: ideally we'd use the Gaussian distribution function
     * to find the histogram occupancy in the window sc..sc+1.
     * However, the distribution function is hard to calculate.
     * Instead, estimate the histogram by taking the density at sc+0.5.
     */
    for (sc = h->min; sc <= h->max; sc++)
    {
        delta = ((float)sc + 0.5) - h->param[GAUSS_MEAN];
        h->expect[sc - h->min] =
                (float) h->total * ((1. / (h->param[GAUSS_SD] * sqrt(2.*3.14159))) *
                                    (exp(-1.*delta*delta / (2. * h->param[GAUSS_SD] * h->param[GAUSS_SD]))));
    }

    /* Calculate the goodness-of-fit (within whole region)
     */
    h->chisq = 0.;
    nbins    = 0;
    for (sc = h->lowscore; sc <= h->highscore; sc++)
        if (h->expect[sc-h->min] >= 5. && h->histogram[sc-h->min] >= 5)
        {
            delta = (float) h->histogram[sc-h->min] - h->expect[sc-h->min];
            h->chisq += delta * delta / h->expect[sc-h->min];
            nbins++;
        }
    /* -1 d.f. for normalization */
    if (nbins > 1)
        h->chip = (float) IncompleteGamma((double)(nbins-1)/2.,
                                          (double) h->chisq/2.);
    else
        h->chip = 0.;
}


// read comma seperated string
std::pair<float *, int> ReadScores(const char * str){
    int n = 0;
    char * pch;
    char * str_copy = new char[strlen(str)+1];
    strcpy(str_copy, str);
    pch = strtok(str_copy, ",");
    while (pch != NULL)
    {
        n++;
        pch = strtok(NULL, ",");
    }
    float * scores = new float[n];
    strcpy(str_copy, str);
    pch = strtok(str_copy, ",");
    for (int i=0; i<n; i++)
    {
        scores[i] = atof(pch);
        pch = strtok(NULL, ",");
    }
    delete [] str_copy;
    return std::make_pair(scores, n);
}


int main(int argc, char *argv[])
{
    if (argc<2)
    {
        printf("Usage: %s <scores file>\n",argv[0]);
        exit(1);
    }
    std::pair<float *, int> ret = ReadScores(argv[1]);
    float mu, lambda;
    EVDMaxLikelyFit(ret.first, NULL, ret.second, &mu, &lambda);

    printf("lamda=%f   mu=%f\n",lambda, mu);
    return 0;
}


