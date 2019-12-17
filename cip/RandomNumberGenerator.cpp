#include "RandomNumberGenerator.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

RandomNumberGenerator::RandomNumberGenerator()
{
}


RandomNumberGenerator::~RandomNumberGenerator()
{
}

/*******************************************************************************
**  Function:  tinv(double u, double df, double delta, int *ifault, int *i)   **
**                                                                            **
**  Authors:   Bruce Schmeise and Huifen Chen, June 1992                      **
**                                                                            **
**  Purpose:   Compute the inverse of the cdf of noncentral t dist.           **
**                                                                            **
**  Input:                                                                    **
**      u    : cdf of noncentral t dist.  (value)                             **
**      df   : degree of freedom of noncentral t dist.   (value)              **
**      delta: noncentral parameter of noncentral t dist.  (value)            **
**                                                                            **
**  Output:                                                                   **
**      tinv_r: the inverse of cdf                                            **
**      ifault: error indicator      (point)                                  **
**      i     : number of iterations (point)                                  **
**                                                                            **
**  Notice: In the support routines, the original APPLIED STATISTICS          **
**         code is in UPPER CASE.  Our modifications are in lower             **
**         case.                                                              **
*******************************************************************************/
double RandomNumberGenerator::tinv(double u, double df, double delta, int *ifault, int *i)
{
	//double ppnd(double, int *);
	//void RT(double, double, double, double *, double *, int *);
	int  itrmax, lower, upper, twobd, ibound;
	double tinv_r, rn, sq, funct, t_l, t_u, dt, tnc, denst, error;

	itrmax = 1000;
	lower = 1;
	upper = 2;
	twobd = 3;
	*i = 0;
	*ifault = 0;

	/*
	**.....get the initial guess.  reference is "Continuous univariate
	**       distributions-2" by johnson and kotz, p 207.
	*/
	tinv_r = 0.;
	rn = ppnd(u, ifault);
	sq = 1. + .5*((delta*delta) - (rn*rn)) / df;
	if (sq >= 0.)
		tinv_r = (delta + rn * sqrt(sq)) / (1. - .5 * (rn*rn) / df);
	if (u > .5 && tinv_r < delta)
		tinv_r = delta;

	RT(tinv_r, df, delta, &tnc, &denst, ifault);
	funct = tnc - u;

	/*
	**....if the function value is close to zero, return.
	*/
	if (fabs((double)funct) < 1.E-6) return(tinv_r);

	/*
	**.....define the bound indicator, ibound.
	**        ibound = lower   ==> only lower bound exists
	**               = upper   ==> only upper bound exists
	**               = twobd   ==> have lower and upper bounds
	*/
	if (funct < 0.) {
		ibound = lower;
		t_l = tinv_r;
	}
	else {
		ibound = upper;
		t_u = tinv_r;
	}

	dt = funct / denst;
	tinv_r = tinv_r - dt;

	for ((*i) = 1; (*i) <= itrmax; (*i)++) {
		RT(tinv_r, df, delta, &tnc, &denst, ifault);
		funct = tnc - u;

		/*
		**........update ibound, lower bound and upper bound
		*/
		if (funct < 0.) {
			if (ibound == upper) {
				ibound = twobd;
				t_l = tinv_r;
			}
			if (tinv_r > t_l) t_l = tinv_r;
			tinv_r = t_l;
		}
		else {
			if (ibound == lower) {
				ibound = twobd;
				t_u = tinv_r;
			}
			if (tinv_r <= t_u)  t_u = tinv_r;
			tinv_r = t_u;
		}

		/*
		**......determine which method to use, either newton or bisection
		*/
		if (ibound == twobd && ((denst*(tinv_r - t_l) - tnc)*
			(denst*(tinv_r - t_u) - tnc) >= 0.)) {
			/*.....if there are lower and upper bounds and next newton
			**     iterate lies outside the bounds, use bisection method.
			*/
			dt = (t_u - t_l) / 2.;
			tinv_r = (t_u + t_l) / 2.;
		}
		else {
			/*...........use newton's method */
			dt = funct / denst;
			tinv_r = tinv_r - dt;
		}

		error = dt / tinv_r;

		if (fabs(error) < 1.E-6 || fabs(funct) < 1.E-6) return(tinv_r);
	}

	return(tinv_r);
}



/*******************************************************************************
**  huifen chen and bruce schmeiser, june 1992                                **
**  update: january 1994.                                                     **
**          correct the computation of density function.                      **
**  purpose: evaluate the cdf and density of non-central t dist.              **
**  input:                                                                    **
**      t    : the point to evaluate (value)                                  **
**      df   : degree of freedom     (value)                                  **
**      delta: non-central parameter (value)                                  **
**  output:                                                                   **
**      tnc: cumulative distribution at point t (point)                       **
**      denst: the density at point t           (point)                       **
**      ifault: the error indicator             (point)                       **
**  reference: ALGORITHM AS 243  APPL. STATIST. (1989),VOL.38,NO. 1           **
*******************************************************************************/
void RandomNumberGenerator::RT(double T, double DF, double DELTA, double *TNC, double *DENST, int *IFAULT)
{
	//double gammln(double);
	//double betai(double, double, double *, double);
	//double anordf(double);
	//void down_recur(double *, double, double *, double *, double *,
	//	double *, double *, double, double *, double *, double *,
	//	double *, double, double *, double *, double, double);
	//void up_recur(double *, double, double *, double *, double *,
	//	double *, double *, double, double *, double *, double *,
	//	double *, double, double *, double *, double, double);
	double A = 0.0, ALBETA = 0.0, B = 0.0, DEL = 0.0, EN = 0.0, ERRMAX = 0.0, PI = 0.0, TMP = 0.0, AMODE = 0.0,
		GEVEN = 0.0, GODD = 0.0, HALF = 0.0, ITRMAX = 0.0, LAMBDA = 0.0, ONE = 0.0, P = 0.0, PP = 0.0, Q = 0.0, QQ = 0.0,
		TT = 0.0, TWO = 0.0, X = 0.0, XEVEN = 0.0, XODD = 0.0, ZERO = 0.0, RT_R = 0.0, ALN_HALF = 0.0, SQRT_TWO = 0.0,
		GAMMLN_A = 0.0, GAMMLN_J1 = 0.0, GAMMLN_B = 0.0, ALBETA2 = 0.0, XX = 0.0, BETA = 0.0, BETA2 = 0.0,
		AA = 0.0, BB = 0.0, EEN = 0.0, XXEVEN = 0.0, XXODD = 0.0, GGEVEN = 0.0, GGODD = 0.0, BBETA = 0.0, BBETA2 = 0.0,
		ITR_LEFT = 0.0;
	(void)RT_R;
	double jplus = 0.;
	int   NEGDEL, J;

	/*
	**  Note - ITRMAX and ERRMAX may be changed to suit one's needs.
	*/
	ITRMAX = 5000.1;
	ERRMAX = 1.E-07;
	ZERO = 0.0;
	HALF = 0.5;
	ONE = 1.0;
	TWO = 2.0;
	ALN_HALF = -.693147181;
	SQRT_TWO = 1.414213562;
	PI = 3.141592654;

	/*
	**.....initialize the summation value
	*/
	*TNC = ZERO;
	*DENST = ZERO;

	/*
	**.....if the degree of freedom is not positive return
	*/
	*IFAULT = 2;
	if (DF <= ZERO)
		return;
	*IFAULT = 0;

	/*
	**.....when t is zero
	*/
	if (T == 0.) {
		/* if delta is less than tweleve, anordf is not negligible; */
		if (DELTA < 12.) *TNC = anordf(-DELTA);
		TMP = -0.5*DELTA*DELTA + gammln((DF + 1.) / 2) - gammln(DF / 2.);
		*DENST = exp(TMP) / sqrt(DF * PI);
		return;
	}

	TT = T;
	DEL = DELTA;
	NEGDEL = 0;

	if (T < ZERO) {
		NEGDEL = 1;
		TT = -TT;
		DEL = -DEL;
	}

	/*
	**     Initialize twin series (Guenther, J. Statist. Computn. Simuln.
	**     vol.6, 199, 1978).
	*/
	X = T * T / (T* T + DF);
	LAMBDA = DEL * DEL;

	/*.....start from the mode of poisson
	*/
	AMODE = LAMBDA / TWO;
	J = (int)AMODE;
	EN = J;
	A = J + HALF;
	B = HALF * DF;
	jplus = (double)(J + 1.);
	XODD = betai(A, B, &GAMMLN_A, X);
	XEVEN = betai(jplus, B, &GAMMLN_J1, X);
	if (DEL == 0.)
		P = .5;
	else
		P = exp(ALN_HALF - AMODE + J * log(AMODE) - GAMMLN_J1);

	Q = (P*DEL / (SQRT_TWO*A)) * exp(GAMMLN_J1 - GAMMLN_A);
	*TNC = P * XODD + Q * XEVEN;

	GAMMLN_B = gammln(B);

	/*.....compute the tmp = logarithm of x**a * (1.-x)**b
	*/
	TMP = 2.* A * log(TT) + B * log(DF) - (A + B) * log(TT*TT + DF);
	ALBETA = gammln(A + B) - gammln(A + 1.) - GAMMLN_B;
	GODD = exp(ALBETA + TMP);
	ALBETA2 = gammln(A + B + HALF) - gammln(A + 1.5) - GAMMLN_B;
	GEVEN = exp(ALBETA2 + TMP) * sqrt(X);

	/*.....compute xx=log(x*(1-x)).  note 1-x = df / (tt**2+df)
	*/
	XX = 2.*log(TT) + log(DF) - 2. * log((TT*TT) + DF);

	/*.....compute beta = godd*a/(x*(1-x)); beta2 = geven*(a+.5)/(x*(1-x)).
	*/
	BETA = exp(ALBETA + TMP + log(A) - XX);
	BETA2 = exp(ALBETA2 + TMP + log(A + HALF) - XX);

	*DENST = P * BETA + Q * BETA2;

	if (J > 0) {
		/*........sum the items left to jth item */
		AA = A;
		BB = B;
		EEN = EN;
		XXODD = XODD;
		XXEVEN = XEVEN;
		GGODD = GODD;
		GGEVEN = GEVEN;
		PP = P;
		QQ = Q;
		BBETA = BETA;
		BBETA2 = BETA2;
		down_recur(&AA, BB, &EEN, &XXODD, &XXEVEN, &GGODD, &GGEVEN, X, &PP, &QQ,
			&BBETA, &BBETA2, AMODE, TNC, DENST, ERRMAX, ITRMAX);
	}

	/*
	**.....sum the items right to jth item
	*/
	ITR_LEFT = ITRMAX - (AMODE - EEN);
	up_recur(&A, B, &EN, &XODD, &XEVEN, &GODD, &GEVEN, X, &P, &Q, &BETA, &BETA2,
		AMODE, TNC, DENST, ERRMAX, ITR_LEFT);

	*IFAULT = 0;
	EN = EN - EEN;
	if (EN > ITRMAX) (*IFAULT) = 1;

	*DENST = *DENST * 2.* (1. - X) * sqrt((X - (X*X)) / DF);

	/*
	**.....if del is less than tweleve, anordf is not negligible
	*/
	if (DEL < 12.) *TNC = *TNC + anordf(-DEL);

	if (NEGDEL) *TNC = ONE - *TNC;
	if (*TNC > 1.) *TNC = 1.;
	if (*TNC < 0.) *TNC = 0.;
	return;
}


/******************************************************************************
** FUNCTION: up_recur(double *a, double b, double *en, double *xodd,          **
**               double *xeven, double *godd, double *geven, double x,        **
**               double *p, double *q, double *beta, double *beta2,           **
**               double amode, double *tnc, double *denst,                    **
**               double errmax, double itrmax)                                **
**                                                                            **
** Purpose: Compute part of the cdf and density of noncentral t               **
**          distribution.  Since both cdf and density are an infinite         **
**          sum, this function approximate the RIGHT-hand-side of the         **
**          sum starting from the mode (amode).                               **
*******************************************************************************/
void RandomNumberGenerator::up_recur(double *a, double b, double *en, double *xodd,
	double *xeven, double *godd, double *geven, double x, double *p,
	double *q, double *beta, double *beta2, double amode, double *tnc,
	double *denst, double errmax, double itrmax)
{
	double half, one, c1, c2, rerror, add;
	half = 0.5;
	one = 1.0;

	/*
	**     Repeat until convergence
	*/
	do {
		*en = (*en) + one;
		*xodd = (*xodd) - (*godd);
		*xeven = (*xeven) - (*geven);
		*p = (*p) * amode / (*en);
		*q = (*q) * amode / ((*en) + half);
		add = (*p) * (*xodd) + (*q) * (*xeven);
		*tnc = *tnc + add;
		*a = (*a) + one;
		c1 = x * ((*a) + b - one) / (*a);
		c2 = x * ((*a) + b - half) / ((*a) + half);
		*godd = (*godd) * c1;
		*geven = (*geven) * c2;
		*beta = (*beta) * ((*a) + b - one) / ((*a) - one);
		*beta2 = (*beta2) * ((*a) + b - half) / ((*a) - half);
		*denst = (*denst) + (*p) * (*beta) + (*q) * (*beta2);
		rerror = fabs(add / (*tnc));
	} while ((rerror > errmax) && (((*en) - amode) <= itrmax));

	return;
}


/*******************************************************************************
** FUNCTION: down_recur(double *a, double b, double *en, double *xodd,        **
**    double *xeven, double *godd, double *geven, double x, double *p,        **
**    double *q, double *beta, double *beta2, double amode,                   **
**    double *tnc, double *denst, double errmax, double itrmax)               **
**                                                                            **
** Purpose: Compute part of the cdf and density of noncentral t               **
**          distribution.  Since both cdf and density are an infinite         **
**          sum, this function approximate the LEFT-hand-side of the          **
**          sum from the firt item up to the mode (amode).                    **
*******************************************************************************/
void RandomNumberGenerator::down_recur(double *a, double b, double *en, double *xodd,
	double *xeven, double *godd, double *geven, double x, double *p, double *q,
	double *beta, double *beta2, double amode, double *tnc, double *denst,
	double errmax, double itrmax)
{
	double half, one, c1, c2, add, rerror;
	half = 0.5;
	one = 1.0;

	/*
	**     Repeat until convergence
	*/
	do {
		c1 = (*a) / (((*a) + b - 1.) * x);
		c2 = ((*a) + half) / (((*a) + b - half) * x);
		*godd = (*godd) * c1;
		*geven = (*geven) * c2;
		*xodd = (*xodd) + (*godd);
		*xeven = (*xeven) + (*geven);
		*p = (*p) * (*en) / amode;
		*q = (*q) * ((*en) + half) / amode;
		add = (*p) * (*xodd) + (*q) * (*xeven);
		*tnc = (*tnc) + add;
		*beta = (*beta) * ((*a) - one) / (((*a) + b - 1.) * x);
		*beta2 = (*beta2) * ((*a) - half) / (((*a) + b - half) * x);
		*denst = (*denst) + (*p) * (*beta) + (*q) * (*beta2);
		rerror = fabs(add / (*tnc));
		*en = (*en) - one;
		*a = (*a) - one;
	} while (((*en) > 0.) && (rerror >= errmax)
		&& ((amode - (*en)) <= itrmax));
	return;
}


/*******************************************************************************
** FUNCTION: anordf(double z)                                                 **
**                                                                            **
** Purpose:  standard normal cdf                                              **
**           ibm scientific subroutine package, 1967, page 78.                **
**                                                                            **
** Input:                                                                     **
**   z:      standard normal random variate                                   **
**                                                                            **
** Output:                                                                    **
** anordf_r: return standard  normal random variate                           **
*******************************************************************************/
double RandomNumberGenerator::anordf(double z)
{
	double az, t, d, p, anordf_r;
	az = fabs(z);
	t = 1. / (1. + .2316419*az);
	d = .3989423 * exp((double)-z * z*.5);
	p = 1. - d * t*((((1.330274*t - 1.821256)*t + 1.781478)*t
		- .3565638)*t + .3193815);
	if (z < 0.)  p = 1. - p;
	anordf_r = p;
	return(anordf_r);
}


/*******************************************************************************
**  FUNCTION:  gammln(double xx)                                              **
**             revised from "Numerical Recipes" by huifen, may, 1992          **
**                                                                            **
**  Purpose:   returns the logarithm of gamma function at xx.  the            **
**             internal arithmetic will be done in double precision.          **
**             if xx is less than one, the reflection formula (6.1.4)         **
**             in "Numerical Recipes" is used.                                **
**                                                                            **
**  Reference: function gammln in "Numerical Recipes" by Press,               **
**             Flannery, Teukolsky and Vetterling                             **
*******************************************************************************/
double RandomNumberGenerator::gammln(double xx)
{
	int   j;
	double cof[6], stp, half, one, fpf, x, tmp, ser, gammln_r, tt;

	cof[0] = 76.18009173;
	cof[1] = -86.50532033;
	cof[2] = 24.01409822;
	cof[3] = -1.231739516;
	cof[4] = .120858003e-2;
	cof[5] = -.536382e-5;
	stp = 2.50662827465;
	half = 0.5;
	one = 1.0;
	fpf = 5.5;
	gammln_r = 0.;

	if ((xx == 1.) || (xx == 2.)) return (gammln_r);

	x = xx - one;
	if (xx < one) x = -x;
	tmp = x + fpf;
	tmp = (x + half)*log(tmp) - tmp;
	ser = one;
	for (j = 0; j <= 5; j++) {
		x = x + one;
		ser = ser + cof[j] / x;
	}
	gammln_r = tmp + log(stp * ser);
	if (xx < one) {
		tmp = 3.141592654*(one - xx);
		tt = sin(tmp);
		tmp = fabs(tmp / tt);
		gammln_r = log(tmp) - gammln_r;
	}
	return (gammln_r);
}


/*******************************************************************************
**  FUNCTION: betai (double a, double b, double *gammln_a, double x)          **
**                                                                            **
**  Purpose:  returns the incomplete beta(a,b) function, ie the cdf           **
**            of beta(a,b) at x.                                              **
**  From   :  "Numberical Recipes" by Press, Flannery, Teukolsky              **
**            and Vetterling.                                                 **
*******************************************************************************/
double RandomNumberGenerator::betai(double a, double b, double *gammln_a, double x)
{
	//double betacf(double, double, double);
	//double gammln(double);
	double betai_r, bt;
	if ((x < 0.) || (x > 1.))
		printf("\n bad argument x in betai; x is between 0. and 1.");
	*gammln_a = gammln(a);
	if ((x == 0.) || (x == 1))
		bt = 0.;
	else {
		bt = exp(gammln(a + b) - (*gammln_a) - gammln(b) +
			a * log(x) + b * log(1. - x));
	}

	if (x < ((a + 1.) / (a + b + 2.))) {
		betai_r = bt * betacf(a, b, x) / a;
		return (betai_r);
	}
	else {
		betai_r = 1. - bt * betacf(b, a, (1. - x)) / b;
		return (betai_r);
	}
}


/*******************************************************************************
**  FUNCTION: betacf(double a, double b, double x)                            **
**                                                                            **
**  Purpose:  Returns continued fraction for incomplete beta function,        **
**            used by BETAI                                                   **
**  From   :  "Numberical Recipes" by Press, Flannery, Teukolsky              **
**            and Vetterling.                                                 **
*******************************************************************************/
double RandomNumberGenerator::betacf(double a, double b, double x)
{
	int    itmax, m, em, tem;
	double eps, am, bm, az, qab, qap, qam, bz, betacf_r, d, ap, bp,
		app, bpp, aold;

	am = 1.;
	bm = 1.;
	az = 1.;
	qab = a + b;
	qap = a + 1.;
	qam = a - 1.;
	bz = 1. - qab * x / qap;
	itmax = 100;
	eps = 3.e-7;

	for (m = 1; m <= itmax; m++) {
		em = m;
		tem = em + em;
		d = em * (b - m)*x / ((qam + tem)*(a + tem));
		ap = az + d * am;
		bp = bz + d * bm;
		d = -(a + em)*(qab + em)*x / ((a + tem)*(qap + tem));
		app = ap + d * az;
		bpp = bp + d * bz;
		aold = az;
		am = ap / bpp;
		bm = bp / bpp;
		az = app / bpp;
		bz = 1.;
		if (fabs(az - aold) < (eps*fabs(az))) {
			betacf_r = az;
			return(betacf_r);
		}
	}
	printf("\n a or b too big, or itmax too small for function betacf");

	betacf_r = az;
	return(betacf_r);
}


/*******************************************************************************
**  FUNCTION: ppnd(double p, int *ifault)                                     **
**                                                                            **
**  Purpose:  Compute the inversed distribution function of a                 **
**            standard normal distribution that is set to the                 **
**            value of z satisfying P(Z<z) = p. therefore,                    **
**            0 < p < 1. ifault = 0 means that all is well.                   **
**                                                                            **
**  Input:                                                                    **
**            p:      the percentage point to evaluate                        **
**                                                                            **
**  Output:                                                                   **
**            ppnd:   the inverse cdf                                         **
**            ifault: error indicator                                         **
*******************************************************************************/
double RandomNumberGenerator::ppnd(double p, int *ifault)
{
	double ppnd1,
		q,
		zero = 0.,
		split = 0.42,
		half = 0.5,
		one = 1.;
	double r,
		a0 = 2.50662823884,
		a1 = -18.61500062529,
		a2 = 41.39119773534,
		a3 = -25.44106049637,
		b1 = -8.47351093090,
		b2 = 23.08336743743,
		b3 = -21.06224101826,
		b4 = 3.13082909833,
		c0 = -2.78718931138,
		c1 = -2.29796479134,
		c2 = 4.85014127135,
		c3 = 2.32121276858,
		d1 = 3.54388924762,
		d2 = 1.63706781897;

	*ifault = 0;
	q = p - half;
	if (fabs(q) > split) {
		r = p;
		if (q > zero) r = one - p;
		if (r <= zero) {
			*ifault = 1;
			ppnd1 = zero;
			return(ppnd1);
		}
		else {
			r = sqrt(-log(r));
			ppnd1 = (((c3 * r + c2) * r + c1) * r + c0);
			ppnd1 = ppnd1 / ((d2 * r + d1) * r + one);
			if (q < 0) ppnd1 = -ppnd1;
			return(ppnd1);
		}
	}
	else {
		r = q * q;
		ppnd1 = q * (((a3 * r + a2) * r + a1) * r + a0);
		ppnd1 = ppnd1 / ((((b4 * r + b3) * r + b2) * r + b1) * r + one);
		return(ppnd1);
	}
}


/******************************************************************************/
/*  purpose:                                                                  */
/*      generate random numbers                                               */
/*                                                                            */
/*  input variables:                                                          */
/*      *iseed: random number seed                                            */
/*                                                                            */
/*  output:                                                                   */
/*      *iseed: new random number seed                                        */
/*      random number                                                         */
/*                                                                            */
/******************************************************************************/
double RandomNumberGenerator::u16807d(int *iseed)
{
	*iseed = (int)fmod(*iseed * 16807., 2147483647.);
	return (*iseed / 2147483647.);
}


/******************************************************************************/
/*  purpose:                                                                  */
/*      compute gamma cdf with alpha = a and beta (scale) evalued at          */
/*      y = x*beta                                                            */
/*  input variables:                                                          */
/*      a: alpha                                                              */
/*      x: cdf evalued at y where y = x*beta                                  */
/*  output:                                                                   */
/*      the probability that gamma is less than or equal to y= x*beta         */
/*  note:                                                                     */
/*      to evaluate gamma(a, beta) cdf at y                                   */
/*      call gamma(a, y/beta) to get the desired probability                  */
/******************************************************************************/
double RandomNumberGenerator::gammp(double a, double x)
{
	double gp;
	//double gcf(double, double);
	//double gser(double, double);

	if (x < 0. || a <= 0.) printf("invalid argument in routine gammp\n");
	if (x < (a + 1.)) {
		gp = gser(a, x);
		return gp;
	}
	else {
		gp = gcf(a, x);
		return (1. - gp);
	}
}


double RandomNumberGenerator::gser(double a, double x)
{
	int n;
	double sum, del, ap, gln, gammser;
	//double gammln(double);

	gln = gammln(a);
	if (x <= 0.) {
		if (x < 0.) printf("x less than 0 in routine gser\n");
		gammser = 0.;
		return gammser;
	}
	else {
		ap = a;
		del = sum = 1. / a;
		for (n = 1; n <= 100; n++) {
			++ap;
			del *= x / ap;
			sum += del;
			if (fabs(del) < fabs(sum)*(3.0e-7)) {
				gammser = sum * exp(-x + a * log(x) - gln);
				return gammser;
			}
		}
		printf("a too large, 100 too small in routine gser\n");
		exit(1);
	}
}


double RandomNumberGenerator::gcf(double a, double x)
{
	int i;
	double an, b, c, d, del, h, gln, gammcf;
	//double gammln(double);

	gln = gammln(a);
	b = x + 1. - a;
	c = 1. / (1.0e-30);
	d = 1. / b;
	h = d;

	for (i = 1; i <= 100; i++) {
		an = -i * (i - a);
		b += 2.;
		d = an * d + b;
		if (fabs(d) < 1.0e-30) d = 1.0e-30;
		c = b + an / c;
		if (fabs(c) < 1.0e-30) c = 1.0e-30;
		d = 1. / d;
		del = d * c;
		h *= del;
		if (fabs(del - 1.) < 1.0e-30) break;
	}
	if (i > 100) printf(" a too large, 100 too small in gcf\n");
	gammcf = exp(-x + a * log(x) - gln)*h;

	return gammcf;
}

/******************************************************************************/
/*  Purpose: Genereate a random variate from gamma distributions,             */
/*             denoted by Gamma(alpha, beta), whose shape and scale           */
/*             parameters are alpha and beta, respectively.  The              */
/*               density function f(x) is proportional to                     */
/*                       x**(alpha-1) * exp{-x/beta}.                         */
/*    Author: Huifen Chen, Industrial Engineering, Chung-Yuan Christian Univ. */
/*            TAIWAN.  Email: huifen@cycu.edu.tw                              */
/*    Input:                                                                  */
/*        u: probability.                                                     */
/*        alpha: shape parameter, positive.                                   */
/*        beta:  scale parameter, positive.                                   */
/*    Output:                                                                 */
/*        rgamma: the quantile of gamma(alpha, beta).                         */
/*        ifault: error indicator:                                            */
/*                0 = no error,                                               */
/*                1 = infeasible parameter,                                   */
/*                2 = error from chi-square inverse.                          */
/*    Routine called:                                                         */
/*        gammln:  computing the logarithm of gamma function                  */
/*        ppchi2:  computing percentage point of chi-square distributions     */
/******************************************************************************/
double RandomNumberGenerator::gammainv(double u, double alpha, double beta, int *ifault)
{
	int   ierror;
	double v, g, chi2, x;
	//double gammln(double);
	//double ppchi2(double, double, double, int *);

	*ifault = 0;
	/* ...check the feasibility of alpha and beta. */
	if ((alpha <= 0.) || (beta <= 0.))  *ifault = 1;

	/* ...special case: exponential distribution */
	if (alpha == 1.) {
		x = beta * (-log((double)1. - u));
		return(x);
	}

	/* ...compute the gamma(alpha, beta) inverse.                   *
	 *    ...compute the chi-square inverse with 2*alpha degrees of *
	 *       freedom, which is equivalent to gamma(alpha, 2).       */
	v = 2.0 * alpha;
	g = gammln(alpha);

	chi2 = ppchi2(u, v, g, &ierror);
	if (ierror != 0)  *ifault = 2;

	/* ...transfer chi-square to gamma. */
	x = beta * chi2 / 2.0;
	return(x);
}


/******************************************************************************/
/*  Algorithm AS 91   Appl. Statist. (1975) Vol.24, P.35                      */
/*    To evaluate the percentage points of the chi-squared                    */
/*    probability distribution function.                                      */
/*    p must lie in the range 0.000002 to 0.999998, /$ percentage $/          */
/*    v must be positive,   /$ degrees of freedom $/                          */
/*    g must be supplied and should be equal to                               */
/*      ln(gamma(v/2.0))                                                      */
/*                                                                            */
/*    Incorporates the suggested changes in AS R85 (vol.40(1),                */
/*    pp.233-5, 1991) which should eliminate the need for the limited         */
/*    range for p above, though these limits have not been removed            */
/*    from the routine.                                                       */
/*    If IFAULT = 4 is returned, the result is probably as accurate as        */
/*    the machine will allow.                                                 */
/*                                                                            */
/*    Auxiliary routines required: PPND = AS 111 (or AS 241) and              */
/*                                 GAMMAD = AS 239.                           */
/******************************************************************************/
double RandomNumberGenerator::ppchi2(double p, double v, double g, int *ifault)

{
	int    i, maxit = 500, if1;
	double aa, e, ppch, zero, half, one,
		two, three, six, pmin, pmax, c1, c2, c3, c4, c5, c6, c7,
		c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19,
		c20, c21, c22, c23, c24, c25, c26, c27, c28, c29, c30,
		c31, c32, c33, c34, c35, c36, c37, c38, a, b, c, ch, p1, p2,
		q, s1, s2, s3, s4, s5, s6, t, x, xx;
	//double ppnd(double, int *);
	//double gammad(double, double, int *);


	aa = 0.6931471806;
	e = 0.0000005;
	pmin = 0.000002;
	pmax = 0.999998;
	zero = 0.0;
	half = 0.5;
	one = 1.0;
	two = 2.0;
	three = 3.0;
	six = 6.0;
	c1 = 0.01;
	c2 = 0.222222;
	c3 = 0.32;
	c4 = 0.4;
	c5 = 1.24;
	c6 = 2.2;
	c7 = 4.67;
	c8 = 6.66;
	c9 = 6.73;
	c10 = 13.32;
	c11 = 60.0;
	c12 = 70.0;
	c13 = 84.0;
	c14 = 105.0;
	c15 = 120.0;
	c16 = 127.0;
	c17 = 140.0;
	c18 = 1175.0;
	c19 = 210.0;
	c20 = 252.0;
	c21 = 2264.0;
	c22 = 294.0;
	c23 = 346.0;
	c24 = 420.0;
	c25 = 462.0;
	c26 = 606.0;
	c27 = 672.0;
	c28 = 707.0;
	c29 = 735.0;
	c30 = 889.0;
	c31 = 932.0;
	c32 = 966.0;
	c33 = 1141.0;
	c34 = 1182.0;
	c35 = 1278.0;
	c36 = 1740.0;
	c37 = 2520.0;
	c38 = 5040.0;

	/* ...test arguments and initialise  */
	ppch = -one;
	*ifault = 1;
	if ((p < pmin) || (p > pmax))
		return(ppch);
	*ifault = 2;
	if (v <= zero)
		return(ppch);
	*ifault = 0;
	xx = half * v;
	c = xx - one;

	/* ...starting approximation for small chi-squared */

	if (v >= (-c5 * log(p))) {
		/*....starting approximation for v less than or equal to 0.32 */
		if (v > c3) {
			/* call to algorithm AS 111 - note that p has been tested above.
			 AS 241 could be used as an alternative.   */
			x = ppnd(p, &if1);

			/* starting approximation using Wilson and Hilferty estimate */

			p1 = c2 / v;
			ch = v * pow((x * sqrt(p1) + one - p1), 3.0);

			/* starting approximation for p tending to 1  */

			if (ch > (c6 * v + six))
				ch = -two * (log(one - p) - c * log(half * ch) + g);
		}
		else {
			ch = c4;
			a = log(one - p);
			do {
				q = ch;
				p1 = one + ch * (c7 + ch);
				p2 = ch * (c9 + ch * (c8 + ch));
				t = -half + (c7 + two * ch) / p1 - (c9 + ch * (c10 +
					three * ch)) / p2;
				ch = ch - (one - exp(a + g + half * ch + c * aa) *
					p2 / p1) / t;
			} while (fabs((double)q / ch - one) > c1);
		}
	}
	else {
		ch = pow((p * xx * exp(g + xx * aa)), (one / xx));
		if (ch < e) {
			ppch = ch;
			return(ppch);
		}
	}

	/*....call to algorithm AS 239 and calculation of seven term
		  Taylor series    */
	for (i = 1; i <= maxit; i++) {
		q = ch;
		p1 = half * ch;
		p2 = p - gammad(p1, xx, &if1);
		if (if1 != 0) {
			*ifault = 3;
			return(ch);
		}
		else {
			t = p2 * exp(xx * aa + g + p1 - c * log(ch));
			b = t / ch;
			a = half * t - b * c;
			s1 = (c19 + a * (c17 + a * (c14 + a * (c13 + a * (c12 +
				c11 * a))))) / c24;
			s2 = (c24 + a * (c29 + a * (c32 + a * (c33 + c35 *
				a)))) / c37;
			s3 = (c19 + a * (c25 + a * (c28 + c31 * a))) / c37;
			s4 = (c20 + a * (c27 + c34 * a) + c * (c22 + a * (c30 +
				c36 * a))) / c38;
			s5 = (c13 + c21 * a + c * (c18 + c26 * a)) / c37;
			s6 = (c15 + c * (c23 + c16 * c)) / c38;
			ch = ch + t * (one + half * t * s1 - b * c * (s1 - b *
				(s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))));
			if (fabs((double)q / ch - one) > e) {
				ppch = ch;
				return(ppch);
			}
		}
	}
	*ifault = 4;
	ppch = ch;
	return(ppch);
}


/******************************************************************************/
/*  Revision: replace auxiliary functions ALNGAM by gammln and                */
/*              ALNORM by ndtr, since ALNORM is written in                    */
/*              FORTRAN 66.  huifen, 11/1/94.                                 */
/*                                                                            */
/*    ALGORITHM AS239  APPL. STATIST. (1988) VOL. 37, NO. 3                   */
/*    Computation of the Incomplete Gamma Integral                            */
/*    /$ integral of function f(t)= t**(p-1) * exp{-t} in [0,x],              */
/*         where x >= 0 and p > 0.    fen, 11/1/94            $/              */
/*                                                                            */
/*    Auxiliary functions required: ALNGAM = logarithm of the gamma           */
/*    function (AS245), and ALNORM = algorithm AS66                           */
/******************************************************************************/
double RandomNumberGenerator::gammad(double X, double P, int *ifault)
{
	double PN1, PN2, PN3, PN4, PN5, PN6, TOL, OFLO,
		XBIG, ARG, C, RN, A, B, ONE, ZERO, GAMMADV,
		AN, TWO, ELIMIT, PLIMIT, THREE, NINE, d;
	double min;
	//void ndtr(double, double *, double *);
	//double gammln(double);

	ZERO = 0.0;
	ONE = 1.0;
	TWO = 2.0;
	OFLO = pow(10.0, 35.0);  /* Original code: OFLO = pow(10.0,37.0);    */
	THREE = 3.0;             /* has overflow problem.  fen 3/31/2000     */
	NINE = 9.0;
	TOL = pow(10.0, -14.0);
	XBIG = pow(10.0, 8.0);
	PLIMIT = 1000.0;
	ELIMIT = -88.0;
	GAMMADV = ZERO;

	/*... Check that we have valid values for X and P  */

	if ((P <= ZERO) || (X < ZERO)) {
		*ifault = 1;
		return(GAMMADV);
	}
	*ifault = 0;
	if (X == ZERO)
		return(GAMMADV);

	/*... Use a normal approximation if P > PLIMIT    */

	if (P > PLIMIT) {
		PN1 = THREE * sqrt(P) * (pow((X / P), (ONE / THREE)) + ONE /
			(NINE * P) - ONE);
		ndtr(PN1, &P, &d);
		return(P);
	}

	/*... If X is extremely large compared to P then set GAMMAD = 1  */
	if (X > XBIG) {
		GAMMADV = ONE;
		return(GAMMADV);
	}
	if ((X <= ONE) || (X < P)) {
		/*....Use Pearson's series expansion.
			  (Note that P is not large enough to force overflow in gammln).
			  No need to test IFAULT on exit since P > 0.   */
		ARG = P * log(X) - X - gammln(P + ONE);
		C = ONE;
		GAMMADV = ONE;
		A = P;
		do {
			A = A + ONE;
			C = C * X / A;
			GAMMADV = GAMMADV + C;
		} while (C > TOL);
		ARG = ARG + log(GAMMADV);
		GAMMADV = ZERO;
		if (ARG >= ELIMIT) GAMMADV = exp(ARG);
	}
	else {
		/*.... Use a continued fraction expansion */
		ARG = P * log(X) - X - gammln(P);
		A = ONE - P;
		B = A + X + ONE;
		C = ZERO;
		PN1 = ONE;
		PN2 = X;
		PN3 = X + ONE;
		PN4 = X * B;
		GAMMADV = PN3 / PN4;
		while (1) {
			A = A + ONE;
			B = B + TWO;
			C = C + ONE;
			AN = A * C;
			PN5 = B * PN3 - AN * PN1;
			PN6 = B * PN4 - AN * PN2;
			if (fabs((double)PN6) > ZERO) {
				RN = PN5 / PN6;
				min = (TOL < (TOL * RN)) ? TOL : (TOL * RN);
				if (fabs((double)GAMMADV - RN) <= min)
					break;
				GAMMADV = RN;
			}
			PN1 = PN3;
			PN2 = PN4;
			PN3 = PN5;
			PN4 = PN6;
			if (fabs((double)PN5) >= OFLO) {
				/*.... Re-scale terms in continued fraction if terms are large */
				PN1 = PN1 / OFLO;
				PN2 = PN2 / OFLO;
				PN3 = PN3 / OFLO;
				PN4 = PN4 / OFLO;
			}
		}
		ARG = ARG + log(GAMMADV);
		GAMMADV = ONE;
		if (ARG >= ELIMIT)
			GAMMADV = ONE - exp(ARG);
	}
	return(GAMMADV);
}

/******************************************************************************/
/*     compute the cumulative distribution function P(X <= x)                 */
/*       of a standard normal random variable X.                              */
/*     input                                                                  */
/*       x: normal value                                                      */
/*     output                                                                 */
/*       p: P(X <= x)                                                         */
/*       d: standard normal density at x                                      */
/*     reference                                                              */
/*       c. hastings, approximations for digital computers.                   */
/*       princeton university press, princeton, new jersey, usa, 1955.        */
/*       this code is a minor modification, by bruce schmeiser,               */
/*       of the code in the ibm scientific subroutine package,                */
/*       1967, page 78.                                                       */
/******************************************************************************/
void RandomNumberGenerator::ndtr(double x, double *p, double *d)
{
	double t, ax;

	ax = fabs((double)x);
	t = 1.0 / (ax * .2316419 + 1.0);
	*d = exp(-(x * x * .5)) * .3989423;
	*p = 1.0 - *d * t * ((((t * 1.330274 - 1.821256) * t
		+ 1.781478) * t - .3565638) * t + .3193815);
	if (x < 0.0) *p = 1 - *p;
	return;
}

/******************************************************************************/
/*  purpose:                                                                  */
/*      generate observations of a steady-state AR(1) process                 */
/*  input:                                                                    */
/*      i : current observation number                                        */
/*      *iseed: random number seed                                            */
/*      xmean: mean of AR(1) process                                          */
/*      xsd: papameter of AR(1) process                                       */
/*      phi: parameter of AR(1) process                                       */
/*  output:                                                                   */
/*      *x: new observation                                                   */
/******************************************************************************/
double RandomNumberGenerator::ar1(int i, int *iseed, double xmean, double xsd, double phi, double *x)
{
	int ifault;
	double p, z, c;
	//double u16807d(int *);
	//double ppnd(double, int *);

	p = u16807d(iseed);
	z = ppnd(p, &ifault);
	if (xsd < 0.) {
		xsd = -xsd;
	}
	if (i == 0) {
		*x = xmean + xsd * z;
	}
	else {
		c = 1. - phi * phi;
		*x = xmean + phi * ((*x) - xmean) + xsd * sqrt(c)*z;
	}

	return *x;
}

//generate one observation from a AR(1) time series with normal marginal distribution
double RandomNumberGenerator::generator(double xmean, double xsd, float phi, double x, int *iseed)
{
	double c, z, p;
	int ifault;

	p = u16807d(iseed);
	z = ppnd(p, &ifault);

	if (xsd < 0)
	{
		xsd = -xsd;
	}

	if (fabs(x - xmean) > 10. * xsd)
	{
		x = xmean + (xsd * z);
	}
	else
	{
		c = 1. - phi * phi;
		if (c < 0.)
		{
			if (phi <= -1.)
			{
				x = xmean - (x - xmean);
			}
		}
		else
		{
			x = xmean + phi * (x - xmean) + (xsd * sqrt(c) * z);
		}
	}
	return (x);
}

