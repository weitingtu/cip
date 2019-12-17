#pragma once
class RandomNumberGenerator
{
public:
	RandomNumberGenerator();
	~RandomNumberGenerator();

	static double ppnd(double p, int *ifault);
	static double tinv(double u, double df, double delta, int *ifault, int *i);
	static void RT(double T, double DF, double DELTA, double *TNC, double *DENST, int *IFAULT);
	static void up_recur(double *a, double b, double *en, double *xodd,
		double *xeven, double *godd, double *geven, double x, double *p,
		double *q, double *beta, double *beta2, double amode, double *tnc,
		double *denst, double errmax, double itrmax);
	static void down_recur(double *a, double b, double *en, double *xodd,
		double *xeven, double *godd, double *geven, double x, double *p, double *q,
		double *beta, double *beta2, double amode, double *tnc, double *denst,
		double errmax, double itrmax);
	static double anordf(double z);
	static double gammln(double xx);
	static double betai(double a, double b, double *gammln_a, double x);
	static double betacf(double a, double b, double x);
	static double u16807d(int *iseed);
	static double gammp(double a, double x);
	static double gser(double a, double x);
	static double gcf(double a, double x);
	static double gammainv(double u, double alpha, double beta, int *ifault);
	static double ppchi2(double p, double v, double g, int *ifault);
	static double gammad(double X, double P, int *ifault);
	static void ndtr(double x, double *p, double *d);
	static double ar1(int i, int *iseed, double xmean, double xsd, double phi, double *x);

	static double generator(double xmean, double xsd, float phi, double x, int *iseed);
};

