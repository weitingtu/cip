#include "IdealCoverageValue.h"
#include "RandomNumberGenerator.h"
#include <math.h>
#include <cmath>

int mm1var(int n, double tau, double vxbar, int *ier, double *rho,
	double *vx, double *arate, double *srate)
	/******************************************************************/
	/*  purpose: for m/m/1 queue                                      */
	/*  input:                                                        */
	/*    n: no. of observations (waiting time)                       */
	/*    tau: traffic intensity (arrival rate / service rate)        */
	/*    vxbar: variance of x bar                                    */
	/*  output:                                                       */
	/*    rho: correlation of waiting time of lag h, h=1,...,n-1      */
	/*    vx: variance of x                                           */
	/*    arate: arrival rate                                         */
	/*    srate: service rate                                         */
	/*  method  : use Simpson's rule to calculate the first n         */
	/*            autocorrelations                                    */
	/*  formula : D.J. Daley, on page 697, J. Austral. Math. Soc.     */
	/*            1968.                                               */
	/*            divide the integral on equation (34) on p. 697 into */
	/*            two integrals, left and right. Then apply Simpson's */
	/*            rule.                                               */
	/******************************************************************/
{
	int nint1 = 1000;
	int nint2 = 10000;
	int ih, inter;
	double fract = .95;
	double coef, a, f, t, realn, width;
	double coefw2, coefw4;
	double tstart, con;

	/* .....edit input parameters */
	*ier = 0;
	if (n < 1) {
		*ier = -1;
	}
	if (tau <= 0. || tau >= 1.) {
		*ier = -2;
	}
	if (vxbar <= 0.) {
		*ier = -3;
	}
	if (*ier != 0) {
		return 0;
	}
	/* .....initialize */
	a = 4. * tau / (1. + tau) / (1. + tau);
	for (ih = 1; ih <= n; ih++) {
		rho[ih] = 0.;
	}
	/* .....the left part of the interval (0,a) */
	width = fract * a / nint1;
	coefw2 = width * 2. / 3.;
	coefw4 = width * 4. / 3.;
	coef = coefw2;
	for (inter = 1; inter <= nint1 - 1; inter++) {
		t = width * inter;
		f = t * sqrt(t * (a - t)) / pow(1. - t, 3.);
		if (coef == coefw4) {
			coef = coefw2;
		}
		else {
			coef = coefw4;
		}
		for (ih = 1; ih <= n; ih++) {
			f = f * t;
			rho[ih] = rho[ih] + coef * f;
		}
	}
	/* .....the last point in the left part (the first part) */
	t = width * nint1;
	f = t * sqrt(t * (a - t)) / pow(1. - t, 3.);
	for (ih = 1; ih <= n; ih++) {
		f = f * t;
		rho[ih] = rho[ih] + width * f / 3.;
	}
	/* .....the first point in the right part ( the second part) */
	tstart = t;
	width = (1. - fract) * a / nint2;
	f = t * sqrt(t * (a - t)) / pow(1. - t, 3.);
	for (ih = 1; ih <= n; ih++) {
		f = f * t;
		rho[ih] = rho[ih] + width * f / 3.;
	}
	/* .....the right part of the interval (0,a) */
	coefw2 = width * 2. / 3.;
	coefw4 = width * 4. / 3.;
	coef = coefw2;
	for (inter = 1; inter <= nint2 - 1; inter++) {
		t = tstart + width * inter;
		f = t * sqrt(t * (a - t)) / pow(1. - t, 3.);
		if (coef == coefw4) {
			coef = coefw2;
		}
		else {
			coef = coefw4;
		}
		for (ih = 1; ih <= n; ih++) {
			f = f * t;
			rho[ih] = rho[ih] + coef * f;
		}
	}
	con = pow(1. - tau, 3.) * (1. + tau) /
		(6.28318531 * pow(tau, 3.) * (2. - tau));
	*vx = 0.;
	realn = (float)n;
	for (ih = 1; ih <= n; ih++) {
		rho[ih] = rho[ih] * con;
		*vx = *vx + (1. - ih / realn) * rho[ih];
	}
	//*vx = n * vxbar / (*vx * 2. + 1.);
	*arate = tau * sqrt(tau*(2. - tau) / *vx) / (1. - tau);
	*srate = *arate / tau;
	return 0;
}


IdealCoverageValue::IdealCoverageValue(double xmean, double eta, double xsd, float phi,
	double arate, double srate,
	const std::vector<double>& d):
	_xmean(xmean),
	_eta(eta),
	_xsd(xsd),
	_phi(phi),
	_arate(arate),
	_srate(srate),
	_data(d)
{
}


IdealCoverageValue::~IdealCoverageValue()
{
}

double IdealCoverageValue::run() const
{
	size_t n = _data.size();
	double a = ((n - 1) / _phi - n + pow(_phi, n - 1)) / pow(1 / _phi - 1, 2);
	double sampleMeanVar = pow(_xsd, 2 ) / ((1 - pow(_phi, 2 )) * (n * n)) * (n + 2 * a);
	double sampleMeanSte = sqrt(sampleMeanVar);

	double p = (1 + _eta) / 2;
	int ifault = 0;
	double idealHalfWidth = RandomNumberGenerator::ppnd( p, &ifault) * sampleMeanSte;
	double dataMean = 0.0;
	for (size_t i = 0; i < _data.size(); ++i)
	{
		dataMean += _data[i];
	}
	dataMean /= _data.size();

	double delta = 2 * _normalCDF( abs(dataMean - _xmean) / sampleMeanSte ) - 1;

	return delta;
}

double IdealCoverageValue::run_mm1() const
{
	double dataMean = 0.0;
	for (size_t i = 0; i < _data.size(); ++i)
	{
		dataMean += _data[i];
	}
	dataMean /= _data.size();

	double tau   = _arate / _srate;
	double alpha = _arate;
	double nu    = _srate;
	double varT = pow(tau, 3) * (2 - tau) / (pow(alpha, 2) * pow(1 - tau, 2)) + 1 / (nu * nu);
	//int mm1var(int n, float tau, float vxbar, int *ier, float *rho, float *vx, float *arate, float *srate);
	// mm1var input
	int n = _data.size();
	double vxbar = 1.0;
	// mm1var output
	int ier     = 0; 
	double *rho = new double[n+1];
	double vx    = 0.0;
	double arate = 0.0;
	double srate = 0.0;
	mm1var( n, tau, vxbar, &ier, rho, &vx, &arate, &srate);
	delete [] rho;

	double dataVar = 1.0 / n * varT * (1 + 2 * vx);
	//double dataVar = 0.0;
	//for (size_t i = 0; i < _data.size(); ++i)
	//{
	//	dataVar += (_data[i] - dataMean) * (_data[i] - dataMean);
	//}
	//dataVar /= _data.size();

	double sampleMeanSte = sqrt(dataVar);

	double p = (1 + _eta) / 2;
	int ifault = 0;
	double idealHalfWidth = RandomNumberGenerator::ppnd( p, &ifault) * sampleMeanSte;

	double delta = 2 * _normalCDF( abs(dataMean - _xmean) / sampleMeanSte ) - 1;

	return delta;
}

double IdealCoverageValue::_normalCDF(double value) const
{
	return 0.5 * erfc(-value * sqrt(0.5));
}
