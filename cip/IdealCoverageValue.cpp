#include "IdealCoverageValue.h"
#include "RandomNumberGenerator.h"
#include <math.h>
#include <cmath>



IdealCoverageValue::IdealCoverageValue(double xmean, double eta, double xsd, float phi, const std::vector<double>& d):
	_xmean(xmean),
	_eta(eta),
	_xsd(xsd),
	_phi(phi),
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
	double dataVar = 0.0;
	for (size_t i = 0; i < _data.size(); ++i)
	{
		dataVar += (_data[i] - dataMean) * (_data[i] - dataMean);
	}
	dataVar /= _data.size();

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
