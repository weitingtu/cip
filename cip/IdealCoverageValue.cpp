#include "IdealCoverageValue.h"
#include "RandomNumberGenerator.h"
#include <math.h>
#include <cmath>



IdealCoverageValue::IdealCoverageValue(double data_mean, double eta, double xsd, float phi, const std::vector<double>& o):
	_data_mean(data_mean),
	_eta(eta),
	_xsd(xsd),
	_phi(phi),
	_observation(o)
{
}


IdealCoverageValue::~IdealCoverageValue()
{
}

double IdealCoverageValue::run() const
{
	size_t n = _observation.size();
	double a = ((n - 1) / _phi - n + pow(_phi, n - 1)) / pow(1 / _phi - 1, 2);
	double sampleMeanVar = pow(_xsd, 2 ) / ((1 - pow(_phi, 2 )) * (n * n)) * (n + 2 * a);
	double sampleMeanSte = sqrt(sampleMeanVar);

	double p = (1 + _eta) / 2;
	int ifault = 0;
	double idealHalfWidth = RandomNumberGenerator::ppnd( p, &ifault) * sampleMeanSte;
	double observationMean = 0.0;
	for (size_t i = 0; i < _observation.size(); ++i)
	{
		observationMean += _observation[i];
	}
	observationMean /= _observation.size();

	double delta = 2 * _normalCDF( abs(observationMean - _data_mean) / sampleMeanSte ) - 1;

	return delta;
}

double IdealCoverageValue::run_mm1() const
{
	double observationMean = 0.0;
	for (size_t i = 0; i < _observation.size(); ++i)
	{
		observationMean += _observation[i];
	}
	observationMean /= _observation.size();
	double observationVar = 0.0;
	for (size_t i = 0; i < _observation.size(); ++i)
	{
		observationVar += (_observation[i] - observationMean) * (_observation[i] - observationMean);
	}
	observationVar /= _observation.size();

	double sampleMeanSte = sqrt(observationVar);

	double p = (1 + _eta) / 2;
	int ifault = 0;
	double idealHalfWidth = RandomNumberGenerator::ppnd( p, &ifault) * sampleMeanSte;

	double delta = 2 * _normalCDF( abs(observationMean - _data_mean) / sampleMeanSte ) - 1;

	return delta;
}

double IdealCoverageValue::_normalCDF(double value) const
{
	return 0.5 * erfc(-value * sqrt(0.5));
}
