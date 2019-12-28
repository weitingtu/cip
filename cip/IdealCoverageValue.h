#pragma once
#include <vector>

class IdealCoverageValue
{
public:
	IdealCoverageValue(double data_mean, double eta, double xsd, float phi, double arate, double srate, const std::vector<double>& o);
	~IdealCoverageValue();

	double run() const;
	double run_mm1() const;

private:
	double _normalCDF(double value) const;

	double _xmean;
	double _eta;
	double _xsd;
	float _phi;
	double _arate;
	double _srate;
	std::vector<double> _data;
};

