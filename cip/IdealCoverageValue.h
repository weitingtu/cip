#pragma once
#include <vector>

class IdealCoverageValue
{
public:
	IdealCoverageValue(double data_mean, double eta, double xsd, float phi, const std::vector<double>& o);
	~IdealCoverageValue();

	double run() const;
	double run_mm1() const;

private:
	double _normalCDF(double value) const;

	double _data_mean;
	double _eta;
	double _xsd;
	float _phi;
	std::vector<double> _observation;
};

