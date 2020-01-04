#pragma once
#include <vector>

class IdealCoverageValue
{
public:
	IdealCoverageValue(double data_mean, double eta, double xsd, float phi, double arate, double srate, const std::vector<double>& o);
	~IdealCoverageValue();

	double run();
	double run_mm1();

	double run(size_t n);
	double run_mm1(size_t n);

	double get_ideal_half_length() const { return _ideal_half_length;  }
private:
	double _normalCDF(double value) const;

	double _xmean;
	double _eta;
	double _xsd;
	float _phi;
	double _arate;
	double _srate;
	std::vector<double> _data;
	double _ideal_half_length;
};

