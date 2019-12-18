#pragma once
#include <vector>
#include "RandomNumberGenerator.h"

class BisectionSearch
{
public:
	BisectionSearch(double tolerance);
	~BisectionSearch();

	double run_skart(const RandomNumberGenerator::Parameter& p, bool precReq, double alpha, double hrstar);
	double run_asap3(const RandomNumberGenerator::Parameter& p, bool RelPrec, double alpha, double r_star );
	double get_eta() const { return _eta; }
	double get_data_mean() const { return _data_mean;  }
	const std::vector<double>& get_observation() const { return _observation; }

private:
	double _a;
	double _b;
	double _eta; // 1 - alpha
	double _psi;
	double _data_mean;
	std::vector<double> _observation;
	const double _tolerance;
};
