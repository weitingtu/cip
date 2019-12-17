#pragma once
#include <vector>

class BisectionSearch
{
public:
	BisectionSearch();
	~BisectionSearch();

	double run_skart(double xmean, double xsd, float phi, int iseed);
	double run_asap3(double xmean, double xsd, float phi, int iseed);
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
