#pragma once
#include <vector>
#include "RandomNumberGenerator.h"

class ASAP3
{
public:
	ASAP3(const RandomNumberGenerator::Parameter& p);
	~ASAP3();

	void procedure(bool RelPrec, double alpha, double r_star);

	double get_CIlb() const { return _CIlb; }
	double get_CIub() const { return _CIub; }
	double get_data_mean() const { return _data_mean;  }
	double get_xmean() const { return _xmean;  }
	const std::vector<double>& get_data() const { return _data; }
	const std::vector<double>& get_observation() const { return _observation; }
private:
	std::vector<double> generate(int n);
	void SWILK(bool& INIT, const std::vector<double>& X, int N, int N1, int N2, std::vector<double>& A, double& W, double& PW, int& IFAULT) const;

	RandomNumberGenerator::Parameter _parameter;
	std::vector<double> _data;
	std::vector<double> _observation;
	double _CIlb;
	double _CIub;
	double _data_mean;
	double _xmean;
};

