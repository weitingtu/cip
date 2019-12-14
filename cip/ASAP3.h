#pragma once
#include <vector>

class ASAP3
{
public:
	ASAP3(double xsd, float phi, int iseed);
	~ASAP3();

	void procedure(double alpha);

	double get_CIlb() const { return _CIlb; }
	double get_CIub() const { return _CIub; }
	double get_data_mean() const { return _data_mean;  }
	const std::vector<double>& get_data() const { return _data; }
	const std::vector<double>& get_observation() const { return _observation; }
private:
	std::vector<double> generate(int n);
	void SWILK(bool& INIT, const std::vector<double>& X, int N, int N1, int N2, std::vector<double>& A, double& W, double& PW, int& IFAULT) const;

	double _x;
	double _xmean;
	double _xsd;
	float _phi;
	int _iseed;
	std::vector<double> _data;
	std::vector<double> _observation;
	double _CIlb;
	double _CIub;
	double _data_mean;
};

