#pragma once
#include <vector>

class ASAP3
{
public:
	ASAP3(double xsd, float phi, int iseed);
	~ASAP3();

	void procedure(double alpha);

private:
	std::vector<double> generate(int n);

	double _x;
	double _xsd;
	float _phi;
	int _iseed;
	std::vector<double> _data;
};

