#include "BisectionSearch.h"
#include "Skart.h"
#include "ASAP3.h"
#include <math.h>

BisectionSearch::BisectionSearch(double tolerance):
	_a(0.0),
	_b(1.0),
	_eta(0.5),
	_psi(0.0),
	_data_mean(0.0),
	_data(),
	_observation(),
	_tolerance(tolerance)
{
}


BisectionSearch::~BisectionSearch()
{
}

double BisectionSearch::run_skart(const RandomNumberGenerator::Parameter& p, bool precReq, double alpha, double hrstar)
{
	_eta = 1 - alpha;
	while (true)
	{
    	Skart s(p);
		s.skart_procedure( precReq, alpha, hrstar );

		double CIlb = s.get_CIlb(); // L
		double CIub = s.get_CIub(); // U
		double theta = s.get_data_mean();
		int C = 0;
		if (theta >= CIlb && theta <= CIub)
		{
			C = 1;
		}
		else
		{
			C = 0;
		}

		if (1 == C)
		{
			_b = _eta;
		}
		else
		{
			_a = _eta;
		}
		_eta = (_a + _b) / 2;
		if ((_b - _a) <= _tolerance)
		{
			_psi = _eta;
			_data_mean = s.get_data_mean();
			_data = s.get_data();
			_observation = s.get_observation();
			printf("psi = %f\n", _psi);
			return _psi;
		}
	}
}

double BisectionSearch::run_asap3( const RandomNumberGenerator::Parameter& p, bool RelPrec, double alpha, double r_star )
{
	_eta = 1 - alpha;
	while (true)
	{
    	ASAP3 s(p);
		s.procedure( RelPrec, alpha, r_star );

		double CIlb = s.get_CIlb(); // L
		double CIub = s.get_CIub(); // U
		double theta = s.get_data_mean();
		int C = 0;
		if (theta >= CIlb && theta <= CIub)
		{
			C = 1;
		}
		else
		{
			C = 0;
		}

		if (1 == C)
		{
			_b = _eta;
		}
		else
		{
			_a = _eta;
		}
		_eta = (_a + _b) / 2;
		if ((_b - _a) <= _tolerance)
		{
			_psi = _eta;
			_data_mean = s.get_data_mean();
			_data = s.get_data();
			_observation = s.get_observation();
			printf("psi = %f\n", _psi);
			return _psi;
		}
	}
}
