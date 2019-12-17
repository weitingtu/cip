#include "BisectionSearch.h"
#include "Skart.h"
#include "ASAP3.h"
#include <math.h>

BisectionSearch::BisectionSearch():
	_a(0.0),
	_b(1.0),
	_eta(0.5),
	_psi(0.0),
	_data_mean(0.0),
	_observation(),
	_tolerance(0.0015)
{
}


BisectionSearch::~BisectionSearch()
{
}

double BisectionSearch::run_skart(double xmean, double xsd, float phi, int iseed)
{
	while (true)
	{
    	Skart s(xmean, xsd, phi, iseed);
		s.skart_procedure( 1 - _eta );

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
			_observation = s.get_observation();
			printf("psi = %f\n", _psi);
			return _psi;
		}
	}
}

double BisectionSearch::run_asap3(double xmean, double xsd, float phi, int iseed)
{
	while (true)
	{
    	ASAP3 s(xmean, xsd, phi, iseed);
		s.procedure( 1 - _eta );

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
			_observation = s.get_observation();
			printf("psi = %f\n", _psi);
			return _psi;
		}
	}
}
