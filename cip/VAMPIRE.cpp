#include "VAMPIRE.h"
#include "Skart.h"
#include "BisectionSearch.h"
#include "IdealCoverageValue.h"
#include <vector>
#include <math.h>

VAMPIRE::VAMPIRE()
{
}


VAMPIRE::~VAMPIRE()
{
}

double VAMPIRE::run(CIP_TYPE cip_type, 
	const RandomNumberGenerator::Parameter& parameter,
	double alpha,
	bool precReq,
	double hrstar,
	bool RelPrec,
	double r_star,
	double tolerance)
{
	size_t size = 400;
	std::vector<int> iseeds(size, 0);
	for (size_t i = 0; i < size; ++i)
	{
		iseeds.at(i) = i + 1;
	}

	std::vector<double> psi_v(size, 0.0);
	std::vector<double> delta_v(size, 0.0);
	for (size_t i = 0; i < iseeds.size(); ++i)
	{
		int iseed = iseeds[i];
		RandomNumberGenerator::Parameter p = parameter;
		p.ar1.iseed = iseed;
		p.mm1.iseed = iseed;
		BisectionSearch b(tolerance);
		double psi = 0.0;
		if(CIP_TYPE::SKART == cip_type)
		{ 
     		psi = b.run_skart(parameter, precReq, alpha, hrstar );
		}
		else
		{
		    psi = b.run_asap3(parameter, RelPrec, alpha, r_star );
		}

		double delta = 0.0;
		const RandomNumberGenerator::AR1_Parameter& ar1 = parameter.ar1;
		IdealCoverageValue icv(b.get_data_mean(), b.get_eta(), ar1.xsd, ar1.phi, b.get_data());
		if (RandomNumberGenerator::Type::AR1 == parameter.type)
		{
		    delta = icv.run();
		}
		else if (RandomNumberGenerator::Type::MM1 == parameter.type)
		{ 
		    delta = icv.run_mm1();
		}

		psi_v.at(i) = psi;
		delta_v.at(i) = delta;
	}

	double criteria = 0.0;
	for (size_t i = 0; i < psi_v.size(); ++i)
	{
		criteria += pow(psi_v.at(i) - delta_v.at(i), 2);
	}
	criteria /= psi_v.size();

	return criteria;
}
