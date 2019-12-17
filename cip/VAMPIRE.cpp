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

double VAMPIRE::run(CIP_TYPE cip_type)
{
	size_t size = 400;
	std::vector<int> iseeds(size, 0);
	for (size_t i = 0; i < size; ++i)
	{
		iseeds.at(i) = i + 1;
	}

	double xmean = 100.0;
	double xsd = 2.1;
	float phi = 0.8f;
	
	std::vector<double> psi_v(size, 0.0);
	std::vector<double> delta_v(size, 0.0);
	for (size_t i = 0; i < iseeds.size(); ++i)
	{
		int iseed = iseeds[i];
		BisectionSearch b;
		double psi = 0.0;
		if(CIP_TYPE::SKART == cip_type)
		{ 
     		psi = b.run_skart(xmean, xsd, phi, iseed);
		}
		else
		{
		    psi = b.run_asap3(xmean, xsd, phi, iseed);
		}

		IdealCoverageValue icv(b.get_data_mean() ,b.get_eta(), xsd, phi, b.get_observation());
		double delta = icv.run();

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
