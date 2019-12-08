#include "VAMPIRE.h"
#include "Skart.h"
#include "BisectionSearch.h"
#include "IdealCoverageValue.h"
#include <vector>



VAMPIRE::VAMPIRE()
{
}


VAMPIRE::~VAMPIRE()
{
}

double VAMPIRE::run()
{
	//size_t size = 400;
	size_t size = 1;
	std::vector<int> iseeds(size, 0);
	for (size_t i = 0; i < size; ++i)
	{
		iseeds.at(i) = i + 1;
	}

	double xsd = 2.1;
	float phi = 0.8f;
	
	std::vector<double> psi_v(size, 0.0);
	std::vector<double> delta_v(size, 0.0);
	for (size_t i = 0; i < iseeds.size(); ++i)
	{
		int iseed = iseeds[i];
		BisectionSearch b;
		double psi = b.run_skart(xsd, phi, iseed);

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
