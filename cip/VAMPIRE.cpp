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
     		psi = b.run_skart(p, precReq, alpha, hrstar );
		}
		else
		{
		    psi = b.run_asap3(p, RelPrec, alpha, r_star );
		}

	    double dataMean = 0.0;
    	for (size_t i = 0; i < b.get_data().size(); ++i)
    	{
	    	dataMean += b.get_data()[i];
	    }
	    dataMean /= b.get_data().size();

		double delta = 0.0;
		const RandomNumberGenerator::AR1_Parameter& ar1 = p.ar1;
		const RandomNumberGenerator::MM1_Parameter& mm1 = p.mm1;
		int A = 0;
		int B = b.get_data().size();
		int n = B / 2;
		double r = 0.15;
		while (true)
		{
		    IdealCoverageValue icv(b.get_xmean(), 1 - alpha, ar1.xsd, ar1.phi, mm1.arate, mm1.srate, b.get_data());
			if (RandomNumberGenerator::Type::AR1 == p.type)
			{
				delta = icv.run(n);
			}
			else if (RandomNumberGenerator::Type::MM1 == p.type)
			{
				delta = icv.run_mm1(n);
			}

			int C = 0;
			if (icv.get_ideal_half_length() > r * dataMean)
			{
				C = 0;
			}
			else
			{
				C = 1;
			}

			if (C == 0)
			{
				A = n;
			}
			else
			{
				B = n;
			}
			n = (A + B) / 2;
			if ((B - A) <= 10)
			{
				break;
			}
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
