#include <stdio.h>
#include "VAMPIRE.h"
#include "Skart.h"
#include "ASAP3.h"
#include "RandomNumberGenerator.h"
#include "BisectionSearch.h"

int main()
{
	bool run_vampire = true;
	bool run_bisection = false;
	bool run_skart = true;
	bool run_asap3 = false;

	RandomNumberGenerator::Parameter parameter;
	parameter.type = RandomNumberGenerator::Type::AR1;
	parameter.ar1.xmean = 100.0;
	parameter.ar1.xsd = 10.0;
	parameter.ar1.phi = 0.001f;
	parameter.ar1.x   = 100;
	parameter.ar1.iseed = 1;
	parameter.mm1.arate = 0.8;
	parameter.mm1.srate = 1.0;
	parameter.mm1.waitq = -1.0;
	parameter.mm1.iseed = 1;

	double alpha = 0.05;

	bool precReq = false;
	double hrstar = 0.0375;
	bool RelPrec = true;
	double r_star = 0.0;

	double tolerance = 0.015;

    if (run_vampire)
	{
		VAMPIRE vampire;
		if (run_skart)
		{
			vampire.run(VAMPIRE::CIP_TYPE::SKART,
				parameter,
				alpha,
				precReq,
				hrstar,
				RelPrec,
				r_star,
				tolerance);
		}
		else if (run_asap3)
		{
			vampire.run(VAMPIRE::CIP_TYPE::ASAP3,
				parameter,
				alpha,
				precReq,
				hrstar,
				RelPrec,
				r_star,
				tolerance);
		}
		return 0;
    }

	for (int i = 1; i <= 400; ++i)
	{
	    parameter.ar1.iseed = i;
	    parameter.mm1.iseed = i;

		if (run_bisection)
		{
			BisectionSearch b(tolerance);
			if (run_skart)
			{
				b.run_skart(parameter, precReq, alpha, hrstar);
			}
			else if (run_asap3)
			{
				b.run_asap3(parameter, RelPrec, alpha, r_star);
			}
		}
		else if (run_skart)
		{
			Skart s(parameter);
			s.skart_procedure(precReq, alpha, hrstar);
		}
		else if (run_asap3)
		{
			parameter.ar1.xsd = 20.0;
			ASAP3 a(parameter);
			a.procedure(RelPrec, alpha, r_star);
		}
	}

	return 0;
}