#include <stdio.h>
#include "VAMPIRE.h"
#include "Skart.h"
#include "ASAP3.h"
#include "RandomNumberGenerator.h"
#include "BisectionSearch.h"

int main()
{
	bool run_vampire = false;
	bool run_bisection = true;
	bool run_skart = true;
	bool run_asap3 = false;

	RandomNumberGenerator::Parameter parameter;
	parameter.type = RandomNumberGenerator::Type::AR1;
	parameter.ar1.xmean = 100.0;
	parameter.ar1.xsd = 2.1;
	parameter.ar1.phi = 0.8f;
	parameter.ar1.x   = 100;
	parameter.ar1.iseed = 1;
	parameter.mm1.arate = 0.0;
	parameter.mm1.srate = 0.0;
	parameter.mm1.waitq = -1.0;
	parameter.mm1.iseed = 1;

	double alpha = 0.2;

	bool precReq = false;
	double hrstar = 0.15;
	bool RelPrec = false;
	double r_star = 0.15;

	double tolerance = 0.00015;

	for (int i = 1; i <= 400; ++i)
	{
	    parameter.ar1.iseed = i;
	    parameter.mm1.iseed = i;

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
		}
		else if (run_bisection)
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