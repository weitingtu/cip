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

	double xmean = 100.0;
	double xsd = 2.1;
	float phi = 0.8f;
	int iseed = 1;
	double alpha = 0.2;

	if (run_vampire)
	{
		VAMPIRE vampire;
		if (run_skart)
		{
			vampire.run(VAMPIRE::CIP_TYPE::SKART);
		}
		else if (run_asap3)
		{
			vampire.run(VAMPIRE::CIP_TYPE::ASAP3);
		}
	}
	else if (run_bisection)
	{
		BisectionSearch b;
		if (run_skart)
		{
		    b.run_skart(xmean, xsd, phi, iseed);
		}
		else if (run_asap3)
		{
		    b.run_asap3(xmean, xsd, phi, iseed);
		}
	}
	else if (run_skart)
	{
		Skart s(xmean, xsd, phi, iseed);
		s.skart_procedure(alpha);
	}
	else if (run_asap3)
	{
		xsd = 20;
		ASAP3 a(xmean, xsd, phi, iseed);
		a.procedure(alpha);
	}

	return 0;
}