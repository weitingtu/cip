#include <stdio.h>
#include "VAMPIRE.h"
#include "Skart.h"
#include "ASAP3.h"
#include "RandomNumberGenerator.h"
#include "BisectionSearch.h"

int main()
{
	//VAMPIRE vampire;
	//vampire.run(VAMPIRE::CIP_TYPE::SKART);
	//vampire.run(VAMPIRE::CIP_TYPE::ASAP3);

	double xsd = 2.1;
	float phi = 0.8f;
	int iseed = 1;

	//BisectionSearch b;
	//b.run_skart(xsd, phi, iseed);
	//b.run_asap3(xsd, phi, iseed);

    //Skart s(xsd, phi, iseed);
	double alpha = 0.2;
	//s.skart_procedure( alpha);

	xsd = 20;
	ASAP3 a(xsd, phi, iseed);
	a.procedure(alpha);

	return 0;
}