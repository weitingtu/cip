#include <stdio.h>
#include "VAMPIRE.h"
#include "Skart.h"
#include "RandomNumberGenerator.h"

int main()
{
	//VAMPIRE vampire;
	//vampire.run();

	double xsd = 2.1;
	float phi = 0.8f;
	int iseed = 0;

    Skart s(xsd, phi, iseed);
	double alpha = 0.2;
	s.skart_procedure( alpha);

	return 0;
}