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
	int iseed = 1;

    Skart s(xsd, phi, iseed);
	double alpha = 0.2;
	s.skart_procedure( alpha);

	//double x = 10.0;
	//std::vector<double> data;
	//for (size_t i = 0; i < 10; ++i)
	//{
	//	x = RandomNumberGenerator::generator(10.5, xsd, phi, x, &iseed);
	//	printf("%f\n", x);
	//	data.push_back(x);
	//}

	//bool r = s.vonNuemannTest(data, alpha);
	//printf("r = %s\n", r ? "true" : "false");

	return 0;
}