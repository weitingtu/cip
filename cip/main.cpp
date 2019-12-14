#include <stdio.h>
#include "VAMPIRE.h"
#include "Skart.h"
#include "ASAP3.h"
#include "RandomNumberGenerator.h"

int main()
{
	//VAMPIRE vampire;
	//vampire.run(VAMPIRE::CIP_TYPE::SKART);
	//vampire.run(VAMPIRE::CIP_TYPE::ASAP3);

	double xsd = 2.1;
	float phi = 0.8f;
	int iseed = 1;

    //Skart s(xsd, phi, iseed);
	double alpha = 0.2;
	//s.skart_procedure( alpha);

	xsd = 20;
	ASAP3 a(xsd, phi, iseed);
	a.procedure(alpha);

	//std::vector<double> data;
	//double x = 10.0;
	//xsd = 20;
	//phi = 0.8f;
	//while(data.size() < 10)
	//{
	//	x = RandomNumberGenerator::generator(10.5, xsd, phi, x, &iseed);
	//	data.push_back(x);
	//}

	//for (size_t i = 0; i < data.size(); ++i)
	//{
	//	printf("data %f\n", data.at(i));
	//}

	//double a1 = 0.0;
	//double b1 = 0.0;
	//for (size_t k = 1; k < data.size(); ++k)
	//{
	//	a1 += data.at(k) * data.at(k - 1);
	//	b1 += data.at(k - 1) * data.at(k - 1);
	//}
	//double phi_hat = a1 / b1;
	//printf("phi hat = %f\n", phi_hat);

	return 0;
}