#pragma once
#include "RandomNumbergenerator.h"

class VAMPIRE
{
public:
	enum class CIP_TYPE {
		SKART,
		ASAP3
	};
	VAMPIRE();
	~VAMPIRE();

	double run(CIP_TYPE cip_type,
		const RandomNumberGenerator::Parameter& parameter,
		double alpha,
		bool precReq,
		double hrstar,
		bool RelPrec,
		double r_star,
		double tolerance);
};

