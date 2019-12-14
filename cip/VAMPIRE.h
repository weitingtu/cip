#pragma once
class VAMPIRE
{
public:
	enum class CIP_TYPE {
		SKART,
		ASAP3
	};
	VAMPIRE();
	~VAMPIRE();

	double run(CIP_TYPE cip_type);
};

