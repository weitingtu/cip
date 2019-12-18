#pragma once
#include <string>
#include <vector>
#include "RandomNumberGenerator.h"

class Skart
{
	enum precType {
		absolute,
		relative
	};
public:
	Skart(const RandomNumberGenerator::Parameter& p);
	~Skart();

	void skart_procedure(bool precReq, double alpha, double hrstar );
	void skart_procedure(
		const std::string& model,
		bool precReq,
		precType precisionType,
		double alpha,
		double hrstar);

	double get_CIlb() const { return _CIlb; }
	double get_CIub() const { return _CIub; }
	double get_data_mean() const { return _data_mean;  }
	double get_xmean() const { return _xmean;  }
	const std::vector<double>& get_data() const { return _data; }
	const std::vector<double>& get_observation() const { return _observation; }

	bool vonNuemannTest(const std::vector<double>& myData, double alpha) const;
private:
	double SkewnessFun(const std::vector<double>& myData, int cutOff = 0) const;
	double SkewnessAdj(const std::vector<double>& nonspacedbatch, double kPrime, double m, double alpha,
		double sampleMean, double sampleVar, double standardError, int w, double& G1, double& G2) const;

	std::vector<double> runSimulation(const std::string& model, 
		const std::vector<double>&givenData,
		int batchsize, int batchcount, int spacerLength = 0);

private:
	RandomNumberGenerator::Parameter _parameter;
	std::vector<double> _data;
	std::vector<double> _observation;
	double _CIlb;
	double _CIub; 
	double _data_mean;
	double _xmean;
};

