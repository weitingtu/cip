#pragma once
#include <string>
#include <vector>

class Skart
{
	enum precType {
		absolute,
		relative
	};
public:
	Skart(double xsd, float phi, int iseed);
	~Skart();

	void skart_procedure(double alpha );
	void skart_procedure(
		const std::string& model,
		bool precReq,
		precType precisionType,
		double alpha,
		double hrstar);

	double get_CIlb() const { return _CIlb; }
	double get_CIub() const { return _CIub; }
	double get_data_mean() const { return _data_mean;  }
	const std::vector<double>& get_observation() const { return _observation; }

private:
	double SkewnessFun(const std::vector<double>& myData, int cutOff = 0) const;
	bool vonNuemannTest(const std::vector<double>& myData, double alpha) const;
	double SkewnessAdj(const std::vector<double>& nonspacedbatch, double kPrime, double m, double alpha,
		double sampleMean, double sampleVar, double standardError, int w, double& G1, double& G2) const;

	std::vector<double> runSimulation(const std::string& model, 
		const std::vector<double>&givenData,
		int batchsize, int batchcount, int spacerLength = 0);

private:
	double _xsd;
	float _phi;
	int _iseed;
	std::vector<double> _data;
	std::vector<double> _observation;
	double _CIlb;
	double _CIub; 
	double _data_mean;
};

