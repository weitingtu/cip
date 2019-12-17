#include "Skart.h"
#include "RandomNumberGenerator.h"
#include <algorithm>
#include <math.h>

Skart::Skart(double xmean, double xsd, float phi, int iseed) :
	_x(100.0),
	_xmean(xmean),
	_xsd(xsd),
	_phi(phi),
	_iseed(iseed),
	_data(),
	_observation(),
	_CIlb(0.0),
	_CIub(0.0),
	_data_mean()
{
}

Skart::~Skart()
{
}

void Skart::skart_procedure(double alpha )
{
	//skart_procedure("MM1", true, precType::relative, alpha, 0.15);
	skart_procedure("MM1", false, precType::relative, alpha, 0.15);
}

void Skart::skart_procedure(
	const std::string& model,
	bool precReq,
	precType precisionType,
	double alpha,
	double hrstar
)
{
	int i = 0;
	int j = 0; 
	int k = 1280;   // initial batch count
	int kPrime = 0; // current batch count
	int m = 0;      // batch size
	double alphaRan = 0.2; // randomness test size
	int dStar = 0;  // maximum number of batches allowed in spacer
	int d = 0;      // current number of batches in a spacer
	int w = 0;      // spacer size
	int b = 0;      // number of times the batch count has been deflated
	std::vector<double> nonspacedbatch(k, 0.0);
	std::vector<double> spacedBatch(k, 0.0);
	double sampleMean = 0.0;
	double sampleVar = 0.0;
	double standardError = 0.0;
	double Skewness = 0.0;
	double halfLength = 0.0;
	double halfLengthStar = 0.0;
	double A = 0.0;
	double G1 = 0.0;
	double G2 = 0.0;
	
	printf("model = %s\n", model.c_str());
	printf("precision requirement = %s\n", precReq ? "true" : "false");
	if (precisionType == absolute)
	{
		printf("precision Type = absolute\n");
	}
	else
	{
		printf("precision Type = relative\n");
	}
	printf("coverage probability = %f%%\n", (1 - alpha) * 100);
	printf("hrstar = %f%%\n", hrstar * 100);
	printf("\n");

	//generates initial sample of size 1,280
	std::vector<double> initialData(1, 0.0);
	std::vector<double> data = runSimulation(model, initialData, 1, k);
	printf("initial sample size = %d * %d = %d %d\n", 1, k, 1 * k, __LINE__);

	//computes the sample skewness of the last 1024 observations in the sample: 1,024=1,280-256
	Skewness = SkewnessFun(data, 256);

	// sets the initial value of the batch size based on the computed sample skewness 
	if (Skewness <= 4)
	{
		m = 1;
		printf("initial sample size = %d * %d = %d %d\n", m, k, m * k, __LINE__);
	}
	else
	{
		m = 16;
		// collects data: 20,480=1,280*16
		data = runSimulation(model, data, m, k);
		printf("initial sample size = %d * %d = %d %d\n", m, k, m * k, __LINE__);
	}
	printf("\n");

	// computes the batch means
	for (int j = 0; j < k; ++j)
	{
		for (int i = 0; i < m; ++i)
		{
			nonspacedbatch[j] += data[j * m + i] / m;
		}
	}

	printf("randomness test:\n");
	bool randomnessStatus = false;
	while (!randomnessStatus)
	{
		// sets the maximum number of batches allowed in a spacer
		Skewness = SkewnessFun(nonspacedbatch, (int)(k * 0.2));
		dStar = 10;
		if (abs(Skewness) > 0.5)
		{
			dStar = 3;
		}

		randomnessStatus = vonNuemannTest(nonspacedbatch, alphaRan);
		if (!randomnessStatus)
		{
			while (!randomnessStatus)
			{
				d += 1;
				kPrime = 0;
				spacedBatch = std::vector<double>((int)floor(k / (d + 1)), 0.0);
				for (j = d; j < k; j += (d + 1))
				{
					spacedBatch[kPrime] = nonspacedbatch[j];
					kPrime += 1;
				}
				randomnessStatus = vonNuemannTest(spacedBatch, alphaRan);
				if (d = dStar)
				{
					break;
				}
			}
			if(!randomnessStatus && d == dStar)
			{
				m = (int)ceil(sqrt(2) * m);
				k = (int)ceil(0.9 * k);
				data = runSimulation(model, data, m, k);
				printf("new sample size = %d * %d = %d %d\n", m, k, m * k, __LINE__);
				d = 0;
				b += 1;
				for (j = 0; j < k; ++j)
				{
					nonspacedbatch.at(j) = 0;
					for (i = 0; i < m; ++i)
					{
						nonspacedbatch.at(j) += data.at(j * m + i) / m;
					}
				}
			}
		}
		else
		{
			kPrime = k;
		}
	}

	// skips the warm-up period
	w = d * m;
	data = runSimulation(model, data, m, kPrime, w);
	printf("new sample size = %d * %d = %d %d\n", m, kPrime, m * kPrime, __LINE__);

    kPrime = (int)ceil(kPrime * pow(1.0 / 0.9, b));
	m = std::max(m, (int)floor(data.size() / kPrime));
	nonspacedbatch = std::vector<double>(kPrime, 0.0);
	if (m * kPrime > (int) data.size())
	{
		data = runSimulation(model, data, m, kPrime);
		printf("new sample size = %d * %d = %d %d\n", m, kPrime, m * kPrime, __LINE__);
	}

	 // computes the current set of truncated, nonspaced batch means
	for (j = 0; j < kPrime; ++j)
	{
		for (i = 0; i < m; ++i)
		{
			nonspacedbatch.at(j) += data.at(j * m + i) / m;
		}
	}

    // computes the sample mean and variance
	for (j = 0; j < kPrime; ++j)
	{
		sampleMean += nonspacedbatch.at(j) / kPrime;
	}
	for (j = 0; j < kPrime; ++j)
	{
		sampleVar += pow(nonspacedbatch.at(j) - sampleMean, 2) / (kPrime - 1);
	}

	// computes the correlation adjustment and the CI's standard error
	double sigma = 0.0;
	double phiHat = 0.0;
	for (j = 0; j < kPrime - 1; ++j)
	{
		sigma += (nonspacedbatch.at(j) - sampleMean) * (nonspacedbatch.at(j + 1) - sampleMean);
	}
	phiHat = sigma / (sampleVar * (kPrime - 1));
	A = (1 + phiHat) / (1 - phiHat);
	standardError = sqrt(A * sampleVar / kPrime);
	if (precReq && precisionType == precType::relative)
	{
		halfLengthStar = hrstar * sampleMean;
	}
	else
	{
		halfLengthStar = hrstar;
	}
	halfLength = SkewnessAdj(nonspacedbatch, kPrime, m, alpha, sampleMean, sampleVar, standardError, w, G1, G2);

	double kPrimeNew = 0.0;
	double maxVal = 0.0;
	printf("\n");
	if (precReq)
	{
		printf("precision requirement:\n");
	}
	while (halfLength > halfLengthStar && precReq)
	{
		kPrimeNew = ceil((pow(halfLength / halfLengthStar, 2)) * kPrime);
		if (kPrime == 1024)
		{
			maxVal = std::max(1.02, (pow(halfLength / halfLengthStar, 2)));
			m = (int)ceil(std::min(maxVal, (double)2) * m);
		}
		kPrime = std::min((int)kPrimeNew, 1024);
		nonspacedbatch = std::vector<double>(kPrime - 1, 0.0);
		data = runSimulation(model, data, m, kPrime);
		printf("new sample size = %d * %d = %d %d\n", m, kPrime, m * kPrime, __LINE__);
		for (j = 0; j < kPrime; ++j)
		{
			nonspacedbatch.at(j) = 0;
			for (i = 0; i < m; ++i)
			{
				nonspacedbatch.at(j) += data.at(j * m + i) / m;
			}
		}
		sampleVar = 0;
		sampleMean = 0;
		// computes the sample mean and variance
		for (j = 0; j < kPrime; ++j)
		{
			sampleMean += nonspacedbatch.at(j) / kPrime;
		}
		for (j = 0; j < kPrime; ++j)
		{
			sampleVar += pow(nonspacedbatch.at(j) - sampleMean, 2) / (kPrime - 1);
		}

		// computes the correlation adjustment and the CI's standard error
		sigma = 0;
		for (j = 0; j < kPrime - 1; ++j)
		{
			sigma += (nonspacedbatch.at(j) - sampleMean) * (nonspacedbatch.at(j + 1) - sampleMean);
		}
		phiHat = sigma / (sampleVar * (kPrime - 1));
		A = (1 + phiHat) / (1 - phiHat);
		standardError = sqrt(A * sampleVar / kPrime);
		if (precReq && precisionType == precType::relative)
		{
			halfLengthStar = hrstar * sampleMean;
		}
		else
		{
			halfLengthStar = hrstar;
	    }
		halfLength = SkewnessAdj(nonspacedbatch, kPrime, m, alpha, sampleMean, sampleVar, standardError, w, G1, G2);
	}

	double CIlb = 0.0;
	double CIub = 0.0; 
	printf("sample mean %f half length %f\n", sampleMean, halfLength);
	if (precReq)
	{
		CIlb = sampleMean - halfLength;
		CIub = sampleMean + halfLength;
	}
	else // 'no precision requirement
	{
		CIlb = sampleMean - G1 * standardError;
		CIub = sampleMean - G2 * standardError;
	}

	printf("\n");
	printf("warm-up period = %d\n", w);
	printf("final sample size = %d\n", kPrime * m);
	printf("\n");
	printf("sample mean = %.2f\n", sampleMean);
	printf("sample variance = %.2f\n", sampleVar);
	printf("\n");
	printf("cil = %.2f\n", CIlb);
	printf("ciu = %.2f\n", CIub);

	_CIlb = CIlb;
	_CIub = CIub;
	// computes the data means
	_data_mean = 0.0;
	for (size_t i = 0; i < _data.size(); ++i)
	{
		_data_mean += _data[i];
	}
	_data_mean /= _data.size();
	_observation = data;
}

double Skart::SkewnessFun(const std::vector<double>& myData, int cutOff) const
{
	double mean = 0.0;
	double mom2 = 0.0;
	double mom3 = 0.0;

	int dataLength = (int)myData.size();
	for (int j = cutOff; j < dataLength; ++j)
	{
		mean += myData[j];
	}
	int size = dataLength - cutOff;
	mean /= size;

	for (int j = cutOff; j < dataLength; ++j)
	{
		mom2 += pow(myData[j] - mean, 2);
		mom3 += pow(myData[j] - mean, 3);
	}

	double skewness = (pow(size, 2.5) * mom3) / (pow(mom2, 1.5) * (size - 1)) * (size - 2);
	return skewness;
}

bool Skart::vonNuemannTest(const std::vector<double>& myData, double alpha) const
{
	int j = 0;
	int dataSize = myData.size();

	// Compute the von Neumann test statistic
	double sampleMean = 0.0;
	double sampleVar  = 0.0;

	for (int j = 0; j < dataSize; ++j)
	{
		sampleMean += myData[j];
	}
	sampleMean /= dataSize;

	for (int j = 0; j < dataSize; ++j)
	{
		sampleVar += pow(myData[j] - sampleMean, 2);
	}

	// Compute the mean square successive difference
	double succDiff = 0;
	for (int j = 0; j < dataSize - 1; ++j)
	{
		succDiff += pow(myData[j] - myData[j + 1], 2);
	}
	double ranTestStat = 1 - (succDiff / (2 * sampleVar));

	// take alpha=0.2 for the 2-sided test
	int ifault = 0;
	double p = RandomNumberGenerator::ppnd(1 - alpha / 2, &ifault);
	if ( abs( ranTestStat ) > p * sqrt( ( dataSize - 2 ) / ( pow( dataSize, 2 ) - 1 ) ) )
	{
		return false;
	}

	return true;
}

double Skart::SkewnessAdj(const std::vector<double>& nonspacedbatch, double kPrime, double m, double alpha,
	double sampleMean, double sampleVar, double standardError, int w, double& G1, double& G2) const
{
	int j = 0;
	std::vector<double> spacedBatch((int)kPrime, 0.0);

	// computes spaced batch means to calculate the skewness adjustment
	int dPrime = (int)(ceil(w / m)); // number of batches in a spacer
	int kDoublePrime = 0;
	spacedBatch = std::vector<double>((int)floor(kPrime / (dPrime + 1)), 0.0);
	for (int j = dPrime; j < kPrime; j += (dPrime + 1))
	{
		spacedBatch[kDoublePrime] = nonspacedbatch[j];
		++kDoublePrime;
	}

	double spacedSampleMean = 0.0;
	double spacedSampleVar  = 0.0;
	double spacedMom3 = 0.0;
	for (int j = 0; j < kDoublePrime; ++j)
	{
		spacedSampleMean += spacedBatch[j] / kDoublePrime;
	}
	for (int j = 0; j < kDoublePrime; ++j)
	{
		spacedSampleVar += pow(spacedBatch[j] - spacedSampleMean, 2) / (kDoublePrime - 1);
		spacedMom3 += (kDoublePrime * pow(spacedBatch[j] - spacedSampleMean, 3)) / ((kDoublePrime - 1) * (kDoublePrime - 2));
	}

	double beta = (spacedMom3 / pow(spacedSampleVar, 1.5)) / (6 * sqrt(kDoublePrime));
	int ifault = 0;
	int i = 0;
	double t1 = RandomNumberGenerator::tinv(1 - alpha / 2, kDoublePrime - 1, 0, &ifault, &i);

	double t2 = -t1;
	if ((1 + 6 * beta * (t1 - beta)) < 0)
	{
		G1 = pow(2 * beta, -1) * ( -(pow(abs(1 + 6 * beta * (t1 - beta)), 1.0 / 3)) - 1);
	}
	else
	{
		G1 = pow(2 * beta, -1) * ( pow(1 + 6 * beta * (t1 - beta), 1.0 / 3) - 1);
	}
	if ((1 + 6 * beta * (t2 - beta)) < 0)
	{
		G2 = pow(2 * beta, -1) * ( -(pow(abs(1 + 6 * beta * (t2 - beta)), 1.0 / 3)) - 1);
	}
	else
	{
		G2 = pow(2 * beta, -1) * ( pow(1 + 6 * beta * (t2 - beta), 1.0 / 3) - 1);
	}
	return std::max(standardError*abs(G1), standardError * abs(G2));
}

std::vector<double> Skart::runSimulation( const std::string& model,
	const std::vector<double>&givenData,
	int batchsize,   // 資料分批時每一批的個數
	int batchcount,  // 有幾批
	int spacerLength // 批與批之間的間隔
)
{
	while( (int) _data.size() < batchsize * batchcount )
	{
		_x = RandomNumberGenerator::AR1_generator(_xmean, _xsd, _phi, _x, &_iseed);
		_data.push_back(_x);
	}

	return std::vector<double>(_data.begin(), _data.begin() + batchsize * batchcount);
}
