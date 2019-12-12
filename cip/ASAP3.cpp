#include "ASAP3.h"
#include "RandomNumberGenerator.h"
#include <Algorithm>



ASAP3::ASAP3(double xsd, float phi, int iseed) :
	_xsd(xsd), // 2.1
	_phi(phi), // 0.8
	_iseed(iseed),
	_data()
{
}


ASAP3::~ASAP3()
{
}

void ASAP3::procedure(double alpha)
{
	// step [0]
	int index = 1;
	int m = 16;    // batch size
	int k = 256;   // batch count
	int n = m * k; // sample size
	int k_prime = k - 4;
	alpha = 0.1;
	double nominal_coverage = 1 - alpha;
	double alpha_arp = 0.01;
	double delta = 0.1;
	double omega = 0.18421;
	bool MVTestPassed = false;
	// relative precision
	bool RelPrec = true;
	double r_star = 0.15;
	double H_star = 0.0;

	while (true)
	{
		// step [1]
		std::vector<double> data = generate(n);
		std::vector<double> batch_means(k, 0.0);
		for (int k_i = 0; k_i < k; ++k_i)
		{
			double batch_mean = 0.0;
			for (int m_i = 0; m_i < m; ++m_i)
			{
				batch_mean += data.at(k_i * m + m_i);
			}
			batch_mean /= m;
			batch_means.at(k_i) = batch_mean;
		}

		double truncated_grand_mean = 0.0;
		for (size_t i = 4 * m; i < data.size(); ++i)
		{
			truncated_grand_mean += data.at(i);
		}
		truncated_grand_mean /= (data.size() - 4 * m);

		if (!MVTestPassed)
		{
			// step [2]
			if (true)
			{
				// go to [1]
				continue;
			}
			else
			{
				MVTestPassed = true;
			    // go to [3]
			}
		}

		// step [3]
		// ...

		if (true)
		{
			// go to [1]
			continue;
		}
		else
		{
			// go to [4]
		}

        // step [4]
		
		// step [5]

		if (true)
		{
			break;
		}
		else
		{
			if (true)
			{
				// go to [1]
			}
			else
			{
				// go to [1]
			}
		}
	}
}

std::vector<double> ASAP3::generate( int n )
{
	while((int)_data.size() < n )
	{
		_data.push_back(RandomNumberGenerator::generator(10.5, _xsd, _phi, 10, &_iseed));
	}

	return _data;
}

/*

        Algorithm AS 181.2   Appl.Statist.  (1982) Vol. 31, No. 2

        Calculates the algebraic polynomial of order nored - 1 with
        array of coefficients c.Zero order coefficient is c(1)

*/

double POLY(const std::vector<double>& c, int nord, double x)
{
	double poly_r = c.at(0);

	if (1 == nord)
	{
		return poly_r;
	}

	double p = x * c.at(nord - 1);

	if (2 == nord)
	{
		poly_r = poly_r + p;
		return poly_r;
	}

	int n2 = nord - 2;
	int j = n2 + 1;

	for (int i = 0; i < n2; ++i)
	{
		p = (p + c.at(j - 1)) * x;
		j = j - 1;
	}

	poly_r = poly_r + p;
	return poly_r;
}

/*
         Algorithm AS66 Applied Statistics(1973) vol22 no.3

       Evaluates the tail area of the standardised normal curve
       from x to infinity if upper is .true. or
       from minus infinity to x if upper is .false.

*/

double ALNORM(double x, bool upper)
{
	double zero = 0.0;
	double one  = 1.0;
	double half = 0.5;

	double con = 1.28;
	double z = 0.0;
	double y = 0.0;

	double ltone  = 7.0;
	double utzero = 18.66;
	double p = 0.398942280444;
	double q = 0.39990348504;
	double r = 0.398942280385;
	double a1 = 5.75885480458;
	double a2 = 2.62433121679;
	double a3 = 5.92885724438;
	double b1 = -29.8213557807;
	double b2 = 48.6959930692;
	double c1 = -3.8052e-8;
	double c2 = 3.98064794e-4;
	double c3 = -0.151679116635;
	double c4 = 4.8385912808;
	double c5 = 0.742380924027;
	double c6 = 3.99019417011;
	double d1 = 1.00000615302;
	double d2 = 1.98615381364;
	double d3 = 5.29330324926;
	double d4 = -15.1508972451;
	double d5 = 30.789933034;

	bool up = upper;
	z = x;

	double alnorm_r = zero;
	if (0 == z)
	{
		up = !up;
		z = -z;
	}
	if (z <= ltone || (up && z <= utzero))
	{
		y = half * z * z;
		if (z > con)
		{
			alnorm_r = exp(-y) / (z + c1 + d1 / ( z + c2 + d2 / ( z + c3 + d3 / (z + c4 + d4 / (z + c5 + d5 / ( z + c6))))));
		}
		else
		{
			alnorm_r = half - z * (p - q * y / (y + a1 + b1 / (y + a2 + b2 / (y + a3))));
		}
	}
	else
	{
    	alnorm_r = zero;
	}

	if (!up)
	{
		alnorm_r = one - alnorm_r;
	}
	return alnorm_r;
}

int SIGN(int i, int j)
{
	if (j >= 0)
	{
		return abs(i);
	}
	return -abs(i);
}

double SIGN(double i, double j)
{
	if (j >= 0)
	{
		return abs(i);
	}
	return -abs(i);
}

void SWILK(bool& INIT, const std::vector<double>& X, int N, int N1, int N2, std::vector<double>& A, double& W, double& PW, int& IFAULT)
{
	/*
		ALGORITHM AS R94 APPL.STATIST. (1995) VOL.44, NO.4

		Calculates the Shapiro - Wilk W test and its significance level

	*/

	double SUMM2 = 0.0, SSUMM2 = 0.0, FAC = 0.0, RSN = 0.0, AN = 0.0, AN25 = 0.0, A1 = 0.0, A2 = 0.0, DELTA = 0.0, RANGE = 0.0;
	double SA = 0.0, SX = 0.0, SSX = 0.0, SSA = 0.0, SAX = 0.0, ASA = 0.0, XSX = 0.0, SSASSX = 0.0, W1 = 0.0, Y = 0.0, XX = 0.0, XI;
	double GAMMA = 0.0, M = 0.0, S = 0.0, LD = 0.0, BF = 0.0, Z90F = 0.0, Z95F = 0.0, Z99F = 0.0, ZFM = 0.0, ZSD = 0.0, ZBAR;
	/*
		Auxiliary routines
	*/
	int NCENS = 0, NN2 = 0, I1 = 0;

	std::vector<double>	C1{ 0.0E0, 0.221157E0, -0.147981E0, -0.207119E1, 0.4434685E1, -0.2706056E1 };
	std::vector<double>	C2{ 0.0E0, 0.42981E-1, -0.293762E0, -0.1752461E1, 0.5682633E1, -0.3582633E1 };
	std::vector<double> C3{ 0.5440E0, -0.39978E0, 0.25054E-1, -0.6714E-3 };
	std::vector<double> C4{ 0.13822E1, -0.77857E0, 0.62767E-1, -0.20322E-2 };
	std::vector<double> C5{ -0.15861E1, -0.31082E0, -0.83751E-1, 0.38915E-2 };
	std::vector<double> C6{ -0.4803E0, -0.82676E-1, 0.30302E-2 };
	std::vector<double> C7{ 0.164E0, 0.533E0 };
	std::vector<double> C8{ 0.1736E0, 0.315E0 };
	std::vector<double>	C9{ 0.256E0, -0.635E-2 };
	std::vector<double> G{ -0.2273E1, 0.459E0 };
	double Z90 = 0.12816E1;
	double Z95 = 0.16449E1;
	double Z99 = 0.23263E1;

	double ZM = 0.17509E1;
	double ZSS = 0.56268E0;
	double BF1 = 0.8378E0;
	double XX90 = 0.556E0;
	double XX95 = 0.622E0;
	double ZERO = 0.0E0;
	double ONE = 1.0E0;
	double TWO = 2.0E0;
	double THREE = 3.0E0;
	double SQRTH = 0.70711E0;
	double QTR = 0.25E0;
	double TH = 0.375E0;
	double SMALL = 1E-19;
	double PI6 = 0.1909859E1;
	double STQR = 0.1047198E1;
	bool UPPER = true;

	PW = ONE;
	if (W >= ZERO)
	{
		W = ONE;
	}
	AN = N;
	IFAULT = 3;
	NN2 = N / 2;
	if (N2 < NN2)
	{
		return;
	}
	IFAULT = 1;
	if (N < 3)
	{
		return;
	}
	/*
		If INIT is false, calculates coefficients for the test
	*/
	if (!INIT)
	{
		if (N == 3)
		{
			A.at(0) = SQRTH;
		}
		else
		{
			AN25 = AN + QTR;
			SUMM2 = ZERO;
			for (int I = 1; I <= N2; ++I)
			{
				int ifault = 0;
				A.at(I - 1) = RandomNumberGenerator::ppnd((I - TH) / AN25, &ifault);
				SUMM2 = SUMM2 + A.at(I - 1) * A.at(I - 1);
			}
			SUMM2 = SUMM2 * TWO;
			SSUMM2 = sqrt(SUMM2);
			RSN = ONE / sqrt(AN);
			A1 = POLY(C1, 6, RSN) - A.at(0) / SSUMM2;
			/*
				Normalize coefficients
			*/
			if (N > 5)
			{
				I1 = 3;
				A2 = -A.at(2 - 1) / SSUMM2 + POLY(C2, 6, RSN);
				FAC = sqrt((SUMM2 - TWO * pow(A.at(0), 2) - TWO * pow(A.at(1), 2)) / (ONE - TWO * pow(A1, 2) - TWO * pow(A2, 2)));
				A.at(1 - 1) = A1;
				A.at(2 - 1) = A2;
			}
			else
			{
				I1 = 2;
				FAC = sqrt((SUMM2 - TWO * pow(A.at(0), 2)) / (ONE - TWO * pow(A1, 2)));
				A.at(0) = A1;
			}
			for (int i = I1 - 1; i < NN2; ++i)
			{
				A.at(i) = -A.at(i) / FAC;
			}
		}
		INIT = true;
	}
	if (N1 < 3)
	{
		return;
	}
	NCENS = N - N1;
	IFAULT = 4;
	if (NCENS < 0 || (NCENS > 0 && N < 20))
	{
		return;
	}
	IFAULT = 5;
	DELTA = (double)(NCENS) / AN;
	if (DELTA > 0.8)
	{
		return;
	}
	/*
		If W input as negative, calculate significance level of - W
	*/
	if (W < ZERO)
	{
		W1 = ONE + W;
		IFAULT = 0;
		// GOTO 70
	}
	else
	{
		/*
			Check for zero range
		*/
		IFAULT = 6;
		RANGE = X.at(N1 - 1) - X.at(1 - 1);
		if (RANGE < SMALL)
		{
			return;
		}
		/*
			Check for correct sort order on range - scaled X
		*/
		IFAULT = 7;
		XX = X.at(1 - 1) / RANGE;
		SX = XX;
		SA = -A.at(1 - 1);
		int J = N - 1;
		for (int I = 2; I <= N1; ++I)
		{
			XI = X.at(I - 1) / RANGE;
			if ((XX - XI) > SMALL)
			{
				//PRINT *, ' ANYTHING'
			}
			SX = SX + XI;
			if (I != J)
			{
				SA = SA + SIGN(1, I - J) * A.at(std::min(I, J) - 1);
			}
			XX = XI;
			J = J - 1;
		}
		IFAULT = 0;
		if (N > 5000)
		{
			IFAULT = 2;
		}
		/*
			Calculate W statistic as squared correlation
			between data and coefficients
		*/
		SA = SA / N1;
		SX = SX / N1;
		SSA = ZERO;
		SSX = ZERO;
		SAX = ZERO;
		J = N;
		for (int I = 1; I <= N1; ++I)
		{
			if (I != J)
			{
				ASA = SIGN(1, I - J) * A.at(std::min(I, J) - 1) - SA;
			}
			else
			{
				ASA = -SA;
			}
			XSX = X.at(I - 1) / RANGE - SX;
			SSA = SSA + ASA * ASA;
			SSX = SSX + XSX * XSX;
			SAX = SAX + ASA * XSX;
			J = J - 1;
		}
		/*
			W1 equals(1 - W) claculated to avoid excessive rounding error
			for W very near 1 (a potential problem in very large samples)
		*/
		SSASSX = sqrt(SSA * SSX);
		W1 = (SSASSX - SAX) * (SSASSX + SAX) / (SSA * SSX);
	}
	// 70
	W = ONE - W1;
	/*
		Calculate significance level for W(exact for N = 3)
	*/

	if (N == 3)
	{
		PW = PI6 * (asin(sqrt(W)) - STQR);
		return;
	}
	Y = log(W1);
	XX = log(AN);
	M = ZERO;
	S = ONE;
	if (N <= 11)
	{
		GAMMA = POLY(G, 2, AN);
		if (Y >= GAMMA)
		{
			PW = SMALL;
			return;
		}
		Y = -log(GAMMA - Y);
		M = POLY(C3, 4, AN);
		S = exp(POLY(C4, 4, AN));
	}
	else
	{
		M = POLY(C5, 4, XX);
		S = exp(POLY(C6, 3, XX));
	}

	if (NCENS > 0)
	{
		/*
			Censoring by proportion NCENS / N.Calculate mean and sd
			of normal equivalent deviate of W.
		*/
		LD = -log(DELTA);
		BF = ONE + XX * BF1;
		Z90F = Z90 + BF * pow(POLY(C7, 2, pow(XX90, XX)), LD);
		Z95F = Z95 + BF * pow(POLY(C8, 2, pow(XX95, XX)), LD);
		Z99F = Z99 + BF * pow(POLY(C9, 2, XX), LD);
	    /*
	        Regress Z90F, ..., Z99F on normal deviates Z90, ..., Z99 to get
		    pseudo - mean and pseudo - sd of z as the slope and intercept
	    */
		ZFM = (Z90F + Z95F + Z99F) / THREE;
		ZSD = (Z90*(Z90F - ZFM) + Z95 * (Z95F - ZFM) + Z99 * (Z99F - ZFM)) / ZSS;
		ZBAR = ZFM - ZSD * ZM;
		M = M + ZBAR * S;
		S = S * ZSD;
	}
	PW = ALNORM((Y - M) / S, UPPER);
	return;
}