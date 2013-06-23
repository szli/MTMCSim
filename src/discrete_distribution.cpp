#include "discrete_distribution.h"

namespace MTMCSim{
discrete_distribution::discrete_distribution(unsigned int n_, const std::vector<double> &prob_)
: n(n_), prob(prob_), mean(0), stdev(0),hasValidTable(false)
{
	genTab();

	
	// cleared to zero in the init list
	for (unsigned int i = 0; i<n; i++)
		mean += prob[i] * i;

	// cleared to zero in the init list
	for (unsigned int i = 0; i<n ; i++)
		stdev += (i - mean) * (i - mean) * prob[i]; // variance

	stdev = sqrt(stdev);
}

void discrete_distribution::setDistribution(unsigned int n_, const std::vector<double> &prob_)
{
	n = n_;
	prob = prob_;
	mean = 0;
	stdev = 0;

	genTab();

	for (unsigned int i = 0; i<n; i++)
		mean += prob[i] * i;


	for (unsigned int i = 0; i<n ; i++)
		stdev += (i - mean) * (i - mean) * prob[i]; // variance

	stdev = sqrt(stdev);
}

void discrete_distribution::setDistributionWOGenTab(unsigned int n_, const std::vector<double> &prob_)
{
	n = n_;
	prob = prob_;
	mean = 0;
	stdev = 0;

	for (unsigned int i = 0; i<n; i++)
		mean += prob[i] * i;


	for (unsigned int i = 0; i<n ; i++)
		stdev += (i - mean) * (i - mean) * prob[i]; // variance

	stdev = sqrt(stdev);
}

discrete_distribution::discrete_distribution(unsigned int n_, const std::vector<double> &prob_, const std::vector<double> &val_)
: n(n_), prob(prob_), val(val_), mean(0), stdev(0),hasValidTable(false)
{
	genTab();

	for (unsigned int i = 0; i<n; i++)
		mean += prob[i] * val[i];

	// cleared to zero in the init list
	for (unsigned int i = 0; i<n ; i++)
		stdev += (val[i] - mean) * (val[i] - mean) * prob[i]; // variance

	stdev = sqrt(stdev);

}

void discrete_distribution::genTab()
{
	if ((n == 0) || (prob.empty()))
		throw runtime_error("Probability distribution not set properly!");

	p.resize(n);
	ialt.resize(n);

	double rn = 1.0/n; 
	double px = 0;

	for (unsigned int i = 0; i<n; i++)
	{
		p[i] = 0.0;
		ialt[i] = -1;
	}

	for (unsigned int i = 0; i<n; i++)
	{
		double pmin = 1.0;
		double pmax = -1.0;
		int jmin = -1;
		int jmax = -1;

		for (unsigned int j = 0; j<n; j++)
		{
			if (ialt[j]<0)
			{
				px = prob[j] - p[j];
				if (pmax <= px)
				{
					pmax = px;
					jmax = j;
				}

				if (pmin >= px)
				{
					pmin = px;
					jmin = j;
				}
			}
		}

		if ((jmax < 0) || (jmin<0)) 
			throw runtime_error("RNG method fails!");

		ialt[jmin] = jmax;
		px = rn - pmin;
		p[jmax] += px;
		px *= n;
		p[jmin] = px;
	} // end for i 
	hasValidTable = true;
}


void discrete_distribution::loadTab(string& tabStr)
{
	p.resize(n);
	ialt.resize(n);

	stringstream sstab(tabStr);
	for (unsigned int i = 0; i<n ;i++)
		sstab >> ialt[i];

	for (unsigned int i = 0; i<n; i++)
		sstab >> p[i];
	hasValidTable = true;
}

void discrete_distribution::tabToStr(string& tabStr)
{
	if (!hasValidTable)
		throw runtime_error("Has not set the table yet!!");

	stringstream sstab;

	for (unsigned int i = 0; i<n; i++)
		sstab << ialt[i] << " ";

	for (unsigned int i = 0; i<n; i++)
		sstab << p[i] << " ";

	tabStr = sstab.str();
}


}