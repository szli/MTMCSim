#ifndef _DISCRETE_DISTRIBUTION_H
#define _DISCRETE_DISTRIBUTION_H
#include "StdHeaders.h"


namespace MTMCSim{

/** Discrete Distribution, generates random number according to given pmf from uniform RVs.
*
*	Follows distribution style of Boost::random library. 
*	copy safe.
*/
class discrete_distribution {
public:
	discrete_distribution()
	{}

	/** Constructor.
	*
	*	Alphabet {0,1,..., n-1}
	*/
	discrete_distribution(unsigned int n_, const vector<double>& prob_);

	/** Constructor.
	*	@param n_ alphabet size
	*	@param prob_ pmf
	*	@param val_ discrete alphabet values, the size should equal to n_
	*/
	discrete_distribution(unsigned int n_, const vector<double>& prob_, const vector<double>& val_);

	void setDistribution(unsigned int n_, const vector<double>& prob_);

	void setDistributionWOGenTab(unsigned int n_, const vector<double>& prob_);

	/** Random number generation operator.
	*
	*	@param engine a uniform pseudorandom number generator that satisfies boost::random concept
	*/
	template<typename Engine> double operator()(Engine & engine) const;

	double getMean() const
	{ return mean;}

	/** Standard deviation.
	*/
	double getStdev() const
	{ return stdev;}

	/** Number of Uniform RV need to generate one RV.
	*/
	// To generate one RV with desired distribution, 
	// two uniformly distributed RV is needed
	int consumeUniRV() const
	{ return 2; } 
	
	void loadTab(string& tabStr);

	void tabToStr(string& tabStr);

private:
	void genTab();
	unsigned int n;
	vector<double> prob;
	vector<double> val;
	double mean, stdev;
	vector<double> p; // for internal algorithm use
	vector<int> ialt; // for internal algorithm use

	bool hasValidTable;
};


template<typename Engine> double discrete_distribution::operator()(Engine & uniEng) const
{
	if (!hasValidTable)
		throw runtime_error("Has not set the table yet!!");
	int i = (int)(n*uniEng());
	if (uniEng() < p[i])
		i = ialt[i];
	return val.empty() ? i : val[i];
}


}

#endif