#ifndef _INCLUDED_NBTOBINPDFCONV_H
#define _INCLUDED_NBTOBINPDFCONV_H

#include "StdHeaders.h"

namespace MTMCSim{

/** Converts a pmf on alphabet with size nbSize to log_2(nbSize) binary sources with a joint pdf. Also computes conditional entropies
*	and conditional probabilities.
*/
class nbToBinPDFConv{

public:
	nbToBinPDFConv(const vector<double>& nbPdf_, int nbSize_);

	nbToBinPDFConv(double* nbPdf_, int nbSize_);

	nbToBinPDFConv()
	{ }

	void setPdf(const std::vector<double> &nbPdf_, int nbSize_);
	
	// not implemented
	//vector<double>& getCondZeroProb(int nKnown); // given nKnown bits, the zero prob of the nKnown+1 bit, returns a vector that is indexed by the given bits (dec experssion)

	

	/** P(b_(len) = 0 | b0,b1,..., b_{len-1} ), b0, b1, ..., b_{len-1} is given by int* knownBits.
	* For generating init LLR.
	*/
	double getZeroPGiven(int* knownBits, int len) const;

	/** H(b0|Y=y), H(b1|b0, Y=y), ... \  H(bi-1|bi-2...b0, Y=y) */
	double getOneCondEn(int index) const
	{
		return condEntro[index];
	}


private:
	int nbSize;

	int nBits;

	vector<unsigned int> QtoBinMap;///< Symbol in Q-alphabet --> binary bits (expressed in terms of unsiged integer)

	vector<vector<double> > condZeroP;///< p(b1=0|b0), p(b2=0|b0,b1), ...., (i,j)th item: P(b_i+1 = 0 | (b0,b1,..., bi) = j(decimal form) ) 

	vector<vector<double> > incrJoPdf; ///< pdf of b0, (b0,b1), (b0,b1,b2), ...(b0,b1,...b_i-1)

	void setQtoBinMap();
	

	void convert(); 

	/** a pdf on discrete alphabet {0, 1, ..., nbSize - 1}.
	* In this class, consider one odf only, in application  
	* use p(x|y) to initialize each object when y is fixed to some value.
	*/
	vector<double> nbPdf;

	void init(); ///< called in constructor, will call convert().

	vector<double> entro; ///<  H(b0|Y=y), H(b0,b1|Y=y), H(b0,b1,b2|Y=y), ..., H(b0, ..., bi-1|Y=y), actually this class does not care about y, just view P(x|y) as a pdf of x for some fixed y

	vector<double> condEntro; ///< H(b0|Y=y), H(b1|b0, Y=y), ... H(bi-1|bi-2...b0, Y=y)

	void compEntro();


};




}

#endif