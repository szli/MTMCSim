#include "nbToBinPDFConv.h"
//#include "BasicOp.h"

namespace MTMCSim{
nbToBinPDFConv::nbToBinPDFConv(const std::vector<double> &nbPdf_, int nbSize_)
{
	nbSize = nbSize_;

	nbPdf = nbPdf_;

	init();

	

}

void nbToBinPDFConv::setPdf(const std::vector<double> &nbPdf_, int nbSize_)
{
	nbSize = nbSize_;

	nbPdf = nbPdf_;

	init();
}





void nbToBinPDFConv::init()
{
	nBits = static_cast<int>(log((double)nbSize) / log(2.0) + 0.5); 

	condZeroP.resize(nBits-1);
	incrJoPdf.resize(nBits);
	entro.resize(nBits);
	condEntro.resize(nBits);
	


	for (int i = 0; i<nBits-1; i++)
	{
		condZeroP[i].resize((1<<(i+1)),0);

		incrJoPdf[i].resize((1<<(i+1)),0);
	}
	incrJoPdf[nBits-1].resize(nbSize,0);


	QtoBinMap.resize(nbSize,0);

	setQtoBinMap();

	convert();

	compEntro();
}

nbToBinPDFConv::nbToBinPDFConv(double* nbPdf_, int nbSize_)
{
	nbSize = nbSize_;
	nbPdf.resize(nbSize);
	for (int i=0; i<nbSize; i++)
		nbPdf[i] = nbPdf_[i];

	init();
}

void nbToBinPDFConv::setQtoBinMap()
{
	// use unsigned integer to represent bit mapping of element in GF(Q)
	// by default, use binary expression of the decimal number between 0-255
	// the lowerest bit (rightmost) is the first binary source
	for (unsigned int i = 0; i<(unsigned int)nbSize; i++)
		QtoBinMap[i] = i;
}

void nbToBinPDFConv::convert()
{

	// compute incrJoPdf
	for (int i = 0; i<nBits; i++)
	{
		for (int j = 0; j<nbSize; j++)
		{
			// QtoBinMap[j] & ((1<< (i+1))-1) is the last i+1 bits (in integer form) of the bit vector representation of j
			// it  represents the bits b0, b1, ..., bi
			incrJoPdf[i][QtoBinMap[j] & ((1<< (i+1))-1)] += nbPdf[j];
		}
	}

	// compute condZeroP
	for (int i=0; i<nBits-1; i++)
	{

		for (int j=0; j<(1<<(i+1)); j++)
		{
			//P(b_i+1 = 0 | (b0,b1,..., bi) = j(decimal form) ) = P(b_i+1 = 0, (b0,b1,...,b1) = j)
			//  / P( (b0,b1,...,bi) = j) = P( (b0,b1,..., bi, b_i+1) = j) / P((b0,b1,...bi) = j)
			condZeroP[i][j] = incrJoPdf[i+1][j] / incrJoPdf[i][j];
		}
	}
}


double nbToBinPDFConv::getZeroPGiven(int* knownBits, int len) const // lower bits should comes first in knownBits, i.e., b0 = knownBits[0]
{
	if (knownBits == NULL) // no conditioning, returns P(b0=0)
		return incrJoPdf[0][0];

	unsigned int index = 0; // index = b0, b1,...b_i-1 if len = i
	for (int i = 0; i<len; i++)
	{
		index = index << 1;
		index = index | knownBits[len-1 -i]; // set the lowest bit of index to knownBits[len-1-i], this will be moved to high bits by <<
		
	}

	return condZeroP[len-1][index];
}




void nbToBinPDFConv::compEntro()
{
	double ln2 = log(2.0);

	for (int i = 0; i<nBits; i++)
	{
		entro[i] = 0;
		for (int j = 0; j < (1<<(i+1)); j++)
		{
			if (incrJoPdf[i][j] < EPS)
				continue;
			entro[i] += - incrJoPdf[i][j] * log(incrJoPdf[i][j]) /ln2;
		}
	}


	condEntro[0] = entro[0];

	for (int i = 1; i<nBits; i++)
		condEntro[i] = entro[i] - entro[i-1];

}



}
