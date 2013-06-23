#ifndef _PDFXY_H_INCLUDED
#define _PDFXY_H_INCLUDED
#include "StdHeaders.h"

namespace MTMCSim{

/** Handles joint pmf for two sources. Reads pmf, compute entropies. (We do not distinguish pmf and pdf here.) copy safe */
class pdfXY{
public:
	pdfXY()
	{}
	/** Construct from marginal pdf of X and the conditional pdf P(Y|X).
	*	@param xMarPdf_ the ith element = p(x = i)
	*	@param yxCondPdf the ith row, jth col:P( Y = j|X = i),  yxCondPdf_ should have same number of cols in each row.
	*/
	pdfXY(vector<double>& xMarPdf_, vector<vector<double> >& yxCondPdf_)
	{ setPdf(xMarPdf_, yxCondPdf_);}

	/** Construct from joint pdf P(X,Y) */
	pdfXY(vector<vector<double> >& joPdf_)
	{ setPdf(joPdf_); }

	/** Q-ary symmetric correlation model.*/
	pdfXY(double agreeProb, int size)
	{ setPdf(agreeProb, size); }

	/** Pdf given in a file.*/
	pdfXY(string pdffn, int xSize, int ySize)
	{ setPdf(pdffn, xSize, ySize); }

	double getJoH() const
	{ return joH;}

	double getHX() const
	{ return HX;}

	double getHY() const
	{ return HY;}

	double getCondHYX() const
	{ return condHYX;}

	double getCondHXY() const
	{ return condHXY;}

	double getIXY() const
	{ return IXY;}

	/** The returned reference to object has the same life time as the pdfXY object. */
	const vector<vector<double> >& getJoPdf() const
	{ return joPdf; }

	/** The returned reference to object has the same life time as the pdfXY object. */
	const vector<double>& getXMarPdf() const
	{ return xMarPdf;}

	/** The returned reference to object has the same life time as the pdfXY object. */
	const vector<vector<double> >& getYXCondPdf() const
	{ return yxCondPdf;}


	/** yxCondPdf should have same number of cols in each row. */
	void setPdf(vector<double>& xMarPdf_, vector<vector<double> >& yxCondPdf_);

	void setPdf(vector<vector<double> >& joPdf_);

	void setPdf(double agreeProb, int size);

	void setPdf(string pdffn, int xSize, int ySize);
	

private:
	unsigned int abXSize, abYSize;
	double joH, HX, HY, condHYX, condHXY, IXY;
	vector<vector<double> > joPdf;//(i,j)th: P(X=i,Y=j)
	vector<double> xMarPdf;
	vector<double> yMarPdf;
	vector<vector<double> > yxCondPdf; //ith row, jth col:P( Y = j|X = i), 

	void compEntro();


};

}
#endif