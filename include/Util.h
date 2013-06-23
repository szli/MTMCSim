#ifndef _UTIL_H_INCLUDED
#define _UTIL_H_INCLUDED

#include "StdHeaders.h"
#include "DispUI.h"


namespace MTMCSim{
namespace Util{

	inline string intToStr(int v)
	{
	stringstream temp;
	temp << v;
	return temp.str();

	}

	inline string douToStr(double v)
	{
		stringstream temp;
		temp << v;
		return temp.str();
	}
	
	/** get integer parameter from the paremeter map, returns whether the given parameter name is in the map. */
	bool getIntPara(const paraMap& para, const string& name, int& data);
	/** get double parameter from the paremeter map, returns whether the given parameter name is in the map. */
	bool getDouPara(const paraMap& para, const string& name, double& data);
	/** get string parameter from the paremeter map, returns whether the given parameter name is in the map. */
	bool getStrPara(const paraMap& para, const string& name, string& data);
	
	inline void paraMissingWarning(DispUI& ui, const string& name)
	{ ui.dispStr("Warning! Parameter " + name + " Missing! Using default value"); }

	/** add a computed result of int type to the computed result map.*/
	inline void addCompRes(compResMap& compRes, const int value, const string& name)
	{ compRes[name] = intToStr(value); }
	/** add a computed result of double type to the computed result map.*/
	inline void addCompRes(compResMap& compRes, const double value,const string& name)
	{	compRes[name] = douToStr(value); }
	/** add a computed result of string type to the computed result map.*/
	inline void addCompRes(compResMap& compRes, const string& value,const string& name)
	{	compRes[name] = value;	}

	/** Put the parameters and computed results into a formatted string whose format is given by a text file. */
	string getFormattedStr( const paraMap& inputPara, const compResMap& compRes, string formatFN);

	/** Binary Entropy function. Make sure 0<= p <= 1. No check inside the function. */
	inline double bEntro(double p)
	{
	if ((p < EPS)||(1-p < EPS))
		return 0;
	else
		return (p * log((double)1/p) + (1-p) * log((double)1/(1-p))) / log(2.0);
	}

	void getConf(const string& confFN, string& workingDir,string&relative_profileFN,string& tempFN,string& formatFN, int& queueSize, int& nThreads);
}
}

#endif