#ifndef _TOPPARA_H_INCLUDED
#define _TOPPARA_H_INCLUDED
#include "SimuPara.h"
#include "Util.h"

namespace MTMCSim{
class topConf:public SimuPara
{ 
public:
	string workingDir,relative_profileFN, tempFN, formatFN;
	int queueSize, nThreads;
	string simuName;
	int mulMode;///< 0 if not in multi-run mode, nonzero if in multi-run mode
	int noSimu;///< 1, do not run simu, just display configurations, 0, run simu
	
	shared_ptr<SimuPara> clone() const
	{ return shared_ptr<topConf>(new topConf(*this));}
	

	void loadPara(const paraMap& para, DispUI& ui)
	{
		using namespace Util;
		if (!getIntPara(para, "nThreads", nThreads))
			paraMissingWarning(ui, "nThreads");

		if (!getIntPara(para, "queueSize", queueSize))
			paraMissingWarning(ui, "queueSize");

		if (!getStrPara(para,"profile file name", relative_profileFN))
			paraMissingWarning(ui,"profile file name");
		string dir = relative_profileFN.substr(0,relative_profileFN.find('/'));

		if (!getStrPara(para,"format file name",  formatFN))
		{
			paraMissingWarning(ui,"format file name");
			formatFN = dir + "/Format.txt";
		}
		
		if (!getStrPara(para,"temp file name", tempFN))
		{
			paraMissingWarning(ui,"temp file name");
			tempFN = dir + "/tempFile.txt";
		}

		if (!getStrPara(para,"simulation name", simuName))
			paraMissingWarning(ui,"simulation name");
		if (!getIntPara(para, "multiple run mode", mulMode))
			paraMissingWarning(ui,"multiple run mode");
		if (!getIntPara(para, "no simulation", noSimu))
			paraMissingWarning(ui,"no simulation");
		
		


	}
};

}





#endif
