#ifndef _MRSIMUFACTORY_H_INCLUDED
#define _MRSIMUFACTORY_H_INCLUDED
#include "StdHeaders.h"
#include "SimuPara.h"
#include <boost/shared_ptr.hpp>
#include "MultiRunSimuBase.h"
#include "TextSaveUI.h"
#include "DOSDispUI.h"

using std::map;
namespace MTMCSim{

/** Factory class for classes derived from MultiRunSimuBase, i.e., simulations that supports multi-run mode. \ Singleton. 
*/
class MRSimuFactory{
public:
	typedef shared_ptr<MultiRunSimuBase> (*createFcn)(SimuPara&, TextSaveUI&,  TextSaveUI& , DispUI& );
	static MRSimuFactory& instance()
	{	static MRSimuFactory theInstance;
		return theInstance; }
	
	void regSimulator(string name, createFcn);
	shared_ptr<MultiRunSimuBase> create(string name, SimuPara& para, TextSaveUI& resFile_, TextSaveUI& tempFile_, DispUI& dispUI_);


private:
	map<string,createFcn> simuMap;

};
}

#endif