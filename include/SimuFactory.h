#ifndef _SIMUFACTORY_H_INCLUDED
#define _SIMUFACTORY_H_INCLUDED
#include "StdHeaders.h"
#include "SimuPara.h"
#include <boost/shared_ptr.hpp>
#include "MTMCSimBase.h"
#include "TextSaveUI.h"
#include "DOSDispUI.h"

using std::map;
namespace MTMCSim{

/** Factory class for classes derived from MTMCSimBase, the simulators. \ Singleton. 
*/
class SimuFactory{
public:
	typedef shared_ptr<MTMCSimBase> (*createFcn)(SimuPara&, TextSaveUI&,  TextSaveUI& , DispUI& );
	static SimuFactory& instance()
	{	static SimuFactory theInstance;
		return theInstance; }
	
	void regSimulator(string name, createFcn);
	shared_ptr<MTMCSimBase> create(string name, SimuPara& para, TextSaveUI& resFile_, TextSaveUI& tempFile_, DispUI& dispUI_);


private:
	//static SimuFactory theInstance;
	map<string,createFcn> simuMap;

};
}

#endif