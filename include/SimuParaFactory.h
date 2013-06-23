#ifndef _SIMUPARAFACTORY_H_INCLUDED
#define _SIMUPARAFACTORY_H_INCLUDED
#include "SimuPara.h"
#include <boost/shared_ptr.hpp>
#include "StdHeaders.h"

namespace MTMCSim{

/** Factory class for classes derived from SimuPara, the simulation parameter objects. \ Singleton. 
*/
class SimuParaFactory{
	public:
	typedef shared_ptr<SimuPara> (*createFcn)();
	static SimuParaFactory& instance()
	{ static SimuParaFactory theInstance;
		return theInstance; }
	
	void regSimuPara(string name, createFcn);
	shared_ptr<SimuPara> create(string name);


private:
	
	map<string,createFcn> simuParaMap;
	
};
}
#endif
