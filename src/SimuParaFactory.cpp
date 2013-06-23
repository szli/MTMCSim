#include "SimuParaFactory.h"


namespace MTMCSim{
void SimuParaFactory::regSimuPara(string name, createFcn creator){
	simuParaMap.insert(std::pair<string,createFcn>(name, creator));
	
	
}

shared_ptr<SimuPara> SimuParaFactory::create(string name){
	map<string, createFcn>::iterator it = simuParaMap.find(name);
	if (it != simuParaMap.end())
		return ((*it).second)();
	else
	{
		cout << ("Simulator " + name + " Not found!") << endl;
		return shared_ptr<SimuPara>();
	}
}


//SimuParaFactory SimuParaFactory::theInstance;


}