#include "SimuFactory.h"


namespace MTMCSim{
void SimuFactory::regSimulator(string name, createFcn creator){
	simuMap.insert(std::pair<string,createFcn>(name, creator));

}

shared_ptr<MTMCSimBase> SimuFactory::create(string name, SimuPara& para, TextSaveUI& resFile_, TextSaveUI& tempFile_, DispUI& dispUI_){
	map<string, createFcn>::iterator it = simuMap.find(name);
	if (it != simuMap.end())
		return ((*it).second)(para,resFile_,tempFile_,dispUI_);
	else
	{
		dispUI_.dispStr("Simulator " + name + " Not found!");
		return shared_ptr<MTMCSimBase>();
	}
}


//SimuFactory SimuFactory::theInstance;


}