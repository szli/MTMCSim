#include "MRSimuFactory.h"

namespace MTMCSim{
void MRSimuFactory::regSimulator(string name, createFcn creator){
	simuMap.insert(std::pair<string,createFcn>(name, creator));

}

shared_ptr<MultiRunSimuBase> MRSimuFactory::create(string name, SimuPara& para, TextSaveUI& resFile_, TextSaveUI& tempFile_, DispUI& dispUI_){
	map<string, createFcn>::iterator it = simuMap.find(name);
	if (it != simuMap.end())
		return ((*it).second)(para,resFile_,tempFile_,dispUI_);
	else
	{
		dispUI_.dispStr("Simulator " + name + " Not found!");
		return shared_ptr<MultiRunSimuBase>();
	}
}


}