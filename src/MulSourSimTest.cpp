#include "TextFileInput.h"
#include "MulSourSim.h"
#include "Util.h"
#include "StdHeaders.h"
#include <boost/shared_ptr.hpp>
#include "MTMCScheduler.h"
#include "DOSDispUI.h"

using namespace MTMCSim;

void runMulSourSimu()
{

	string profilename = "MulSourSim/In3S_1.txt";
	TextFileInput paraIn(profilename);

	//simulation parameter map and computed result map are created here and everywhere else are references to here
	paraMap simuParaMap = paraIn.getParaMap();
	compResMap crMap;

	DOSDispUI simDispUI;

	//simulation parameter object is created here and COPIED into the simulator
	MulSourSim::MulSourPara simuPara;
	simuPara.loadPara(simuParaMap, simDispUI);

	//UI objects are created here and simulation object refer to these
	
	TextSaveUI resFile(simuPara.saveFN);
	TextSaveUI tempFile("MulSourSim/MTMCtempFile.txt");


	MulSourSim simulator(simuPara,resFile,tempFile, simDispUI);
	
	//computed result object is owned by the simulation object and used to create compted result map
	simulator.setCompResMap(crMap);

	// display and save the parameters and computed results
	Util::addCompRes(crMap, profilename, "profile file name");
	
	simDispUI.dispStr(Util::getFormattedStr(simuParaMap, crMap, "MulSourSim/MulSourFormat.txt"));
	resFile.saveStr(Util::getFormattedStr(simuParaMap, crMap, "MulSourSim/MulSourFormat.txt"));

	simulator.singleThreadSim();

	simulator.postProcess();
}

void runMTMulSourSim()
{
	string profilename = "MulSourSim/In3S_1.txt";
	TextFileInput paraIn(profilename);

	//simulation parameter map and computed result map are created here and everywhere else are references to here
	paraMap simuParaMap = paraIn.getParaMap();
	compResMap crMap;

	DOSDispUI simDispUI;

	//simulation parameter object is created here and COPIED into the simulator
	MulSourSim::MulSourPara simuPara;
	simuPara.loadPara(simuParaMap, simDispUI);

	//UI objects are created here and simulation object refer to these
	
	TextSaveUI resFile(simuPara.saveFN);
	TextSaveUI tempFile("MulSourSim/MTMCtempFile.txt");


	shared_ptr<MulSourSim> psimulator(new MulSourSim(simuPara,resFile,tempFile, simDispUI));
	
	//computed result object is owned by the simulation object and used to create compted result map
	psimulator->setCompResMap(crMap);

	// display and save the parameters and computed results
	Util::addCompRes(crMap, profilename, "profile file name");
	
	simDispUI.dispStr(Util::getFormattedStr(simuParaMap, crMap, "MulSourSim/MulSourFormat.txt"));
	resFile.saveStr(Util::getFormattedStr(simuParaMap, crMap, "MulSourSim/MulSourFormat.txt"));

	MTMCScheduler sche(2,8,psimulator);

	sche.start();

}


//int main()
//{

	//runMulSourSimu();
	//runMTMulSourSim();
//}
