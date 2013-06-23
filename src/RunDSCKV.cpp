#include "TextFileInput.h"
#include "DSCKVSim.h"
#include "Util.h"
#include "StdHeaders.h"
#include <boost/shared_ptr.hpp>
#include "MTMCScheduler.h"
#include "Util.h"

using namespace MTMCSim;

void runMTDSCKVSimu()
{
	string workingDir,relative_profileFN, tempFN, formatFN;
	int queueSize, nThreads;

	Util::getConf("MTMCSimconf.txt", workingDir,relative_profileFN, tempFN, formatFN,queueSize, nThreads);

	string profilename = workingDir + relative_profileFN;
	TextFileInput paraIn(profilename);

	//simulation parameter map and computed result map are created here and everywhere else are references to here
	paraMap simuParaMap = paraIn.getParaMap();
	compResMap crMap;

	DOSDispUI simDispUI;

	//simulation parameter object is created here and COPIED into the simulator
	DSCKVSim::DSCKVPara simuPara;
	simuPara.loadPara(simuParaMap, simDispUI);
	simuPara.mulMode = 0;
	//simuPara.profileFN = profilename;
	//simuPara.dispUI = simDispUI;

	//UI objects are created here and simulation object refer to these
	
	TextSaveUI resFile(simuPara.saveFN);
	TextSaveUI tempFile(workingDir +tempFN);


	shared_ptr<DSCKVSim> psimulator(new DSCKVSim(simuPara, resFile, tempFile,simDispUI));


	//computed result object is owned by the simulation object and used to create compted result map
	psimulator->setCompResMap(crMap);

	// display and save the parameters and computed results
	Util::addCompRes(crMap, profilename, "profile file name");
	
	string ourStr = Util::getFormattedStr(simuParaMap, crMap, workingDir +formatFN);
	simDispUI.dispStr(ourStr);
	resFile.saveStr(ourStr);
	tempFile.saveStr(ourStr);

	MTMCScheduler sche(nThreads,queueSize, psimulator);

	sche.start();
}

void runMTDSCKVSimuMulMode()
{
	string workingDir,relative_profileFN, tempFN, formatFN;
	int queueSize, nThreads;

	Util::getConf("MTMCSimconf.txt", workingDir,relative_profileFN, tempFN, formatFN,queueSize, nThreads);

	string profilename = workingDir + relative_profileFN;
	TextFileInput paraIn(profilename);

	//simulation parameter map and computed result map are created here and everywhere else are references to here
	paraMap simuParaMap = paraIn.getParaMap();
	compResMap crMap;

	DOSDispUI simDispUI;

	//simulation parameter object is created here and COPIED into the simulator
	DSCKVSim::DSCKVPara simuPara;
	simuPara.loadPara(simuParaMap, simDispUI);
	simuPara.mulMode = 1;

	//UI objects are created here and simulation object refer to these
	
	TextSaveUI resFile(simuPara.saveFN);
	TextSaveUI tempFile(workingDir + tempFN);

	string mess = "In Multi Simulation Mode. Group simulation starts here.\n";
	simDispUI.dispStr(mess);
	resFile.saveStr(mess);
	tempFile.saveStr(mess);

	MTMCScheduler sche;
	sche.setNRunThreads(nThreads);
	sche.setRNQSize(queueSize);

	while (true)
	{
		shared_ptr<DSCKVSim> psimulator(new DSCKVSim(simuPara,resFile,tempFile, simDispUI));

		//computed result object is owned by the simulation object and used to create compted result map
		psimulator->setCompResMap(crMap);

		// display and save the parameters and computed results
		Util::addCompRes(crMap, profilename, "profile file name");
	
		string ourStr = Util::getFormattedStr(simuParaMap, crMap, workingDir +formatFN);
		simDispUI.dispStr(ourStr);
		resFile.saveStr(ourStr);
		tempFile.saveStr(ourStr);

		sche.setSimulator(psimulator);

		sche.start();

		if (psimulator->quit())
			break;
		else
		{
			simuPara = * (boost::static_pointer_cast<DSCKVSim::DSCKVPara, SimuPara>(psimulator->getNextPara()));
			string mess = "Current simulation terminated. New parameters will be use in the next simulation.\n";
			simDispUI.dispStr(mess);
			resFile.saveStr(mess);
			tempFile.saveStr(mess);

		}
	}

}


/*

int main()
{ 
	runMTDSCKVSimu();
	//runMTDSCKVSimuMulMode();

}
*/
