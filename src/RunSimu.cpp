#include "SimuFactory.h"
#include "SimuParaFactory.h"
#include "SimuPara.h"
#include "MTMCSimBase.h"
#include "MTMCScheduler.h"
#include "Util.h"
#include "TextFileInput.h"
#include "topConf.h"
#include "MRSimuFactory.h"
#include "MultiRunSimuBase.h"
using namespace MTMCSim;


void runMTSingleSimu(topConf& topConfs)
{
    DOSDispUI simDispUI;
    //TextFileInput topConfIn("MTMCSimconf.txt");
    //paraMap topConfMap = topConfIn.getParaMap();
    //topConf topConfs;
    //topConfs.loadPara(topConfMap,simDispUI);

    string profilename = topConfs.workingDir + topConfs.relative_profileFN;
    TextFileInput paraIn(profilename);

    //simulation parameter map and computed result map are created here and everywhere else are references to here
    paraMap simuParaMap = paraIn.getParaMap();
    compResMap crMap;

    string simuName = topConfs.simuName;
    shared_ptr<SimuPara> psimuPara = (SimuParaFactory::instance()).create(simuName);
    psimuPara->loadPara(simuParaMap,simDispUI);
    
    TextSaveUI resFile(psimuPara->saveFN);
    TextSaveUI tempFile(topConfs.workingDir + topConfs.tempFN);

    shared_ptr<MTMCSimBase> psimulator = (SimuFactory::instance()).create(simuName, *psimuPara, resFile,tempFile, simDispUI);

    //computed result object is owned by the simulation object and used to create compted result map
    psimulator->setCompResMap(crMap);

    // display and save the parameters and computed results
    Util::addCompRes(crMap, profilename, "profile file name");
    
    string ourStr = Util::getFormattedStr(simuParaMap, crMap, topConfs.workingDir + topConfs.formatFN);
    simDispUI.dispStr(ourStr);
    resFile.saveStr(ourStr);
    tempFile.saveStr(ourStr);
    if (topConfs.noSimu == 0)

    {
        MTMCScheduler sche(topConfs.nThreads,topConfs.queueSize, psimulator);

        sche.start();
    }
}

void runMTMulSimu(topConf& topConfs)
{
    if (topConfs.noSimu == 1)
        throw runtime_error("noSimu cannot be true in MTMulSimu mode");

    DOSDispUI simDispUI;
    /*TextFileInput topConfIn("MTMCSimconf.txt");
    paraMap topConfMap = topConfIn.getParaMap();
    topConf topConfs;
    topConfs.loadPara(topConfMap,simDispUI);*/

    string profilename = topConfs.workingDir + topConfs.relative_profileFN;
    TextFileInput paraIn(profilename);

    //simulation parameter map and computed result map are created here and everywhere else are references to here
    paraMap simuParaMap = paraIn.getParaMap();
    if ((simuParaMap.find("multiple run mode") != simuParaMap.end()) && (simuParaMap["multiple run mode"] == "0"))
    throw runtime_error("multiple run mode parameter not consisitent!");
    compResMap crMap;

    string simuName = topConfs.simuName;
    shared_ptr<SimuPara> psimuPara = (SimuParaFactory::instance()).create(simuName);
    psimuPara->loadPara(simuParaMap,simDispUI);
    
    TextSaveUI resFile(psimuPara->saveFN);
    TextSaveUI tempFile(topConfs.workingDir + topConfs.tempFN);


    string mess = "In Multi Simulation Mode. Group simulation starts here.\n";
    simDispUI.dispStr(mess);
    resFile.saveStr(mess);
    tempFile.saveStr(mess);

    MTMCScheduler sche;
    sche.setNRunThreads(topConfs.nThreads);
    sche.setRNQSize(topConfs.queueSize);

    while (true)
    {
        shared_ptr<MultiRunSimuBase> psimulator = (MRSimuFactory::instance()).create(simuName, *psimuPara, resFile,tempFile, simDispUI);

        //computed result object is owned by the simulation object and used to create compted result map
        psimulator->setCompResMap(crMap);

        // display and save the parameters and computed results
        Util::addCompRes(crMap, profilename, "profile file name");
    
        string ourStr = Util::getFormattedStr(simuParaMap, crMap, topConfs.workingDir + topConfs.formatFN);
        simDispUI.dispStr(ourStr);
        resFile.saveStr(ourStr);
        tempFile.saveStr(ourStr);

        sche.setSimulator(psimulator);

        sche.start();

        if (psimulator->quit())
            break;
        else
        {
            psimuPara = psimulator->getNextPara();
            string mess = "Current simulation terminated. New parameters will be use in the next simulation.\n";
            simDispUI.dispStr(mess);
            resFile.saveStr(mess);
            tempFile.saveStr(mess);

        }
    }
}

void runMTSimu(const string& topconffn)
{
    DOSDispUI simDispUI;
    TextFileInput topConfIn(topconffn);
    paraMap topConfMap = topConfIn.getParaMap();
    topConf topConfs;
    topConfs.loadPara(topConfMap,simDispUI);
    if (topConfs.mulMode == 0)
        runMTSingleSimu(topConfs);
    else
        runMTMulSimu(topConfs);
}

//Single thread simulation
void runSTSimu(const string& topconffn)
{
    DOSDispUI simDispUI;
    TextFileInput topConfIn(topconffn);
    paraMap topConfMap = topConfIn.getParaMap();
    topConf topConfs;
    topConfs.loadPara(topConfMap,simDispUI);

    string profilename = topConfs.workingDir + topConfs.relative_profileFN;
    TextFileInput paraIn(profilename);

    //simulation parameter map and computed result map are created here and everywhere else are references to here
    paraMap simuParaMap = paraIn.getParaMap();
    compResMap crMap;

    string simuName = topConfs.simuName;
    shared_ptr<SimuPara> psimuPara = (SimuParaFactory::instance()).create(simuName);
    psimuPara->loadPara(simuParaMap,simDispUI);
    
    TextSaveUI resFile(psimuPara->saveFN);
    TextSaveUI tempFile(topConfs.workingDir + topConfs.tempFN);

    shared_ptr<MTMCSimBase> psimulator = (SimuFactory::instance()).create(simuName, *psimuPara, resFile,tempFile, simDispUI);

    //computed result object is owned by the simulation object and used to create compted result map
    psimulator->setCompResMap(crMap);

    // display and save the parameters and computed results
    Util::addCompRes(crMap, profilename, "profile file name");
    
    string ourStr = Util::getFormattedStr(simuParaMap, crMap, topConfs.workingDir + topConfs.formatFN);
    simDispUI.dispStr(ourStr);
    resFile.saveStr(ourStr);
    tempFile.saveStr(ourStr);
    if (topConfs.noSimu == 0)
    {
        psimulator->singleThreadSim();
    }

}
    


int main(int argc, char* argv[] ){

    string topconffn;
    if (argc > 2)
    {
        cout << "Argument error! At most one argument!" << endl;
        exit(-1);
    }
    if (argc == 1)
    {
        cout << "Must specify a top configuration file name! Exiting..." << endl;
        exit(-1);
    }
    else
    {
        topconffn = argv[1];
        if (topconffn.find(".txt") == string::npos)
            topconffn.append(".txt");
    }


    runMTSimu(topconffn);
    //runSTSimu(topconffn);
}
