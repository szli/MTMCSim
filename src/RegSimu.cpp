
#include "SimuItemCreator.h"

#include "SimuParaItemCreator.h"
#include "SimuPara.h"
#include "MTMCSimBase.h"

#include "DSCKVSim.h"
#include "MulStageSim.h"
#include "MulSourSim.h"
#include "MRSimuItemCreator.h"

#include "KVChanSim.h"
#include "MSALDPC.h"

#include "DSCKVSimRate.h"
#include "MulStageSimRate.h"
#include "MulSourSimRate.h"

namespace MTMCSim{

SimuParaItemCreator<DSCKVSim::DSCKVPara> DSCKVParaItem("DSCKVSim");

SimuItemCreator<DSCKVSim, DSCKVSim::DSCKVPara> DSCKVSimItem("DSCKVSim");


SimuParaItemCreator<MulStageSim::MulStagePara> MulStageParaItem("MulStageSim");

SimuItemCreator<MulStageSim, MulStageSim::MulStagePara> MulStageSimItem("MulStageSim");

SimuParaItemCreator<MSALSim::MSALPara> MSALParaItem("MSALSim");

SimuItemCreator<MSALSim, MSALSim::MSALPara> MSALSimItem("MSALSim");

SimuParaItemCreator<MulSourSim::MulSourPara> MulSourParaItem("MulSourSim");

SimuItemCreator<MulSourSim, MulSourSim::MulSourPara> MulSourSimItem("MulSourSim");

MRSimuItemCreator<DSCKVSim, DSCKVSim::DSCKVPara> DSCKVMRSimItem("DSCKVSim");

SimuParaItemCreator<KVChanSim::KVChanPara> KVChanParaItem("KVChanSim");

SimuItemCreator<KVChanSim, KVChanSim::KVChanPara> KVChanSimItem("KVChanSim");

SimuParaItemCreator<DSCKVSimRate::DSCKVRatePara> DSCKVRateParaItem("DSCKVSimRate");

SimuItemCreator<DSCKVSimRate, DSCKVSimRate::DSCKVRatePara> DSCKVSimRateSimItem("DSCKVSimRate");

SimuParaItemCreator<MulStageSimRate::MulStagePara> MSRateParaItem("MulStageSimRate");

SimuItemCreator<MulStageSimRate,MulStageSimRate::MulStagePara> MSRateSimItem("MulStageSimRate");

SimuItemCreator<MulSourSimRate, MulSourSimRate::MulSourParaRate> MSourRateSimItem("MulSourSimRate");

SimuParaItemCreator<MulSourSimRate::MulSourParaRate> MSourRateParaItem("MulSourSimRate");

}