#include "KVChanSim.h"

namespace MTMCSim{



KVChanSim::KVChanSim(KVChanPara& para_,TextSaveUI& resFile_, TextSaveUI& tempFile_, DispUI& dispUI_)
    : nbSize(para_.nbSize), codeLen(para_.codeLen), infoLen(para_.infoLen),
    maxErrFr(para_.maxErrFr), lambda(para_.lambda),
    uniGen(para_.seed), uniIntGen(0,para_.nbSize-1), para(para_), 
resFile(resFile_), tempFile(tempFile_), dispUI(dispUI_),
MultiRunSimuBase(para_.mulMode,para_), channel(para_.SNR, (double)para_.infoLen/(double)para_.codeLen)
{
    nBits = static_cast<int>(log((double)nbSize) / log(2.0) + 0.5);
    BPSKmodem.setBitPerSym(nBits);
    getCompRes();
    
}

void KVChanSim::getCompRes()
{

    time_t timep;
    time (&timep);

    #define _CRT_SECURE_NO_WARNINGS
    compRes.startTime = asctime(localtime(&timep));
}


int KVChanSim::KVLooseBound(const int &cost) // Corollary 5 in KV paper
{
    return((int)sqrt((double) 2 * (infoLen - 1) * cost));
}


int KVChanSim::KVTightBound(const int &cost) // Lemma 1 in KV paper
{
    int delta;
    int term1, term2;
    
    int looseB = KVLooseBound(cost);

    for(delta=0; delta<=looseB; delta++)
    {
        if ((delta + 1) % (infoLen - 1) == 0)
        {
            term1 = (delta + 1) / (infoLen - 1);
        }
        else
            term1 = (int) floor((double) (delta + 1) / (infoLen - 1)) + 1;
        
        term2 = (int) floor((double) delta / (infoLen - 1));

        if ((term1 * (delta + 1 - (infoLen - 1) * term2 / 2)) > cost)
        {
            return(delta);
        }
    }
    return(MAX_VALUE); // Added Nov.4 2009, So that all control paths have a return, this statement should not be executed. 
}

void KVChanSim::oneStep(const std::vector<int> &xSeq,const vector<double>& normRV, ThreadWorkSpace& ws)
{
    int cost = 0, score = 0;

    BPSKmodem.modulate(xSeq, ws.channelIn);
    channel.pass(ws.channelIn,ws.channelOut, normRV);
    BPSKmodem.demGetReli(ws.channelOut, ws.reliMat, channel.getSigma());

    //Hard decision
    vector<int> HDseq(codeLen);

    for (int i=0; i<codeLen; i++)
    {
        double maxP = 0;
        for (int j = 0; j<nbSize; j++)
        {
            if (ws.reliMat[j][i] > maxP)
            {
                maxP = ws.reliMat[j][i];
                HDseq[i] = j;
            }
        }
    }

    int numE = 0;
    for (int i = 0; i<codeLen; i++)
        numE += (HDseq[i]!= xSeq[i]);

    // If want to apply this to long codes, the type of cost should be changed to double
    // long code ususllly uses larger lambda
    // MulMat[i][j] : The jth location the ith field element, Q rows n columns
    cost = 0;
    // Propotional assignment
    /*for (int i = 0; i<nbSize; i++)
        for (int j = 0; j<codeLen;j++)
        {
            ws.MulMat[i][j] = static_cast<int>(floor(para.lambda * ws.reliMat[i][j]));
            cost += (ws.MulMat[i][j] + 1)*(ws.MulMat[i][j]) /2;
        }
        */

    //For Hemitian code 
    int g = (int)(sqrt((double)nbSize)*(sqrt((double)nbSize)-1) / 2);
    cost = KVMul(ws.reliMat, ws.MulMat, (int)para.lambda, 1, infoLen-1 + g );
    score = 0;
    for (int i = 0; i<codeLen; i++)
        score += ws.MulMat[xSeq[i]][i];
    int RHS = KVTightBound(cost);
    ws.isFSucc = (score >RHS);

}

// For Hemitian Code
int KVChanSim::HemitianLB(int cost)
{
    int g = (int)(sqrt((double)nbSize)*(sqrt((double)nbSize)-1) / 2);
    return((int)sqrt((double) 2 * (infoLen - 1 + g) * cost) + g);
}

int KVChanSim::KVMul(vector<vector<double> >& reliMat, vector<vector<int> >& MulMat, int para, int mode, int yweight)
{
    for (int i = 0; i<nbSize; i++)
        for (int j = 0; j<codeLen; j++)
            MulMat[i][j] = 0;

    int currentCost = 0;
    while (true)
    {
        double max = 0;
        int maxi = 0,maxj = 0;

        for (int i = 0; i<nbSize; i++)
            for (int j = 0; j<codeLen; j++)
                if (reliMat[i][j] > max)
                {
                    max = reliMat[i][j];
                    maxi = i;
                    maxj = j;
                }
        reliMat[maxi][maxj] = reliMat[maxi][maxj] / (MulMat[maxi][maxj] + 2);
        MulMat[maxi][maxj] ++;
        currentCost += MulMat[maxi][maxj];

        if (mode == 0) // limit sum mode, para  = total multiplicity
        {
            if ((--para) == 0)
                break;
        }
        else if (mode == 1) // limit the y degree, para  = max y degree
        {
            if (KVLooseBound(currentCost) / yweight > para)
            {
                currentCost -= MulMat[maxi][maxj];
                MulMat[maxi][maxj] -- ;
                
                break;
            }
        }
    }
    return currentCost;
}

            

    



void KVChanSim::oneThread(shared_ptr<RNObj> pobj, const MCThreadInfo& threadInfo)
{
    ThreadWorkSpace& ws = *(workspaceVec[threadInfo.threadID]);
    
    shared_ptr<XNSeq> pxnseq = boost::dynamic_pointer_cast<XNSeq, RNObj>(pobj);

    oneStep(pxnseq->xSeq, pxnseq->normRV, ws);

}

void KVChanSim::dumpOneResult(const MCThreadInfo& threadInfo)
{
    ThreadWorkSpace& ws = *(workspaceVec[threadInfo.threadID]);

        // update stat
        simuRes.totFr ++;
        
        if (!ws.isFSucc)
            simuRes.errFr ++;

        simuRes.FER = (double)simuRes.errFr / (double)simuRes.totFr;

        // one step output
        if (simuRes.totFr % para.saveTempPerFr == 0)
        {
            ostringstream ress;
            ress << "Frame: " << simuRes.totFr << " Err Fr: " << simuRes.errFr << 
                " FER: " << simuRes.FER << "\n";
            
            resFile.saveStr(ress.str());
            tempFile.saveStr(ress.str());
            
        }

        if (simuRes.totFr % para.dispPerFr == 0)
        {
            ostringstream ress;
            ress <<  "Frame: " << simuRes.totFr << " Err Fr: " << simuRes.errFr << 
                " FER: " << simuRes.FER << "\r";
            dispUI.dispStr(ress.str());
        }

}

void KVChanSim::postProcess()
{
        ostringstream ress;
        ress << "Simulation terminated.\n" << endl;
        ress << "Frame: " << simuRes.totFr << " Err Fr: " << simuRes.errFr << 
            " FER: " << simuRes.FER << "\n";
        ress.flush();
        resFile.saveStr(ress.str());

        ress.str(string());
        ress << "\nSimulation terminated.\n";
        ress <<  "Frame: " << simuRes.totFr << " Err Fr: " << simuRes.errFr << 
            " FER: " << simuRes.FER << endl;
        ress.flush();
        dispUI.dispStr(ress.str());
}

void KVChanSim::prepare(const MCThreadInfo& threadInfo)
{
    for (int i = 0; i<threadInfo.nThread; i++)
    {
        workspaceVec.push_back(shared_ptr<ThreadWorkSpace>
            (new ThreadWorkSpace(codeLen, codeLen*nBits,nbSize)));
    }
}

shared_ptr<RNObj> KVChanSim::genRNObj()
{
    shared_ptr<XNSeq> pseq(new XNSeq(codeLen, codeLen*nBits));
    for (int i = 0; i<codeLen; i++)
    {
        pseq->xSeq[i] = uniIntGen(uniGen);
    }

    //cout << uniGen() << " " << normGen(uniGen);
    boost::uniform_01<> uni01Gen;

    //Box Muller

    int mid = codeLen * nBits / 2;
    double R, U2;
    for (int i = 0; i<mid; i++)
    {
        R = sqrt(-2*log(uni01Gen(uniGen)));
        U2 = uni01Gen(uniGen);
        pseq->normRV[2*i] = R*cos(2*3.1415926535897932 * U2);
        pseq->normRV[2*i+1] = R*sin(2*3.1415926535897932 * U2);
    }
    //the last element if the length is odd.
    if (codeLen*nBits % 2 != 0)
    {
        R = sqrt(-2*log(uni01Gen(uniGen)));
        U2 = uni01Gen(uniGen);
        pseq->normRV[codeLen*nBits-1] = R*cos(2*3.1415926535897932 * U2);
    }


        


    //cout << uni01Gen(uniGen);
    /*for (int i = 0; i<codeLen*nBits; i++)
    {
        pseq->normRV[i] = normGen(uniGen);
        
    }*/

    return pseq;
}



void KVChanSim::KVChanPara::loadPara(const paraMap& para, DispUI& ui)
{
    using namespace Util;
    if (!getIntPara(para, "alphabet size", nbSize))
        paraMissingWarning(ui, "alphabet size");

    if (!getIntPara(para, "code length", codeLen))
        paraMissingWarning(ui, "code length");

    if (!getIntPara(para, "info length", infoLen))
        paraMissingWarning(ui, "info length");

    if (!getDouPara(para, "lambda", lambda))
        paraMissingWarning(ui,"lambda");

    if (!getDouPara(para,"SNR", SNR))
        paraMissingWarning(ui,"SNR");

    if (!getIntPara(para, "seed", seed))
        paraMissingWarning(ui, "seed");

    if (!getIntPara(para, "max error frame", maxErrFr))
        paraMissingWarning(ui, "max error frame");

    if (!getIntPara(para, "display freq", dispPerFr))
        paraMissingWarning(ui, "display freq");

    if (!getIntPara(para, "save freq", saveTempPerFr))
        paraMissingWarning(ui, "save freq");

    if (!getIntPara(para,"multiple run mode", mulMode))
        paraMissingWarning(ui, "multiple run mode");

    if (!getStrPara(para,"save file", saveFN))
        paraMissingWarning(ui,"save file");

    SIM_ASSERT(nbSize > codeLen, "The alphabet size should greater than the code length for Reed-Solomon code");
    SIM_ASSERT(codeLen > infoLen, "The code length should greater than the information length for Reed-Solomon code");


}


void KVChanSim::setCompResMap(compResMap& compResM_) const
{
    compRes.setMapping(compResM_);
}


void KVChanSim::KVChanCompRes::setMapping(compResMap &crMap) const
{
    Util::addCompRes(crMap, startTime, "Start time");
}
    
bool KVChanSim::checkTerm() 
{
    static vector<int> usedInfoLen;
    if (mulMode == 0) // not multi-run mode
    {
        usedInfoLen.clear();
        return (para.maxErrFr == simuRes.errFr);
    }
    else
    {
        return true;// TO DO
        
    }
}



}
