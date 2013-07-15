#include "DSCKVSim.h"

namespace MTMCSim{

extern const string default_pdfdir = "PDFfiles/";

DSCKVSim::DSCKVSim(DSCKVPara& para_,TextSaveUI& resFile_, TextSaveUI& tempFile_, DispUI& dispUI_)
    : nbSize(para_.nbSize), codeLen(para_.codeLen), infoLen(para_.infoLen),
    maxErrFr(para_.maxErrFr), lambda(para_.lambda),
uniGen(para_.seed), uni01Gen(uniGen, uni_01), para(para_), 
resFile(resFile_), tempFile(tempFile_), dispUI(dispUI_),
MultiRunSimuBase(para_.mulMode,para_)
{
    // Initialize the joint pdf object
    if (para_.pdfPara[0] == '#') // generate pdf in C++ program
    {
        if (para_.pdfPara.find("QARY") != para_.pdfPara.npos) // if "QARY" is in pdfPara
        {
            double agreeProb = atof(para_.pdfPara.substr(para_.pdfPara.rfind(' ')+1).c_str());
            nbJoPdf.setPdf(agreeProb, nbSize);
        }
    }
    else // import pdf from a file 
    {
        nbJoPdf.setPdf(para_.pdfPara, nbSize, nbSize); // pdfParaStr is the file name
    }

    getCompRes();

    

    //Init the RNG
    initRNG();
    
    const vector<vector<double> >* pyxCondPdf = 0;
    pdfXY decPdf;
    if (!para_.decPdfPara.empty())
    {
    

        if (para_.decPdfPara[0] == '#') // generate pdf in C++ program
        {
            if (para_.decPdfPara.find("QARY") != para_.decPdfPara.npos) // if "QARY" is in pdfPara
            {
                double agreeProb = atof(para_.decPdfPara.substr(para_.decPdfPara.rfind(' ')+1).c_str());
                decPdf.setPdf(agreeProb, nbSize);
            }
        }
        else // import pdf from a file 
            {
                decPdf.setPdf(para_.decPdfPara, nbSize, nbSize); // pdfParaStr is the file name
            }
    
        pyxCondPdf = &decPdf.getYXCondPdf();
    }
    else
        pyxCondPdf = &nbJoPdf.getYXCondPdf();


    //Initialize the pre-determined multiplicity matrix
    //Q-by-Q matrix
    PreMul.resize(nbSize,vector<int>(nbSize,0));
    //const vector<vector<double> >& yxCondPdf = nbJoPdf.getYXCondPdf();



    for (int i = 0; i<nbSize; i++)
        for (int j = 0; j<nbSize; j++)
            PreMul[i][j] = static_cast<int>(floor(lambda*(*pyxCondPdf)[i][j]));
    // PreMul[i] is the "column" (in KV sense) of multiplicity matrix when X = i

//#define MODIFY_MUL

//  #ifdef MODIFY_MUL
//  //for qary symmetric 0.6 only
//  for (int i = 1; i<nbSize-1; i++)
//  {
//      PreMul[i][i-1] = static_cast<int>(floor(PreMul[i][i-1] * 1.1));
//      PreMul[i][i+1] = static_cast<int>(floor(PreMul[i][i+1] *  1.1));
//  }
//
//  for (int i = 0; i<nbSize; i++)
//      PreMul[i][i] = PreMul[i][i] * 0.9;
//  PreMul[0][1] = static_cast<int>(floor(PreMul[0][1] *  1.1));
//  PreMul[0][255] = static_cast<int>(floor(PreMul[0][255] *  1.1));
//  PreMul[255][254] = static_cast<int>(floor(PreMul[255][254] *  1.1));
//  PreMul[255][0] = static_cast<int>(floor(PreMul[255][1] *  1.1));
//
//
//#endif
}

void DSCKVSim::initRNG()
{
    switch (para.rngTableOpt)
    {
        case 0 :// regularly compute the tables
                xDistr.setDistribution(nbSize, nbJoPdf.getXMarPdf());

                yxCondDistr.resize(nbSize);

                for ( int i = 0; i<nbSize; i++)
                    yxCondDistr[i].setDistribution(nbSize, nbJoPdf.getYXCondPdf()[i]);

                break;
        case 1: // load the table
            {
                xDistr.setDistributionWOGenTab(nbSize, nbJoPdf.getXMarPdf());

                yxCondDistr.resize(nbSize);

                for ( int i = 0; i<nbSize; i++)
                    yxCondDistr[i].setDistributionWOGenTab(nbSize, nbJoPdf.getYXCondPdf()[i]);

                string tablefn;
                if (para.pdfPara[0] == '#')
                    tablefn = default_pdfdir + para.pdfPara.substr(1) + "Q" + Util::intToStr(nbSize) + ".txt";
                else
                {
                    tablefn = para.pdfPara;
                    tablefn.insert(tablefn.find_first_of('.'),"Tab"); // XXXXXTab.txt
                }

                ifstream tabfile(tablefn.c_str());
                string oneline;

                getline(tabfile,oneline); // the first line is for the distribution of X
                xDistr.loadTab(oneline);

                for (int i = 0; i<nbSize; i++)
                {
                    getline(tabfile,oneline);
                    yxCondDistr[i].loadTab(oneline);
                }
                tabfile.close();

                break;
            } // scope needed for variable definition

        case 2: // compute and save the table
            {

                string tablefn;
                if (para.pdfPara[0] == '#')
                    tablefn = default_pdfdir + para.pdfPara.substr(1)   + "Q" + Util::intToStr(nbSize) + ".txt";
                else
                {
                    tablefn = para.pdfPara;
                    tablefn.insert(tablefn.find_first_of('.'),"Tab"); // XXXXXTab.txt
                }

                xDistr.setDistribution(nbSize, nbJoPdf.getXMarPdf());

                yxCondDistr.resize(nbSize);

                for ( int i = 0; i<nbSize; i++)
                    yxCondDistr[i].setDistribution(nbSize, nbJoPdf.getYXCondPdf()[i]);

                

                ofstream otabfile(tablefn.c_str());

                string tabstr;
                xDistr.tabToStr(tabstr);
                otabfile << tabstr << endl;// each row has one table (in fact two, p and ialt)

                for (int i=0; i<nbSize; i++)
                {
                    yxCondDistr[i].tabToStr(tabstr);
                    otabfile << tabstr << endl;
                }

                otabfile.close();
                break;
            }// scope needed for variable definition
    }// end of switch
}
void DSCKVSim::getCompRes()
{
    compRes.HX = nbJoPdf.getHX();
    compRes.HY = nbJoPdf.getHY();
    compRes.HYX = nbJoPdf.getCondHYX();
    compRes.infoLen = infoLen;

    
    int nBits = static_cast<int>(log((double)nbSize) / log(2.0) + 0.5); 

    compRes.txRate = (codeLen - infoLen)/static_cast<double>(codeLen) * nBits;

    compRes.gap = compRes.txRate - compRes.HYX;

    time_t timep;
    time (&timep);

    #define _CRT_SECURE_NO_WARNINGS
    compRes.startTime = asctime(localtime(&timep));
}


double DSCKVSim::KVLooseBound(const double &cost) // Corollary 5 in KV paper
{
    return(sqrt((double) 2 * (infoLen - 1) * cost));
}


double DSCKVSim::KVTightBound(const double &cost) // Lemma 1 in KV paper
{
    double delta;
    double term1, term2;
    
    double looseB = KVLooseBound(cost);

    for(delta=0; delta<=looseB; delta++)
    {
        
        term1 = ceil((delta + 1) / (infoLen - 1));
        
        term2 = (int) floor((double) delta / (infoLen - 1));

        if ((term1 * (delta + 1 - (infoLen - 1) * term2 / 2)) > cost)
        {
            return(delta);
        }
    }
    return(MAX_VALUE); // Added Nov.4 2009, So that all control paths have a return, this statement should not be executed. 
}

void DSCKVSim::oneStep(const std::vector<int> &xSeq, const std::vector<int> &ySeq, ThreadWorkSpace& ws)
{
    double cost = 0, score = 0;

#define KVDSC
#define HermitianKVDSC
#define GSDSC
#define HermitianGSDSC

    if (para.algoMode == "KVDSC")
    {


#ifdef KVDSC
    // KV DSC
    for (int i = 0; i<codeLen; i++) // the ith "column" of M, in fact, here it is stored in a row
        for (int j = 0; j<nbSize; j++)
            cost += PreMul[xSeq[i]][j] * (PreMul[xSeq[i]][j] + 1) / 2;

    for (int i = 0; i<codeLen; i++)
            score += PreMul[xSeq[i]][ySeq[i]];

    double RHS = KVTightBound(cost);
    //cout <<"cost:" << cost << "score:"<< score << "RHS:"<< RHS<<endl;
    //cout << endl;
    ws.isFSucc = (score > RHS);
#endif
    }
    else if (para.algoMode == "HermitianKVDSC")
    {
#ifdef HermitianKVDSC
    for (int i = 0; i<codeLen; i++) // the ith "column" of M, in fact, here it is stored in a row
        for (int j = 0; j<nbSize; j++)
            cost += PreMul[xSeq[i]][j] * (PreMul[xSeq[i]][j] + 1) / 2;

    for (int i = 0; i<codeLen; i++)
            score += PreMul[xSeq[i]][ySeq[i]];

    int RHS = HemitianLB(cost);//, RHS2 = HermitianLB2(cost);
    //cout << RHS << " " << RHS2 << endl;
    //cout << HemitianLB(76) << " " << HermitianLB2(76) << endl;
    ws.isFSucc = (score > RHS);
    
#endif
    }
    else if (para.algoMode == "GSDSC")
    {

#ifdef GSDSC
    // GS DSC
    int numE = 0;
    for (int i = 0; i < codeLen ; i++)
        numE += (xSeq[i]!=ySeq[i]);

    
    ws.isFSucc = (numE < codeLen - sqrt((double)codeLen*infoLen));
#endif
    }
    else if (para.algoMode == "HermitianGSDSC")
    {
#ifdef HermitianGSDSC
    int numE = 0;
    for (int i = 0; i<codeLen; i++)
        numE += (xSeq[i]!=ySeq[i]);
    int g = (int)(sqrt((double)nbSize)*(sqrt((double)nbSize)-1) / 2);
    ws.isFSucc = (numE < codeLen - sqrt((double)codeLen * (infoLen + g - 1)));
#endif
    }
}

double DSCKVSim::HemitianLB(double cost)
{
    int g = static_cast<int>(sqrt((double)nbSize)*(sqrt((double)nbSize)-1) / 2 + 0.5);
    return (sqrt((double) 2 * (infoLen - 1 + g) * cost) + g);
}

double DSCKVSim::HermitianLB2(double cost)
{
    int g = static_cast<int>(sqrt((double)nbSize)*(sqrt((double)nbSize)-1) / 2 + 0.5);

    static vector<int> WG(2*g,0); /// Weierstrass Gap

    static bool isWGComputed = false;

    if (!isWGComputed)
    {
        int i = 0,j = 0;
        int q = (int)sqrt((double)nbSize);
        while (true)
        {
            while (true)
            {
                if ((q*i + (q+1)*j) < 2*g)
                    WG[q*i+(q+1)*j] = 1;
                else
                    break;
                j ++;
            }
            
            if (((++i)*q) >= 2*g)
                break;
        }
    }

    int w = 0;
    int sum = 0;
    while (true)
    {
        for (int j = 0; j< floor((double)w/(infoLen - 1 + g)); j++)
        {
            int index = w-(infoLen - 1 + g) * j;
            sum += index < 2*g ? WG[index] : 1;
        }

        if (sum > cost)
            return w;
        w ++ ;
    }
    return -1;
}
    




void DSCKVSim::oneThread(shared_ptr<RNObj> pobj, const MCThreadInfo& threadInfo)
{
    ThreadWorkSpace& ws = *(workspaceVec[threadInfo.threadID]);
    
    shared_ptr<XYSeq> pxyseq = static_pointer_cast<XYSeq, RNObj>(pobj);

    oneStep(pxyseq->xSeq, pxyseq->ySeq, ws);

}

void DSCKVSim::dumpOneResult(const MCThreadInfo& threadInfo)
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

void DSCKVSim::postProcess()
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

void DSCKVSim::prepare(const MCThreadInfo& threadInfo)
{
    for (int i = 0; i<threadInfo.nThread; i++)
    {
        workspaceVec.push_back(shared_ptr<ThreadWorkSpace>
            (new ThreadWorkSpace()));
    }
}

shared_ptr<RNObj> DSCKVSim::genRNObj()
{
    shared_ptr<XYSeq> pseq(new XYSeq(codeLen));
    for (int i = 0; i<codeLen; i++)
    {
        pseq->xSeq[i] = static_cast<int>(xDistr(uni01Gen));
        pseq->ySeq[i] = static_cast<int>(yxCondDistr[pseq->xSeq[i]](uni01Gen));
    }

    return pseq;
}



void DSCKVSim::DSCKVPara::loadPara(const paraMap& para, DispUI& ui)
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

    if (!getStrPara(para,"pdf para", pdfPara))
        paraMissingWarning(ui,"pdf para");

    if (!getStrPara(para,"Dec Pdf",decPdfPara))
        paraMissingWarning(ui,"Dec Pdf");
    if (decPdfPara == "NA")
        decPdfPara = "";

    if (!getStrPara(para,"algorithm mode", algoMode))
        paraMissingWarning(ui,"algorithm mode");

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

    if (!getIntPara(para, "RNG Table Option", rngTableOpt))
            paraMissingWarning(ui,"RNG Table Option");

}


void DSCKVSim::setCompResMap(compResMap& compResM_) const
{
    compRes.setMapping(compResM_);
}


void DSCKVSim::DSCKVCompRes::setMapping(compResMap &crMap) const
{
    Util::addCompRes(crMap, HX, "HX");
    Util::addCompRes(crMap, HY, "HY");
    Util::addCompRes(crMap, HYX, "HYX");
    Util::addCompRes(crMap, txRate, "Tx Rate");
    Util::addCompRes(crMap, startTime, "Start time");
    Util::addCompRes(crMap,infoLen, "info length");
    Util::addCompRes(crMap, gap,"gap");
}
    
bool DSCKVSim::checkTerm() 
{
    static vector<int> usedInfoLen;
    if (mulMode == 0) // not multi-run mode
    {
        usedInfoLen.clear();
        return (para.maxErrFr == simuRes.errFr);
    }
    else
    {

    // The multi-run mode is in an experimental phase. It runs multiple simulations with different parameters. The implementation
    // here is to find a code rate (info length) that gives FER < 10^-3. If one simulation with one parameter gives result
    // far from this number, it will adjust accordingly. But this requires careful writting of this part of the code and may still
    // need human monitoring, since the way the parameter is adjusted may not be stable (it may oscilate between two parameters, which are 
    // close to what you want.
        shared_ptr<DSCKVPara> npara(boost::static_pointer_cast<DSCKVPara, SimuPara>(nextPara));//nextPara is a copy of a DSCKVPara obj, so cast is safe
        bool termCurrStep = false;
        string message;

        if (simuRes.errFr == 10)
        {
            if (simuRes.FER == 1)
            {
                termCurrStep = true;
                npara->infoLen -= 10;
                message = "\nFER = 1. Quit\n";
            }
            else if (simuRes.FER >= 0.1)
            {
                termCurrStep = true;
                npara->infoLen -= 5;
                message = "\nFER >= 0.1. Quit\n";
            }
            else if (simuRes.FER >= 0.05)
            {
                termCurrStep = true;
                npara->infoLen -= 3;
                message = "\nFER >= 0.05. Quit\n";
            }
            else if (simuRes.FER >= 0.01 )
            {
                termCurrStep = true;
                npara->infoLen -= 2;
                message = "\nFER >= 0.01. Quit\n";
            }

            else if (simuRes.FER >= 0.005 )
            {
                termCurrStep = true;
                npara->infoLen -= 2;
                message = "\nFER >= 0.005. Quit\n";
            }

            else if (simuRes.FER >= 0.0015)
            {
                termCurrStep = true;
                npara->infoLen -= 1;
                message = "\nFER >= 0.0015. Quit\n";
            }

        }

        if ((simuRes.errFr == 50) && (simuRes.FER >= 0.0011))
        {


            termCurrStep = true;
            npara->infoLen -= 1;
            message = "\nFER >= 0.0011. Quit\n";
        }

        if ((simuRes.totFr == 20000) && (simuRes.FER <= 3e-4))
        {
            if (find(usedInfoLen.begin(), usedInfoLen.end(), infoLen) != usedInfoLen.end()) // was used 
            { 
                    termCurrStep = false;
            }
            else
            {
                
            termCurrStep = true;
            npara->infoLen += 1;
            message = "\nFER <= 3e-4. Quit\n";
            }
        }

        if ((simuRes.totFr == 20000) && (simuRes.FER <= 1e-8))
        {
            termCurrStep = true;
            npara->infoLen += 2;
            message = "\nFER == 0. Quit\n";
        }

        if (termCurrStep)
        {
            usedInfoLen.push_back(para.infoLen);
            dispUI.dispStr(message);
            resFile.saveStr(message);
            tempFile.saveStr(message);
            return true;
        }
        else
        {
            if (para.maxErrFr == simuRes.errFr)
            {
                isQuit = true;
                return true;
            }
            else
            {
                isQuit = false;
                return false;
            }
        }
    }




}



}
