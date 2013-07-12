#include "MulSourSim.h"
#include "Util.h"

namespace MTMCSim{

	MulSourSim::MulSourSim(MulSourPara& para_, TextSaveUI& resFile_, TextSaveUI& tempFile_, DispUI& dispUI_)
		: para(para_), resFile(resFile_), tempFile(tempFile_), dispUI(dispUI_),
		numE(para_.numE), codeLen(para_.codeLen), codePeriod(para_.codePeriod),
		pCross(para_.pCross), usedCodes(para_.usedCodes),
		uniGen(para_.seed), ldpcCodes(para_.ldpcFN, para_.usedCodes.begin(), para_.usedCodes.end()),
		maxErrFr(para_.maxErrFr),
		simuRes(para_.numE)
	{

		// compRes is default constructed and everything is clear.
		for (int i = 0; i<numE; i++)
		{
			compRes.pCrossStr += Util::douToStr(pCross[i]) + " " ;
			compRes.codeIDStr += Util::intToStr(usedCodes[i]) + " ";
			double ent = Util::bEntro(pCross[i]);
			compRes.EntStr += Util::douToStr(ent);
			compRes.totEnt += ent;
			compRes.rateStr += Util::douToStr((double)usedCodes[i] / codePeriod) + " ";
			compRes.totRate += (double)usedCodes[i] / codePeriod;
			compRes.nTxBits.push_back(usedCodes[i] * (codeLen/codePeriod));
		}

		time_t timep;
		time (&timep);

		#define _CRT_SECURE_NO_WARNINGS
		compRes.startTime = asctime(localtime(&timep));

		for (int i = 0; i<numE; i++)
			BernolliGenVec.push_back(boost::bernoulli_distribution<>(pCross[i]));


	}

	shared_ptr<RNObj> MulSourSim::genRNObj()
	{
		shared_ptr<ESeq> pseq(new ESeq(numE, codeLen));

		for (int i = 0; i<numE; i++)
			for (int j = 0; j<codeLen; j++)
				(pseq->eseq)[i][j] = BernolliGenVec[i](uniGen);

		return pseq;
	}


	void MulSourSim::prepare(const MCThreadInfo& threadInfo)
	{

		double ln2 = log((double)2);
		for (int i = 0; i<threadInfo.nThread; i++)
		{
			workspaceVec.push_back(ThreadWorkSpace(codeLen, numE, compRes.nTxBits, ldpcCodes));
			for (int j = 0; j<numE; j++)
				for (int k = 0; k < codeLen; k++)
					workspaceVec[i].initLLRs[j][k] = log((1-pCross[j])/pCross[j]) / ln2;
		}		
	}


	void MulSourSim::oneThread(shared_ptr<RNObj> pobj, const MCThreadInfo& threadInfo)
	{
		ThreadWorkSpace& ws = workspaceVec[threadInfo.threadID];
		DLDPCEnDec& codec = ws.ldpcCodec;

		shared_ptr<ESeq> peseq = boost::static_pointer_cast<ESeq,RNObj>(pobj);

		oneStep(peseq->eseq, ws);

		ws.isAllFrameSucc = true;
		
		for (int i = 0; i<numE; i++)
			if (ws.isOneFrameSucc[i] == false)
				ws.isAllFrameSucc = false;
	}

	void MulSourSim::oneStep(const vector<vector<int> >& ESeqs, ThreadWorkSpace& ws)
	{
		//vector<vector<int> > ESeqs = ws.ESeqs;
		vector<vector<double> >& initLLRs = ws.initLLRs;
		vector<vector<int> >& decSeqs = ws.decSeqs;
		vector<vector<int> >& accSyns = ws.accSyns;
		vector<int>& isOneFrameSucc = ws.isOneFrameSucc;
		bool& isAllFrameSucc = ws.isAllFrameSucc;

		DLDPCEnDec& ldpcCodec = ws.ldpcCodec;


		for (int i = 0; i<numE; i++)
			ldpcCodec.encodeOne(ESeqs[i], usedCodes[i], accSyns[i]);

		for (int i = 0; i<numE; i++)
			isOneFrameSucc[i] = ldpcCodec.decodeOne(accSyns[i], usedCodes[i], initLLRs[i], decSeqs[i]);

		ws.nBitErr = 0;
		bool thisBitErr = false;
		for (int i = 0; i<codeLen; i++)
		{
			for (int j = 0; j<numE; j++)
			{
				if (ESeqs[j][i] != decSeqs[j][i])
				{
					thisBitErr = true;
					break;
				}
			}

			ws.nBitErr += thisBitErr;
		}
	}

	void MulSourSim::dumpOneResult(const MCThreadInfo& threadInfo)
	{
		ThreadWorkSpace& ws = workspaceVec[threadInfo.threadID];
		// initialized in default constrcutor 
		simuRes.totFr ++;

		if (!ws.isAllFrameSucc)
		{
			simuRes.errFr ++;
			simuRes.errBit += ws.nBitErr;
		}

		for (int i = 0; i< numE; i++)
			if (!ws.isOneFrameSucc[i])
				simuRes.indiviErrFr[i] ++;


		simuRes.FER = (double)simuRes.errFr / (double)simuRes.totFr;
		simuRes.BER = (double)simuRes.errBit / (double)(simuRes.totFr * codeLen);

		for (int i = 0; i<numE; i++)
			simuRes.indiviSourFER[i] = (double)simuRes.indiviErrFr[i] / (double)simuRes.totFr;

		// one step output
		if (simuRes.totFr % para.saveTempPerFr == 0)
		{
			ostringstream ress;
			ress << "Frame: " << simuRes.totFr << " Err Fr: " << simuRes.errFr << " IndiFER:";

			for (int i = 0; i<numE; i++)
				ress << simuRes.indiviSourFER[i] << " ";

			
			ress  << " FER: " << simuRes.FER << " BER:" << simuRes.BER << "\n" << "IndiErrFr:";
			std::ostream_iterator<int> out_iter(ress, " ");
			std::copy(simuRes.indiviErrFr.begin(), simuRes.indiviErrFr.end(), out_iter);
			ress << endl;
			tempFile.saveStr(ress.str());
			resFile.saveStr(ress.str());
		}

		if (simuRes.totFr % para.dispPerFr == 0)
		{
			ostringstream ress;
			ress <<  "Frame: " << simuRes.totFr << " Err Fr: " << simuRes.errFr <<" FER: " << simuRes.FER <<  "\r";

		/*	for (int i = 0; i<numE; i++)
				ress << simuRes.indiviErrFr[i] << " ";
			
			ress <<
				" FER: " << simuRes.FER << " BER:" << simuRes.BER << "\r";*/
			dispUI.dispStr(ress.str());
		}
	}

	void MulSourSim::postProcess()
	{
		ostringstream ress;
		ress << "Simulation terminated." << endl;
		ress << "Frame: " << simuRes.totFr << " Err Fr: " << simuRes.errFr << 
			" FER: " << simuRes.FER << " BER:" << simuRes.BER << "\n" <<"IndiFER:";
		std::ostream_iterator<double> out_iter(ress, " ");
		std::copy(simuRes.indiviSourFER.begin(), simuRes.indiviSourFER.end(), out_iter);
		ress << endl;
		resFile.saveStr(ress.str());

		ress.str(string());
		ress << "\nSimulation terminated.\n";
		ress <<  "Frame: " << simuRes.totFr << " Err Fr: " << simuRes.errFr << 
			" FER: " << simuRes.FER << " BER:" << simuRes.BER <<  endl;
		ress.flush();
		dispUI.dispStr(ress.str());
	}	

	void MulSourSim::MulSourPara::loadPara(const paraMap &para, MTMCSim::DispUI &ui)
	{
		using namespace Util;

		if (!getIntPara(para,"seed", seed))
			paraMissingWarning(ui,"seed");

		if (!getIntPara(para,"num error", numE))
			paraMissingWarning(ui,"num error");

		if (!getIntPara(para,"code length", codeLen))
			paraMissingWarning(ui,"code length");

		if (!getIntPara(para,"code period", codePeriod))
			paraMissingWarning(ui,"code period");

		pCross.resize(numE);
		usedCodes.resize(numE);

		string key;

		for (int i = 0; i< numE; i++)
		{
			key = "pCross " + intToStr(i);

			if (!getDouPara(para,key,pCross[i]))
				paraMissingWarning(ui,key);

			key = "code ID " + intToStr(i);
			
			if (!getIntPara(para,key,usedCodes[i]))
				paraMissingWarning(ui,key);
		}

		if (!getStrPara(para,"ldpc file",ldpcFN))
			paraMissingWarning(ui,"ldpc file");

		if (!getStrPara(para,"save file",saveFN))
			paraMissingWarning(ui,"save file");

		if (!getIntPara(para,"max error frame", maxErrFr))
			paraMissingWarning(ui,"max error frame");

		if (!getIntPara(para,"display freq", dispPerFr))
			paraMissingWarning(ui,"display freq");

		if (!getIntPara(para,"save freq", saveTempPerFr ))
			paraMissingWarning(ui,"save freq");

	}

	void MulSourSim::MulSourCompRes::setMapping(compResMap &crMap) const
	{
		Util::addCompRes(crMap,pCrossStr,"pCross");
		Util::addCompRes(crMap,codeIDStr, "code IDs");
		Util::addCompRes(crMap,totEnt,"Total Entro");
		Util::addCompRes(crMap,EntStr,"Entro");
		Util::addCompRes(crMap, totRate,  "Total Tx Rate");
		Util::addCompRes(crMap, rateStr, "Tx Rate");
		Util::addCompRes(crMap,startTime,  "Start time");
	}

	bool MulSourSim::checkTerm() 
	{ return (para.maxErrFr == simuRes.errFr);}


	void MulSourSim::OneKnownFrame(const string& seqFN)
	{
		ifstream seqFile(seqFN.c_str());

		shared_ptr<ESeq> peseq(new ESeq(numE, codeLen));

		for (int i = 0; i< numE; i++)
			for (int j = 0; j<codeLen; j++)
				if (!(seqFile >> (peseq->eseq)[i][j]))
					dispUI.dispStr("Input File Error!");
		MCThreadInfo threadInfo(1,0);
		oneThread(genRNObj(), threadInfo);
		dumpOneResult(threadInfo);

		seqFile.close();
	}






	











































		

		


}

