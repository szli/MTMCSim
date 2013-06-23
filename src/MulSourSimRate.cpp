#include "MulSourSimRate.h"
#include "Util.h"

#define DISPLAY
namespace MTMCSim{

	MulSourSimRate::MulSourSimRate(MulSourParaRate& para_, TextSaveUI& resFile_, TextSaveUI& tempFile_, DispUI& dispUI_)
		: para(para_), resFile(resFile_), tempFile(tempFile_), dispUI(dispUI_),
		numE(para_.numE), codeLen(para_.codeLen), codePeriod(para_.codePeriod),
		pCross(para_.pCross), usedCodes(para_.usedCodes),
		uniGen(para_.seed), ldpcCodes(para_.ldpcFN, para_.allCodes.begin(), para_.allCodes.end()),
		maxErrFr(para_.maxErrFr),
		simuRes()
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

		// Initialize result object
		simuRes.aveIndRates.resize(numE,0);
		simuRes.totIndRates.resize(numE,0);
		

	}

	shared_ptr<RNObj> MulSourSimRate::genRNObj()
	{
		shared_ptr<ESeq> pseq(new ESeq(numE, codeLen));

		for (int i = 0; i<numE; i++)
			for (int j = 0; j<codeLen; j++)
				(pseq->eseq)[i][j] = BernolliGenVec[i](uniGen);

		return pseq;
	}


	void MulSourSimRate::prepare(const MCThreadInfo& threadInfo)
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


	void MulSourSimRate::oneThread(shared_ptr<RNObj> pobj, const MCThreadInfo& threadInfo)
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

	void MulSourSimRate::oneStep(const vector<vector<int> >& ESeqs, ThreadWorkSpace& ws)
	{
		//vector<vector<int> > ESeqs = ws.ESeqs;
		vector<vector<double> >& initLLRs = ws.initLLRs;
		vector<vector<int> >& decSeqs = ws.decSeqs;
		vector<vector<int> >& accSyns = ws.accSyns;
		vector<int>& isOneFrameSucc = ws.isOneFrameSucc;
		bool& isAllFrameSucc = ws.isAllFrameSucc;

		DLDPCEnDec& ldpcCodec = ws.ldpcCodec;


		for (int eIdx = 0; eIdx < numE; eIdx ++)
		{

			int tryCodeID = 0;
			for (tryCodeID = usedCodes[eIdx]; tryCodeID<codePeriod; tryCodeID++)
			{

				ldpcCodec.encodeOne(ESeqs[eIdx], tryCodeID, accSyns[eIdx]);
				isOneFrameSucc[eIdx] = ldpcCodec.decodeOne(accSyns[eIdx], tryCodeID, initLLRs[eIdx], decSeqs[eIdx]);

				if (isOneFrameSucc[eIdx])
					break;
			}// if not breaked, at the end tryCodeID = 66
			ws.actCodes[eIdx] = tryCodeID;
		}
	}

	void MulSourSimRate::dumpOneResult(const MCThreadInfo& threadInfo)
	{
		ThreadWorkSpace& ws = workspaceVec[threadInfo.threadID];
		// initialized in default constrcutor of Result object 
		simuRes.totFr ++;
		simuRes.totRateVec.push_back(0);

		simuRes.aveTotRate = 0;
		vector<double> rawIndRates(numE,0);
		for (int i = 0; i<numE; i++)
			{//init in constructor of the simulator class
				rawIndRates[i] = (double)ws.actCodes[i] / (double)codePeriod;
				simuRes.totIndRates[i] += rawIndRates[i]; // sum across samples
				simuRes.aveIndRates[i] = simuRes.totIndRates[i] / simuRes.totFr;
		}

		// Note that the actual sum rate computation is a little tricky. For a given realized sequence, it is possible that the sequence
		// is such that the relationship between k_i and k_i+1 fails, i.e., k_i+1 can > k_i. Thus, we need to be careful when determining the 
		//actual transmission rate for each source. Put in another word, for a given realization, it is possible that we cannot achieve the "best"
		double currSumRate = getSumRate(para.schMat,rawIndRates);
		simuRes.totTotRate += currSumRate;
		simuRes.aveTotRate = simuRes.totTotRate / simuRes.totFr;  // The average total rate up to now
		simuRes.totRateVec[simuRes.totFr-1] = currSumRate;// store the the sum rate of current instance 

		double sum = 0;
		for (int i = 0; i< simuRes.totFr;i++)
		{
			sum += (simuRes.totRateVec[i] - simuRes.aveTotRate) * (simuRes.totRateVec[i] - simuRes.aveTotRate);
		}
		simuRes.stdev = simuRes.totFr == 1? 0 : sqrt(1.0/(simuRes.totFr-1) * sum);
		
#ifdef DISPLAY
		cout << "Required rates (in terms of codeID):" ;

		for (int i= 0; i<numE; i++)
		{
			cout << ws.actCodes[i] << " " ;
		}
		cout << endl;


		bool isIncr = true;
		for (int i = 0; i<numE-1; i++)
			if (ws.actCodes[i] > ws.actCodes[i+1])
				isIncr = false;

		cout << "is Incr = " << isIncr << endl;

#endif


		// one step output
		if (simuRes.totFr % para.saveTempPerFr == 0)
		{
			ostringstream ress;
			ress << "Frame: " << simuRes.totFr << "Ave Rates";

			for (int i = 0; i<numE; i++)
				ress << simuRes.aveIndRates[i] << " ";

			ress  << "ave Tot Rate " << simuRes.aveTotRate << "stdev" << simuRes.stdev << "totRate " << simuRes.totRateVec[simuRes.totFr-1] << "\n";

			tempFile.saveStr(ress.str());
			resFile.saveStr(ress.str());
		}

		if (simuRes.totFr % para.dispPerFr == 0)
		{
			ostringstream ress;
			ress << "Fr: " << simuRes.totFr << "ARates";

			for (int i = 0; i<numE; i++)
				ress << simuRes.aveIndRates[i] << " ";
			ress  << "aTRate " << simuRes.aveTotRate << "stdev" << simuRes.stdev << "tRate " << simuRes.totRateVec[simuRes.totFr-1] <<  "\r";

		/*	for (int i = 0; i<numE; i++)
				ress << simuRes.indiviErrFr[i] << " ";
			
			ress <<
				" FER: " << simuRes.FER << " BER:" << simuRes.BER << "\r";*/
			dispUI.dispStr(ress.str());
		}
	}

	void MulSourSimRate::postProcess()
	{
		ostringstream ress;
		ress << "Frame: " << simuRes.totFr << "Ave Indi Rates";

			for (int i = 0; i<numE; i++)
				ress << simuRes.aveIndRates[i] << " ";

			ress  << "ave Tot Rate " << simuRes.aveTotRate << "stdev" << simuRes.stdev << "\n";
		resFile.saveStr(ress.str());

		//ress.str(string());
	//	ress << "\nSimulation terminated.\n";
	//	ress <<  "Frame: " << simuRes.totFr << " Err Fr: " << simuRes.errFr << 
			//" FER: " << simuRes.FER << " BER:" << simuRes.BER <<  endl;
		//ress.flush();
		dispUI.dispStr(ress.str());
	}	

	void MulSourSimRate::MulSourParaRate::loadPara(const paraMap &para, MTMCSim::DispUI &ui)
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


		string sschMat;
		getStrPara(para,"Schedule Matrix", sschMat);
		int idx = 0;
		stringstream ss;
		vector<int> oneR;
		while (idx <= sschMat.size()-1)
		{
			if (sschMat[idx]!=';')
			{
				ss << sschMat[idx];
				idx ++;
			
			}
			else
			{
				
				int value;
				while (ss >> value)
					oneR.push_back(value);
				schMat.push_back(oneR);
				idx ++;
				ss.str(string());
				ss.clear();
				oneR.clear();
			}
		}// end while

	for (int i = 2; i<= codePeriod; i++)
		allCodes.push_back(i);

					
				

	}

	void MulSourSimRate::MulSourCompRes::setMapping(compResMap &crMap) const
	{
		Util::addCompRes(crMap,pCrossStr,"pCross");
		Util::addCompRes(crMap,codeIDStr, "code IDs");
		Util::addCompRes(crMap,totEnt,"Total Entro");
		Util::addCompRes(crMap,EntStr,"Entro");
		Util::addCompRes(crMap, totRate,  "Total Tx Rate");
		Util::addCompRes(crMap, rateStr, "Tx Rate");
		Util::addCompRes(crMap,startTime,  "Start time");
	}

	bool MulSourSimRate::checkTerm() 
	{ return (para.maxErrFr == simuRes.totFr);}


	void MulSourSimRate::OneKnownFrame(const string& seqFN)
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

	
	double MulSourSimRate::getSumRate(const vector<vector<int> >& schMat, const vector<double>& indRates)
	{
		int nSour = schMat[0].size();
		vector<double> actIndRates(nSour,0);
		for (int i = 0; i<numE; i++)
		{
			for (int j = 0; j<nSour; j++)
				if (schMat[i][j] == 1 && actIndRates[j]<indRates[i])
					actIndRates[j] = indRates[i];
		}



		double sum = 0;
		for (int i = 0; i<nSour; i++)
			sum += actIndRates[i];
		return sum;
	}


	

	








































		

		


}

