#include "MulStageSim.h"

namespace MTMCSim{

extern const string default_pdfdir; // defined in DSCKVSim.cpp

MulStageSim::MulStageSim(MulStagePara& para_, TextSaveUI& resFile_, TextSaveUI& tempFile_, DispUI& dispUI_)
: nbSize(para_.nbSize), nbCodeLen(para_.nbCodeLen), codePer(para_.codePer),
maxErrFr(para_.maxErrFr), 
ldpcCodes(para_.ldpcFN, para_.usedCodesRaw.begin(), para_.usedCodesRaw.end()),
uniGen(para_.seed), uni01Gen(uniGen, uni_01), para(para_), 
resFile(resFile_), tempFile(tempFile_), dispUI(dispUI_)
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
	
	nbBinConv.resize(nbSize);
	if (!para_.decPdfPara.empty())
	{
		pdfXY decPdf;

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

		for (int i = 0; i<nbSize; i++)
			nbBinConv[i].setPdf(decPdf.getYXCondPdf()[i], nbSize);

	}
	else
	{
		for (int i = 0; i<nbSize; i++)
				nbBinConv[i].setPdf(nbJoPdf.getYXCondPdf()[i], nbSize);
		}



	usedCodesRaw = para_.usedCodesRaw;

	double ln2 = log(2.0);
	nBits = static_cast<int>(log((double)nbSize) / ln2 + 0.5); // Without + 0.5, may go wrong under Linux (size = 4096) because the division gives a value close to but smaller than 12

	compRes.condEntro.resize(nBits, 0);

	const vector<double>& xMarPdf = nbJoPdf.getXMarPdf();

	for (int i = 0; i<nBits; i++)
		for (int j = 0; j<nbSize; j++)
		{
			compRes.condEntro[i] +=  xMarPdf[j] * nbBinConv[j].getOneCondEn(i);
		}

	compRes.HX = nbJoPdf.getHX();
	compRes.HY = nbJoPdf.getHY();
	compRes.HYX = nbJoPdf.getCondHYX();

	compRes.nTxBits.resize(nBits);
	compRes.codeRates.resize(nBits);
	int totTxBits = 0;

	for (int i = 0; i<nBits; i++)
	{
		compRes.nTxBits[i] = usedCodesRaw[i] * (nbCodeLen/codePer);
		totTxBits += compRes.nTxBits[i];
		compRes.codeRates[i] = 1 - (double)compRes.nTxBits[i] / nbCodeLen;
	}

	compRes.totTxRate = (double)totTxBits/nbCodeLen;

	//compRes.codeIDStr is empty by default construction
	for (int i=0; i<nBits; i++)
		compRes.codeIDstr += Util::intToStr(usedCodesRaw[i]) + " ";

	time_t timep;
	time (&timep);

	#define _CRT_SECURE_NO_WARNINGS
	compRes.startTime = asctime(localtime(&timep));

	//ldpcCodes.setCodes(para_.ldpcFN, usedCodesRaw.begin(), usedCodesRaw.end());

	//Init the RNG
	initRNG();
	/*xDistr.setDistribution(nbSize, nbJoPdf.getXMarPdf());

	yxCondDistr.resize(nbSize);

	for ( int i = 0; i<nbSize; i++)
		yxCondDistr[i].setDistribution(nbSize, nbJoPdf.getYXCondPdf()[i]);*/

	//Init simulation result
	simuRes.totFr = 0;
	simuRes.errFr = 0;
	simuRes.FER = 0;
	
	simuRes.binSourBER.resize(nBits,0);
	simuRes.binErrBit.resize(nBits,0);
	simuRes.binErrFrame.resize(nBits,0);
	simuRes.binSourFER.resize(nBits,0);





}


void MulStageSim::initRNG()
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


void MulStageSim::oneThread(shared_ptr<RNObj> pobj, const MCThreadInfo& threadInfo)
{
	ThreadWorkSpace& ws = *(workspaceVec[threadInfo.threadID]);
	DLDPCEnDec& codec = ws.ldpcCodec;

	//ws.initLLR.resize(nbCodeLen);
	//ws.binSourSeqs.resize(nBits, vector<int>(nbCodeLen));
	//ws.decBinSeqs.resize(nBits, vector<int>(nbCodeLen));
	//
	//ws.encSyns.resize(nBits);
	//for (int i = 0; i< nBits; i++)
	//	ws.encSyns[i].resize(compRes.nTxBits[i]);

	//ws.isBSucc.resize(nBits);
	//ws.nBErr.resize(nBits);

	//DLDPCEnDec codec(ldpcCodes);

	// load the random sequences

	shared_ptr<XYSeq> pxyseq = static_pointer_cast<XYSeq, RNObj>(pobj);

	oneStep(pxyseq->xSeq, pxyseq->ySeq, ws, codec);

	ws.isFSucc = true;

	for (int i = 0; i<nBits; i++)
		if (ws.isBSucc[i] == false)
			ws.isFSucc = false;


}

void MulStageSim::postProcess()
{
		ostringstream ress;
		ress << "Simulation terminated." << endl;
		ress << "Frame: " << simuRes.totFr << " Err Fr: " << simuRes.errFr << 
			" FER: " << simuRes.FER << "\n";
		std::ostream_iterator<double> out_iter(ress, " ");
		std::copy(simuRes.binSourFER.begin(), simuRes.binSourFER.end(), out_iter);
		ress << "\n";
		ress.flush();
		resFile.saveStr(ress.str());

		ress.str(string());
		ress << "\nSimulation terminated.\n";
		ress <<  "Frame: " << simuRes.totFr << " Err Fr: " << simuRes.errFr << 
			" FER: " << simuRes.FER << endl;

		std::ostream_iterator<double> out_iter2(ress, " ");
		std::copy(simuRes.binSourFER.begin(), simuRes.binSourFER.end(), out_iter2);
		ress << endl;
		ress.flush();
		dispUI.dispStr(ress.str());
}



void MulStageSim::prepare(const MCThreadInfo& threadInfo)
{
	for (int i = 0; i<threadInfo.nThread; i++)
	{
		workspaceVec.push_back(shared_ptr<ThreadWorkSpace>
			(new ThreadWorkSpace(nbCodeLen, nBits, compRes.nTxBits, ldpcCodes)));
	}
}


shared_ptr<RNObj> MulStageSim::genRNObj()
{
	shared_ptr<XYSeq> pseq(new XYSeq(nbCodeLen));
	for (int i = 0; i<nbCodeLen; i++)
	{
		pseq->xSeq[i] = static_cast<int>(xDistr(uni01Gen));
		pseq->ySeq[i] = static_cast<int>(yxCondDistr[pseq->xSeq[i]](uni01Gen));
	}

	return pseq;
}

	
void MulStageSim::oneStep(const std::vector<int> &xSeq, const std::vector<int> &ySeq, ThreadWorkSpace& ws, DLDPCEnDec& codec)
{
	vector<double>& initLLR = ws.initLLR;
	vector<vector<int> >& binSourSeqs = ws.binSourSeqs;
	vector<vector<int> >& decBinSeqs = ws.decBinSeqs;
	vector<vector<int> >& encSyns = ws.encSyns;
	vector<int>& isBSucc = ws.isBSucc, & nBErr = ws.nBErr;

	// convert YSeq to binary sequences	
	// use binary expression of the decimal number
	for (int i = 0; i<nbCodeLen; i++)
	{
		int YSym = ySeq[i];
		for (int j = 0; j<nBits; j++)
		{
			binSourSeqs[j][i] = YSym % 2;
			YSym = YSym >> 1;
		}
	}

	for (int i = 0; i<nBits; i++)
		codec.encodeOne(binSourSeqs[i], usedCodesRaw[i], encSyns[i]);

	double ln2 = log(2.0);
	vector<int> givB(nBits-1);

	// Decode
	for (int i = 0; i<nBits; i++)
	{
		if (usedCodesRaw[i] == codePer) // decode by matrix inversion, actuallly no decoding is performed
		{
			for (int j = 0; j< nbCodeLen; j++)
				decBinSeqs[i][j] = binSourSeqs[i][j];
			nBErr[i] = 0;
			isBSucc[i] = 1;
			continue;
		}


		double zeroP;
		if (i == 0) // b0
		{
			for (int j = 0; j<nbCodeLen; j++)
			{
				zeroP = nbBinConv[xSeq[j]].getZeroPGiven(0, 0);
				initLLR[j] = log(zeroP/(1-zeroP)) / ln2;

			}
		}

		else 
		{
			for (int j = 0; j<nbCodeLen; j++)
			{
				
				for (int k = 0; k<i; k++)
					givB[k] = decBinSeqs[k][j];


				zeroP = nbBinConv[xSeq[j]].getZeroPGiven(&givB[0], i);
				initLLR[j] = log(zeroP/(1-zeroP)) / ln2;
			}
		}



		isBSucc[i] = codec.decodeOne(encSyns[i], usedCodesRaw[i], initLLR,decBinSeqs[i]);


		nBErr[i] = 0;
		//vector<int> errPos;// record the error bit position
		for (int k = 0; k<nbCodeLen; k++)
		{
			nBErr[i] += (decBinSeqs[i][k] != binSourSeqs[i][k]);
		}
	}// end of decode
}// end of oneStep()

void MulStageSim::dumpOneResult(const MCThreadInfo& threadInfo)
{
	ThreadWorkSpace& ws = *(workspaceVec[threadInfo.threadID]);

		// update stat
		simuRes.totFr ++;
		
		if (!ws.isFSucc)
			simuRes.errFr ++;

		for (int i = 0; i<nBits;i++)
			simuRes.binErrBit[i] += ws.nBErr[i];

		for (int i = 0; i<nBits; i++)
		{
			if (!ws.isBSucc[i])
				simuRes.binErrFrame[i] ++;
			simuRes.binSourFER[i] = (double)simuRes.binErrFrame[i]/(double)simuRes.totFr;
		}
		
		// updata result

		simuRes.FER = (double)simuRes.errFr / (double)simuRes.totFr;

		for (int i = 0; i<nBits;i++)
			simuRes.binSourBER[i] = (double)simuRes.binErrBit[i] / (double) (simuRes.totFr * nbCodeLen);

		// one step output
		if (simuRes.totFr % para.saveTempPerFr == 0)
		{
			ostringstream ress;
			ress << "Frame: " << simuRes.totFr << " Err Fr: " << simuRes.errFr << 
				" FER: " << simuRes.FER << "\n";
			std::ostream_iterator<double> out_iter(ress, " ");
			std::copy(simuRes.binSourFER.begin(), simuRes.binSourFER.end(), out_iter);
			ress << "\n";
			tempFile.saveStr(ress.str());
			resFile.saveStr(ress.str());
		}

		if (simuRes.totFr % para.dispPerFr == 0)
		{
			ostringstream ress;
			ress <<  "Frame: " << simuRes.totFr << " Err Fr: " << simuRes.errFr << 
				" FER: " << simuRes.FER << "\r";
			dispUI.dispStr(ress.str());
		}



	//Termination check is moved to the scheduler
}

void MulStageSim::MulStagePara::loadPara(const paraMap& para, DispUI& ui)
{
	using namespace Util;
	if (!getIntPara(para, "alphabet size", nbSize))
		paraMissingWarning(ui, "alphabet size");

	if (!getIntPara(para, "code length", nbCodeLen))
		paraMissingWarning(ui, "code length");

	if (!getIntPara(para, "code period", codePer))
		paraMissingWarning(ui, "code period");


	nBits = static_cast<int>(log((double)nbSize) / log(2.0) + 0.5); 

	usedCodesRaw.resize(nBits);
	for (int i = 0; i < nBits; i++)
	{
		string key = "code ID " + intToStr(i);
		if (!getIntPara(para,key,usedCodesRaw[i]))
			paraMissingWarning(ui,key);
	}

	if (!getStrPara(para,"pdf para", pdfPara))
		paraMissingWarning(ui,"pdf para");

	if (!getStrPara(para,"Dec Pdf",decPdfPara))
		paraMissingWarning(ui,"Dec Pdf");
	if (decPdfPara == "NA")
		decPdfPara = "";


	if (!getIntPara(para, "seed", seed))
		paraMissingWarning(ui, "seed");

	if (!getStrPara(para,"ldpc file", ldpcFN))
		paraMissingWarning(ui,"ldpc file");

	if (!getIntPara(para, "max error frame", maxErrFr))
		paraMissingWarning(ui, "max error frame");

	if (!getIntPara(para, "display freq", dispPerFr))
		paraMissingWarning(ui, "display freq");

	if (!getIntPara(para, "save freq", saveTempPerFr))
		paraMissingWarning(ui, "save freq");

	if (!getStrPara(para,"save file", saveFN))
		paraMissingWarning(ui,"save file");

	if (!getIntPara(para, "RNG Table Option", rngTableOpt))
			paraMissingWarning(ui,"RNG Table Option");

}




void MulStageSim::oneKnownFrame()// Testing only
{
	MCThreadInfo threadInfo(1,0);

	prepare(threadInfo);

	shared_ptr<XYSeq> pxySeq(new XYSeq(nbCodeLen));


	ifstream seqFile("tempSeq.txt");

	for (int i = 0; i<nbCodeLen; i++)
		seqFile >> pxySeq->xSeq[i];
	for (int i = 0; i<nbCodeLen; i++)
		seqFile >> pxySeq->ySeq[i];

	
	oneThread(pxySeq, threadInfo);
}

void MulStageSim::setCompResMap(compResMap& compResM_) const
{
	compRes.setMapping(compResM_);
}



MulStageSim::ThreadWorkSpace::ThreadWorkSpace(int nbCodeLen, int nBits, const  vector<int>& nTxBits, const DLDPCCodes& ldpcCodes)
: ldpcCodec(ldpcCodes), initLLR(nbCodeLen), binSourSeqs(nBits, vector<int>(nbCodeLen)),
  decBinSeqs(nBits, vector<int>(nbCodeLen)), encSyns(nBits)
{
	//initLLR.resize(nbCodeLen);
	//binSourSeqs.resize(nBits, vector<int>(nbCodeLen));
	//decBinSeqs.resize(nBits, vector<int>(nbCodeLen));
	//
	//encSyns.resize(nBits);
	for (int i = 0; i< nBits; i++)
		encSyns[i].resize(nTxBits[i]);

	isBSucc.resize(nBits);
	nBErr.resize(nBits);

}

void MulStageSim::MulStageCompRes::setMapping(compResMap &crMap) const
{

	std::ostringstream ostr;
	std::ostream_iterator<double> osit(ostr," ");
	std::copy(condEntro.begin(), condEntro.end(), osit);
	ostr.flush();
	Util::addCompRes(crMap, ostr.str(),"CondEntro");
	ostr.str(string());

	Util::addCompRes(crMap, codeIDstr,  "code IDs");

	Util::addCompRes(crMap, HX, "HX");
	Util::addCompRes(crMap, HY, "HY");
	Util::addCompRes(crMap, HYX, "HYX");

	for (unsigned int i = 0; i<codeRates.size(); i++)
	{
		ostr << 1 - codeRates[i] << " ";
	}
	ostr.flush();
	Util::addCompRes(crMap, ostr.str(), "Tx Rate");

	Util::addCompRes(crMap, totTxRate, "Total Tx Rate");

	Util::addCompRes(crMap, startTime, "Start time");
}
	
	



bool MulStageSim::checkTerm() 
{
	return (para.maxErrFr == simuRes.errFr);
}




}//end of namespace
