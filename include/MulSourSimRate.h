#ifndef _MULSOURSIMRATE_H_INCLUDED
#define _MULSOURSIMRATE_H_INCLUDED

#include "MTMCSimBase.h"
#include "RNObj.h"
#include "StdHeaders.h"
#include "SimuPara.h"
#include "CompRes.h"
#include "DLDPCCodes.h"
#include "DLDPCEnDec.h"
#include "TextSaveUI.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/bernoulli_distribution.hpp>

namespace MTMCSim{
/** Multiple source linear correlation model simulation.
*	Simulation only simulates the encoding and decoding of the "error" vectors. 
*/
class MulSourSimRate:public MTMCSimBase{

public:

	/** "Error" sequences.
	*/
	class ESeq:public RNObj{
	public:
		vector<vector<int> > eseq;
		ESeq(int numE, int len):eseq(numE, vector<int>(len))
		{}
	};

	/** The parameters for the simulation. All members are public. Most members are the parameters. */
	class MulSourParaRate:public SimuPara{
	public:
		int numE;
		int codeLen;
		int codePeriod;

		int seed;

		int maxErrFr, dispPerFr, saveTempPerFr;

		vector<double> pCross;

		vector<int> usedCodes;

		string ldpcFN /*, saveFN*/ ;

		vector<vector<int> > schMat;

		vector<int> allCodes;

		void loadPara(const paraMap& para, DispUI& ui);

		shared_ptr<SimuPara> clone() const
		{ return shared_ptr<MulSourParaRate>(new MulSourParaRate(*this)); }
	};

	/** The computed results from input parameters. This does not include the MC simulation results. */
	class MulSourCompRes:public CompRes{
	public:
		string pCrossStr;
		string codeIDStr;
		double totEnt;
		string EntStr;
		double totRate;
		string rateStr;
		string startTime;
		vector<int> nTxBits;
		void setMapping(compResMap &crMap) const;
		MulSourCompRes():totEnt(0), totRate(0)
			{}

	};
	
	/** Work space for each thread. */
	class ThreadWorkSpace
	{
	public:

		vector<vector<int> > decSeqs; ///< Decoded sequences.
		vector<vector<double> > initLLRs; ///< Initial LLRs.
		vector<vector<int> > accSyns;
		vector<int> isOneFrameSucc;
		bool isAllFrameSucc;
		int nBitErr;
		vector<int> actCodes;


		DLDPCEnDec ldpcCodec;
		ThreadWorkSpace(int codeLen, int numE, const vector<int>& nTxBits, const DLDPCCodes& ldpcCodes)
			: decSeqs(numE, vector<int>(codeLen)),
			initLLRs(numE, vector<double>(codeLen)), ldpcCodec(ldpcCodes), accSyns(numE),
			isOneFrameSucc(numE),isAllFrameSucc(false),actCodes(numE,0)
		{
			for (int i = 0; i<numE; i++)
				accSyns[i].resize(codeLen);
		}

	};
	/** The simulation results. */
	class MulSourSimuRes{
	public:
		int totFr;
		vector<double> aveIndRates;
		vector<double> totIndRates; // sum across the simulation samples
		vector<double> totRateVec; // the sum rates over the error sources
		double aveTotRate;// sum accross the error sources
		double stdev;
		double totTotRate; // sum of the sum rate 

		MulSourSimuRes():totFr(0),aveIndRates(0),totIndRates(0),aveTotRate(0),totTotRate(0)
		{
		
		}

	};






	MulSourSimRate(MulSourParaRate& para_, TextSaveUI& resFile_, TextSaveUI& tempFile_, DispUI& dispUI_);




	virtual void oneThread(shared_ptr<RNObj>, const MCThreadInfo&);

	virtual shared_ptr<RNObj> genRNObj();

	virtual void dumpOneResult(const MCThreadInfo&);
	
	virtual bool checkTerm();

	void prepare(const MCThreadInfo& threadInfo);

	void postProcess();

	/** Calls compRes.setMapping(compResM_);.*/
	void setCompResMap(compResMap& compResM_) const
	{
		compRes.setMapping(compResM_);
	}
	
	void OneKnownFrame(const string&);

private:

	void oneStep(const vector<vector<int> >& ESeqs, ThreadWorkSpace& ws);

	double getSumRate(const vector<vector<int> >& schMat, const vector<double>& indRates);

	//internal parameters
	int numE;
	int codeLen;
	int codePeriod;
	vector<double> pCross;
	vector<int> usedCodes;
	int maxErrFr;

	//functional components
	DLDPCCodes ldpcCodes;
	vector<ThreadWorkSpace> workspaceVec;
	//rng
	boost::mt19937 uniGen;
	vector<boost::bernoulli_distribution<> > BernolliGenVec;







	//parameter in and out
	MulSourParaRate para;
	MulSourSimuRes simuRes;
	MulSourCompRes compRes;
	TextSaveUI& resFile;
	TextSaveUI& tempFile;
	DispUI& dispUI;



};



}
#endif
