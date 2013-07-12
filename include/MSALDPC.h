#ifndef _MSALDPC_H_INCLUDED
#define _MSALDPC_H_INCLUDED
#include "StdHeaders.h"
#include "RNObj.h"
#include <boost/shared_ptr.hpp>
#include "pdfXY.h"
//#include "DLDPCCodes.h"
#include "discrete_distribution.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>
#include "nbToBinPDFConv.h"
#include <boost/pointer_cast.hpp>
//#include "DLDPCEnDec.h"
#include "Util.h"
#include "TextSaveUI.h"
#include "DOSDispUI.h"
#include "MCThreadInfo.h"
#include "CompRes.h"
#include "SimuPara.h"
#include "MTMCSimBase.h"
#include "LDPCDec.h"



namespace MTMCSim{
using boost::shared_ptr;
using boost::static_pointer_cast;
/** Simulator for multistage LDPC for two source distributed source coding. */
class MSALSim: public MTMCSimBase{
public:

	//typedef vector<shared_ptr<LDPCDec> > LDPCGroup;
	/** A random number object for this problem is two correlated sequence X and Y. Note X is the side 
	*	information and Y is the source to be transmitted.
	*/
	class XYSeq: public RNObj
	{
		public:
		vector<int> xSeq;// X is side info
		vector<int> ySeq;// encode Y
		XYSeq(int len):xSeq(len), ySeq(len)
		{ }
	};

	class ALDPCInfo{
	public:

		struct OneCodeInfo{
				int nCol, nRow, maxRowDeg, maxColDeg;
				string codeName;
				OneCodeInfo(int nR_, int nC_, int mRD, int mCD, string cN):
				nCol(nC_), nRow(nR_), maxRowDeg(mRD), maxColDeg(mCD), codeName(cN)
				{}
				OneCodeInfo():nCol(0), nRow(0), maxRowDeg(0),maxColDeg(0) {}
		};

		//int nCol, maxRowDeg, maxColDeg,nCodes;
		string ldpcDir;
		int nCodes;
		vector<OneCodeInfo> allCodes;

		vector<int> usedCodeIDs;
		//vector<int> nRowVec;
		//vector<string> codeName;
		//ALDPCInfo(vector<int>& nRow_, int nCol_, int maxRD_, int maxCD_, vector<string>& codeName_, int nCodes_,string& dir):
		//nRowVec(nRow_), nCol(nCol_), maxRowDeg(maxRD_), maxColDeg(maxCD_), codeName(codeName_),nCodes(nCodes_),ldpcDir(dir)
		//{}
		ALDPCInfo()
		{}
	};

	/** The parameters for the simulation. All members are public. Most members are the parameters. */
	class MSALPara: public SimuPara{
	public:
		
		int nbSize; ///< alphabet size.
		int nBits;///< number of bits for each symbol, log_2(nbSize)
		int nbCodeLen;///< code length.
		//int codePer;///< code period.
		//vector<int> usedCodesRaw; ///< codeID of the code that will be used, no need to sorted and unique. 

		vector<int> txRate;// tx rate for each source per block

		/** pdf parameter. #QARY 0.8 means a qary symmetric correlation with agree probability 0.8. Or it could be a filename that does 
		*	not start with #.
		*/
		string pdfPara;
		string decPdfPara;
		int seed; ///< seed for RNG.
		//string ldpcFN; ///< .lad file for LDPC codes.
		//string saveFN; ///< The file  name of the result file. 

		int maxRowDeg, maxColDeg;
		int maxErrFr; ///< maximum error frame (termination condition).
		int dispPerFr;///< The frequency of result display on screen.
		int saveTempPerFr;///< The frequency of writting the temp file.

		int rngTableOpt;///< 0, regularly compute the table (but not save), 1, load the table, 2, compute and save the table

		ALDPCInfo ldpcpara;

		MSALPara(): rngTableOpt(0)
		{}

		/** Load the parameters from a mapping, which is obtained from some form of user interface. The user interface type is irrelevent of the parameter
		*	object. 
		*	@param para a <string,string> (STL) mapping that contains parameters. 
		*	@param ui display user interface to display error message. 
		*/
		void loadPara(const paraMap& para, DispUI& ui);

		shared_ptr<SimuPara> clone() const
		{ return shared_ptr<MSALPara>(new MSALPara(*this)); }
	};

	/** The computed results from input parameters. This does not include the MC simulation results. */
	class MSALCompRes:public CompRes{
	public:
		vector<double> condEntro; ///< H(b0|X), H(b1|b0,X), ... H(b_i-1|b_0, ..., b_i-2, X), X is side info
	
		vector<double> codeRates; ///< k/n

		vector<int> nTxBits; ///< number of txed bits for each sub (binary) source.

		double totTxRate; ///< in bits/symbol.

		string startTime;///< simulation start time.

		/** Transfer the data (results) stored in the object to a mapping between the variable name and variable value (both string).
		*/
		void setMapping(compResMap& crMap) const;

		double HX/** H(X) */,HY/** H(Y) */, HYX/** H(Y|X)*/; 

		string codeIDstr;///<  codeIDs.

		MSALCompRes():totTxRate(0), HX(0), HY(0), HYX(0)
		{}
	};

	class LDPCGroup{
	
		vector<shared_ptr<LDPCDec> > vpLDPCDec;

	public:
		LDPCGroup(ALDPCInfo& ldpcinfo);
		LDPCGroup(){}
		LDPCDec* operator[](int i)
		{ return vpLDPCDec[i].get();}


	};
	/** Work space for each thread. */
	class ThreadWorkSpace{
	public:
		ThreadWorkSpace(int nbCodeLen, int nBits,const  vector<int>& nTxBits, ALDPCInfo& ldpcpara);

		LDPCGroup ldpcGroup;

	//	DLDPCEnDec ldpcCodec; ///< LDPC codec.
		vector<double> initLLR;////< Initial LLR.
		vector<vector<int> > binSourSeqs;///< binary source sequences.
		vector<vector<int> > decBinSeqs;///< decoded sequences
		vector<vector<int> > encSyns;///< encoded accumulated syndromes
		vector<int> isBSucc/** whether the binary source sequence decoded successful, 1 successful 0 otherwise*/, 
			nBErr/** number of erronous bits in each binary source */; // not a good idea to use vector<bool>
		bool isFSucc;///< Whether the frame is decoded correctly
	};

	/** The simulation results. */
	class MSALSimuRes{
	public:
		int totFr/** # of total simulated frames.*/, errFr/** # of error frames. */;
		double FER;
		vector<double> binSourBER; ///< BER for each binary source
		vector<double> binSourFER; ///< FER for each binary source
		vector<int> binErrFrame; ///< # of error frames for each binary source
		vector<int> binErrBit;///< # of total erroneous bits in each binary source
	};



	MSALSim(MSALPara& para_, TextSaveUI& resFile_, TextSaveUI& tempFile_, DispUI& dispUI_);

	void prepare(const MCThreadInfo&);//To do: May consider change the arugment to nThread

	void oneThread(shared_ptr<RNObj>, const MCThreadInfo&);//To do: May consider change the arguement to threadID

	

	shared_ptr<RNObj> genRNObj();//generate randome seq for EACH THREAD

	

	void dumpOneResult(const MCThreadInfo&);

	/** Single thread simulation. */
	//void singleThreadSim();//put into the base class

	/** Simulate one given frame, for testing or careful examination. */
	void oneKnownFrame(); //Testing only

	/** Calls compRes.setMapping(compResM_);.*/
	void setCompResMap(compResMap& compResM_) const;

	void initRNG();

	void postProcess();

	virtual bool checkTerm();

	
private:

	void oneStep(const vector<int>& xSeq, const vector<int>& ySeq, ThreadWorkSpace&);
	//void procThreadRes();

	// ****** Code Parameters ******
	int nbSize;
	int nbCodeLen;
	//int codePer; // period
	//vector<int> usedCodesRaw; // In terms of codeID, could have repetitive numbers, need not be ordered

	// ****** simulation parameters ******
	int maxErrFr; 

	// for implementation 
	pdfXY nbJoPdf; 

	/** length nbSize vector, each one corresponds to 
	* a joint pdf p(X|Y=y) for a fixed y. Each nbToBinPDFConv object convert pmf on large alphabet to 
	* binary sources. */
	vector<nbToBinPDFConv> nbBinConv;

	
	// Computed parameters

	int nBits;

	MSALPara para;
	MSALCompRes compRes;
	MSALSimuRes simuRes;
	TextSaveUI& resFile;
	TextSaveUI& tempFile;
	DispUI& dispUI;


	//LDPC codes
	//DLDPCCodes ldpcCodes;

	//workspaces for each thread, use heap storage (no need to use pointer actually).
	vector<shared_ptr<ThreadWorkSpace> > workspaceVec; // Warning! Copy not safe, the original and new shared pointers will refer to the same object
		
	// random number generator
	boost::mt19937 uniGen;
	boost::uniform_01<> uni_01;
	boost::variate_generator<boost::mt19937&, boost::uniform_01<> > uni01Gen;

	discrete_distribution xDistr; ///< The pmf of X.
	vector<discrete_distribution> yxCondDistr;///< The pmf of Y given X = x.

	
	LDPCGroup vpLDPCDec;

	//ALDPCInfo ldpcInfo;

};

















}
#endif
