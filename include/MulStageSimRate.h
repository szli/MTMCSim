#ifndef _MULSTAGESIMRate_H_INCLUDED
#define _MULSTAGESIMRate_H_INCLUDED
#include "StdHeaders.h"
#include "RNObj.h"
#include <boost/shared_ptr.hpp>
#include "pdfXY.h"
#include "DLDPCCodes.h"
#include "discrete_distribution.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>
#include "nbToBinPDFConv.h"
#include <boost/pointer_cast.hpp>
#include "DLDPCEnDec.h"
#include "Util.h"
#include "TextSaveUI.h"
#include "DOSDispUI.h"
#include "MCThreadInfo.h"
#include "CompRes.h"
#include "SimuPara.h"
#include "MTMCSimBase.h"

using boost::shared_ptr;
using boost::static_pointer_cast;

namespace MTMCSim{

/** Simulator for multistage LDPC for two source distributed source coding. */
class MulStageSimRate: public MTMCSimBase{
public:
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

	/** The parameters for the simulation. All members are public. Most members are the parameters. */
	class MulStagePara: public SimuPara{
	public:
		int nbSize; ///< alphabet size.
		int nBits;///< number of bits for each symbol, log_2(nbSize)
		int nbCodeLen;///< code length.
		int codePer;///< code period.
		vector<int> usedCodesRaw; ///< codeID of the code that will be used, no need to sorted and unique. 
		vector<int> allCodes;
		/** pdf parameter. #QARY 0.8 means a qary symmetric correlation with agree probability 0.8. Or it could be a filename that does 
		*	not start with #.
		*/
		string pdfPara;
		string decPdfPara;
		int seed; ///< seed for RNG.
		string ldpcFN; ///< .lad file for LDPC codes.
		//string saveFN; ///< The file  name of the result file. 


		int maxErrFr; ///< maximum error frame (termination condition).
		int dispPerFr;///< The frequency of result display on screen.
		int saveTempPerFr;///< The frequency of writting the temp file.

		int rngTableOpt;///< 0, regularly compute the table (but not save), 1, load the table, 2, compute and save the table

		MulStagePara(): rngTableOpt(0)
		{}

		/** Load the parameters from a mapping, which is obtained from some form of user interface. The user interface type is irrelevent of the parameter
		*	object. 
		*	@param para a <string,string> (STL) mapping that contains parameters. 
		*	@param ui display user interface to display error message. 
		*/
		void loadPara(const paraMap& para, DispUI& ui);

		shared_ptr<SimuPara> clone() const
		{ return shared_ptr<MulStagePara>(new MulStagePara(*this)); }
	};

	/** The computed results from input parameters. This does not include the MC simulation results. */
	class MulStageCompRes:public CompRes{
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

		MulStageCompRes():totTxRate(0), HX(0), HY(0), HYX(0)
		{}
	};

	/** Work space for each thread. */
	class ThreadWorkSpace{
	public:
		ThreadWorkSpace(int nbCodeLen, int nBits,const  vector<int>& nTxBits, const DLDPCCodes& ldpcCodes);
		DLDPCEnDec ldpcCodec; ///< LDPC codec.
		vector<double> initLLR;////< Initial LLR.
		vector<vector<int> > binSourSeqs;///< binary source sequences.
		vector<vector<int> > decBinSeqs;///< decoded sequences
		vector<vector<int> > encSyns;///< encoded accumulated syndromes
		vector<int> isBSucc/** whether the binary source sequence decoded successful, 1 successful 0 otherwise*/, 
			nBErr/** number of erronous bits in each binary source */; // not a good idea to use vector<bool>
		bool isFSucc;///< Whether the frame is decoded correctly
		vector<int> actCodes;
	};

	/** The simulation results. */
	class MulStageSimuRes{
	public:
		int totFr;
		double totTxRate;
		double aveTxRate;
		vector<double> txRateVec;
		double stdev;

		MulStageSimuRes():totFr(0), totTxRate(0), aveTxRate(0){}
	
	};


	MulStageSimRate(MulStagePara& para_, TextSaveUI& resFile_, TextSaveUI& tempFile_, DispUI& dispUI_);

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

	void oneStep(const vector<int>& xSeq, const vector<int>& ySeq, ThreadWorkSpace&, DLDPCEnDec& );
	//void procThreadRes();

	// ****** Code Parameters ******
	int nbSize;
	int nbCodeLen;
	int codePer; // period
	vector<int> usedCodesRaw; // In terms of codeID, could have repetitive numbers, need not be ordered

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

	MulStagePara para;
	MulStageCompRes compRes;
	MulStageSimuRes simuRes;
	TextSaveUI& resFile;
	TextSaveUI& tempFile;
	DispUI& dispUI;


	//LDPC codes
	DLDPCCodes ldpcCodes;

	//workspaces for each thread, use heap storage (no need to use pointer actually).
	vector<shared_ptr<ThreadWorkSpace> > workspaceVec; // Warning! Copy not safe, the original and new shared pointers will refer to the same object
		
	// random number generator
	boost::mt19937 uniGen;
	boost::uniform_01<> uni_01;
	boost::variate_generator<boost::mt19937&, boost::uniform_01<> > uni01Gen;

	discrete_distribution xDistr; ///< The pmf of X.
	vector<discrete_distribution> yxCondDistr;///< The pmf of Y given X = x.


};

















}
#endif
