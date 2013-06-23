#ifndef _KVCHANSIM_H_INCLUDED
#define _KVCHANSIM_H_INCLUDED
#include "StdHeaders.h"
#include "RNObj.h"
#include <boost/shared_ptr.hpp>
#include "pdfXY.h"
#include "discrete_distribution.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/pointer_cast.hpp>
#include "Util.h"
#include "TextSaveUI.h"
#include "DOSDispUI.h"
#include "MCThreadInfo.h"
#include "CompRes.h"
#include "SimuPara.h"
#include "MTMCSimBase.h"
#include "MultiRunnable.h"
#include "MultiRunSimuBase.h"
#include "Comm.h"

using boost::shared_ptr;
using boost::static_pointer_cast;

namespace MTMCSim
{
/** Koetter-Vardy decoding simulation for channel coding. 
*	The simulation is for channel coding scenario.
*/
class KVChanSim: public MultiRunSimuBase
{
public:
	/** A random number object for this problem is one source sequence and one channel noise sequence (standard normal).
	*/
	class XNSeq: public RNObj
	{
		public:
		vector<int> xSeq;// X is side info
		vector<double> normRV;
		XNSeq(int xlen, int normLen):xSeq(xlen), normRV(normLen)
		{ }
	};

	/** The parameters for the simulation. All members are public. Most members are the parameters. */
	class KVChanPara: public SimuPara{
	public:
		int nbSize; ///< alphabet size.
		int codeLen;///< code length.
		int infoLen;///< information length, k.
		
		double lambda; ///< proportional multiplicity assignment parameter
		
		double SNR;

		int seed; ///< seed.

		int chanSeed;

		int maxErrFr; ///< maximum error frame (termination condition).
		int dispPerFr;///< The frequency of result display on screen.
		int saveTempPerFr;///< The frequency of writting the temp file.
		int mulMode;
		
		/** Load the parameters from a mapping, which is obtained from some form of user interface. The user interface type is irrelevent of the parameter
		*	object. 
		*	@param para a <string,string> (STL) mapping that contains parameters. 
		*	@param ui display user interface to display error message. 
		*/
		void loadPara(const paraMap& para, DispUI& ui);

		shared_ptr<SimuPara> clone() const
		{ return shared_ptr<KVChanPara>(new KVChanPara(*this)); }

	};

	/** The computed results from input parameters. This does not include the MC simulation results. */
	class KVChanCompRes:public CompRes{
	public:
		
		string startTime;///< simulation start time.

		/** Transfer the data (results) stored in the object to a mapping between the variable name and variable value (both string).
		*/
		void setMapping(compResMap& crMap) const;

		KVChanCompRes()
		{}
	};

	/** Work space for each thread. */
	class ThreadWorkSpace{
	public:
		vector<double> channelIn; ///< Modulated signal.
		vector<double> channelOut;///<channel output.
		vector<vector<double> > reliMat;///< reliablity matrix.
		vector<vector<int> > MulMat;///< Multiplicity matrix.
		bool isFSucc;///< Whether the frame is decoded correctly
		ThreadWorkSpace(int codeLen, int bitLen, int nbSize):channelIn(bitLen), channelOut(bitLen), 
			reliMat(nbSize, vector<double>(codeLen)), MulMat(nbSize, vector<int>(codeLen))
		{}
			
	};

	/** The simulation results. */
	class KVChanSimuRes{
	public:
		int totFr/** # of total simulated frames.*/, errFr/** # of error frames. */;
		double FER;
		KVChanSimuRes():totFr(0), errFr(0), FER(0)
		{}
	};

	KVChanSim(KVChanPara& para_, TextSaveUI& resFile_, TextSaveUI& tempFile_, DispUI& dispUI_);

	void prepare(const MCThreadInfo&);

	void oneThread(shared_ptr<RNObj>, const MCThreadInfo&);

	shared_ptr<RNObj> genRNObj();//generate randome seq for EACH THREAD

	void dumpOneResult(const MCThreadInfo&);


	/** Calls compRes.setMapping(compResM_);.*/
	void setCompResMap(compResMap& compResM_) const;



	void postProcess();

	virtual bool checkTerm();


private:

	void oneStep(const vector<int>& xSeq, const vector<double>& normRV, ThreadWorkSpace&);

	int KVLooseBound(const int &cost);

	int KVTightBound(const int &cost);

	void getCompRes();

	void initRNG();

	/** Multiplicity assignment scheme in KV's paper.
	*	If mode == 0, limit sum mode, para is the max sum multiplicity
	*	If mode == 1, limit y-degree, pare is the max y-degree
	*	yweight is the infomation length for RS,  infoLen-1 + g  for Hermitian. 
	*/
	int KVMul(vector<vector<double> >& reliMat, vector<vector<int> >& MulMat, int para, int mode, int yweight); 

	int HemitianLB(int cost);

	int nbSize, nBits;
	int codeLen, infoLen;

	double lambda;

	// ****** simulation parameters ******
	int maxErrFr; 

	// for implementation 
	pdfXY nbJoPdf; 
	
	KVChanPara para;
	KVChanCompRes compRes;
	KVChanSimuRes simuRes;
	TextSaveUI& resFile;
	TextSaveUI& tempFile;
	DispUI& dispUI;

	//workspaces for each thread, use heap storage (no need to use pointer actually).
	vector<shared_ptr<ThreadWorkSpace> > workspaceVec; // Warning! Copy not safe, the original and new shared pointers will refer to the same object
	

	// random number generator
	boost::mt19937 uniGen;
	boost::uniform_smallint<> uniIntGen;
	boost::normal_distribution<> normGen;

	Comm::BaseBPSK BPSKmodem;
	Comm::BaseAWGNReal channel;

};



}

#endif