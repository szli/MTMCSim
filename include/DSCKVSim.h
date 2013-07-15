#ifndef _DSCKV_H_INCLUDED
#define _DSCKV_H_INCLUDED
#include "StdHeaders.h"
#include "RNObj.h"
#include <boost/shared_ptr.hpp>
#include "pdfXY.h"
#include "discrete_distribution.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>
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

using boost::shared_ptr;
using boost::static_pointer_cast;

namespace MTMCSim
{
class DSCKVSim: public MultiRunSimuBase
{
public:
    /** A random number object for this problem is two correlated sequence X and Y. Note X is the side 
    *   information and Y is the source to be transmitted.
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
    class DSCKVPara: public SimuPara{
    public:
        int nbSize; ///< alphabet size.
        int codeLen;///< code length.
        int infoLen;///< information length, k.
        
        double lambda; ///< proportional multiplicity assignment parameter
        /** pdf parameter. #QARY 0.8 means a qary symmetric correlation with agree probability 0.8. Or it could be a filename that does 
        *   not start with #.
        */
        string pdfPara;
        string decPdfPara;
        int seed; ///< seed for RNG.
    //  string saveFN; ///< The file  name of the result file. 

        int maxErrFr; ///< maximum error frame (termination condition).
        int dispPerFr;///< The frequency of result display on screen.
        int saveTempPerFr;///< The frequency of writting the temp file.
        int mulMode;
        string algoMode;
        int rngTableOpt; ///< 0, regularly compute the table (but not save), 1, load the table, 2, compute and save the table

        /** Load the parameters from a mapping, which is obtained from some form of user interface. The user interface type is irrelevent of the parameter
        *   object. 
        *   @param para a <string,string> (STL) mapping that contains parameters. 
        *   @param ui display user interface to display error message. 
        */
        void loadPara(const paraMap& para, DispUI& ui);

        shared_ptr<SimuPara> clone() const
        { return shared_ptr<DSCKVPara>(new DSCKVPara(*this)); }

        DSCKVPara():rngTableOpt(0), algoMode("KVDSC") ///< Default constructor is used to set default values
        {}
    };

    /** The computed results from input parameters. This does not include the MC simulation results. */
    class DSCKVCompRes:public CompRes{
    public:
        
        double txRate; ///< in bits/symbol.
        string startTime;///< simulation start time.
        int infoLen;///< Make it a computed result because it may change during multi-run mode

        /** Transfer the data (results) stored in the object to a mapping between the variable name and variable value (both string).
        */
        void setMapping(compResMap& crMap) const;

        double HX/** H(X) */,HY/** H(Y) */, HYX/** H(Y|X)*/; 

        double gap; 

        DSCKVCompRes():txRate(0), HX(0), HY(0), HYX(0)
        {}
    };

    /** Work space for each thread. */
    class ThreadWorkSpace{
    public:
        bool isFSucc;///< Whether the frame is decoded correctly
    };

    /** The simulation results. */
    class DSCKVSimuRes{
    public:
        int totFr/** # of total simulated frames.*/, errFr/** # of error frames. */;
        double FER;
        DSCKVSimuRes():totFr(0), errFr(0), FER(0)
        {}
    };

    DSCKVSim(DSCKVPara& para_, TextSaveUI& resFile_, TextSaveUI& tempFile_, DispUI& dispUI_);

    void prepare(const MCThreadInfo&);

    void oneThread(shared_ptr<RNObj>, const MCThreadInfo&);

    shared_ptr<RNObj> genRNObj();//generate randome seq for EACH THREAD

    void dumpOneResult(const MCThreadInfo&);


    /** Calls compRes.setMapping(compResM_);.*/
    void setCompResMap(compResMap& compResM_) const;



    void postProcess();

    virtual bool checkTerm();


private:

    void oneStep(const vector<int>& xSeq, const vector<int>& ySeq, ThreadWorkSpace&);

    double KVLooseBound(const double &cost);

    double KVTightBound(const double &cost);

    double HemitianLB(double cost);

    double HermitianLB2(double cost);

    void getCompRes();

    void initRNG();

    int nbSize;
    int codeLen, infoLen;

    double lambda;

    // ****** simulation parameters ******
    int maxErrFr; 

    // for implementation 
    pdfXY nbJoPdf; 
    
    DSCKVPara para;
    DSCKVCompRes compRes;
    DSCKVSimuRes simuRes;
    TextSaveUI& resFile;
    TextSaveUI& tempFile;
    DispUI& dispUI;

    //workspaces for each thread, use heap storage (no need to use pointer actually).
    vector<shared_ptr<ThreadWorkSpace> > workspaceVec; // Warning! Copy not safe, the original and new shared pointers will refer to the same object
    

    // random number generator
    boost::mt19937 uniGen;
    boost::uniform_01<> uni_01;
    boost::variate_generator<boost::mt19937&, boost::uniform_01<> > uni01Gen;

    discrete_distribution xDistr; ///< The pmf of X.
    vector<discrete_distribution> yxCondDistr;///< The pmf of Y given X = x.

    // Pre-determined multiplicity matrix
    vector<vector<int> > PreMul;
};



}

#endif