#ifndef _INCLUDED_DLDPCENDEC_H
#define _INCLUDED_DLDPCENDEC_H
#include "StdHeaders.h"
#include "DLDPCCodes.h"

namespace MTMCSim{

const int MAX_BP_ITER = 100;
/** LDPC encoder and decoder. Should be defined for each working thread. Copy safe.
*/
class DLDPCEnDec{
public:
		DLDPCEnDec()
		{}

		DLDPCEnDec(const DLDPCCodes& ldpcCodes_);

		void setLDPCCodes(const DLDPCCodes& ldpcCodes_);

		/** Encode a binary source sequence to accumulated syndrome using one code specified by codeID.*/
		void encodeOne(const vector<int>& binSourSeq,int codeID, vector<int>& accSyn);

		/** Decode an accumulated syndrome using Belief Propogation, hard output.*/
		bool decodeOne(const vector<int>& accSyn, int codeID, const vector<double>& initLLR, vector<int>& decSeq);

		/** Soft output of final LLR. The life time of the returned object is the same as the DLDPCEnDec object.*/
		const vector<double>& getLLR_overall() const
		{ return LLR_overall;}

private:
	vector<double> LLR_extrinsic, check_LLR, check_LLR_mag, rowTotal, LLR_overall;

	vector<int> synWhole;

	DLDPCCodes ldpcCodes;//To do: Consider to use reference (notice the setLDPCCodes function... if it is a reference, it has to be initialized in the constructor

	bool beliefPropagation(int *ir, int *jc, int m, int n, int nzmax, 
                      const double *LLR_intrinsic, int *syndrome,
                       int *decoded);

};







}
#endif 