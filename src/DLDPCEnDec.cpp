#include "DLDPCEnDec.h"

namespace MTMCSim{


DLDPCEnDec::DLDPCEnDec(const DLDPCCodes &ldpcCodes_)
	:ldpcCodes(ldpcCodes_)
{
	LLR_extrinsic.resize(ldpcCodes.nzmax);
	check_LLR.resize(ldpcCodes.nzmax);
	check_LLR_mag.resize(ldpcCodes.nzmax);
	rowTotal.resize(ldpcCodes.codeLen);
	LLR_overall.resize(ldpcCodes.codeLen);
	synWhole.resize(ldpcCodes.codeLen);// syn len = n
	
}

void DLDPCEnDec::setLDPCCodes(const DLDPCCodes& ldpcCodes_)
{
	ldpcCodes = ldpcCodes_;
	LLR_extrinsic.resize(ldpcCodes.nzmax);
	check_LLR.resize(ldpcCodes.nzmax);
	check_LLR_mag.resize(ldpcCodes.nzmax);
	rowTotal.resize(ldpcCodes.codeLen);
	LLR_overall.resize(ldpcCodes.codeLen);
	synWhole.resize(ldpcCodes.codeLen);// syn len = n
	
}


void DLDPCEnDec::encodeOne(const vector<int>& binSourSeq, int codeID, vector<int>& accSyn)
{
	// first compute the syndome of the rate 0 code, the last code
	// This should be done for every binary source since the source sequence is different

	for (int i = 0; i<ldpcCodes.codeLen; i++)
		synWhole[i] = 0;

	for (int i = 0; i<ldpcCodes.codeLen; i++)
		for (int j = ldpcCodes.jc[i]; j<ldpcCodes.jc[i+1]; j++)
			synWhole[ldpcCodes.irWhole[j]] += binSourSeq[i];// Now it is the syndrome of a rate 0 code

	int synLen = (ldpcCodes.codeLen/ldpcCodes.period) * codeID; //codeID is the synLen per period

	int codeIndex = ldpcCodes.codeIDMapInd[codeID];

	// Accumulate

	int prevIndex, currIndex;

	prevIndex = 0;

	for (int i = 0; i<synLen; i++)
	{
		currIndex = i/codeID * ldpcCodes.period + ldpcCodes.txSeqs[codeIndex][i % codeID];
		
		accSyn[i] = i == 0? 0: accSyn[i-1];

		for (int j = prevIndex; j<= currIndex; j++)
			accSyn[i] += synWhole[j];


		accSyn[i] %= 2;

		prevIndex = currIndex + 1;
	}
}

bool DLDPCEnDec::decodeOne(const vector<int>& accSyn, int codeID, const vector<double>& initLLR, vector<int>& decSeq)
{
	if (codeID == ldpcCodes.nCodes + 1)
		return true;

	int codeIndex = ldpcCodes.codeIDMapInd[codeID];

	int synLen = (ldpcCodes.codeLen / ldpcCodes.period) * codeID;

	// get the syndrome from the accumulated syndrome

	vector<int> syn(synLen);

	int currIndex = 0, prevIndex = 0;

	syn[0] = accSyn[0];

	for (int i = 1; i<synLen; i++)
	{
		syn[i] =  ((accSyn[i] + accSyn[i-1]) % 2);
	}

	return beliefPropagation(&ldpcCodes.ir[codeIndex][0], &ldpcCodes.jc[0],
		synLen, ldpcCodes.codeLen, ldpcCodes.nzmax, &initLLR[0], &syn[0], &decSeq[0]);
}




//From David
//For implementation outline of beliefPropagation(), refer to 
//W. E. Ryan, "An Introduction to LDPC Codes," in CRC Handbook for Coding 
//and Signal Processing for Recording Systems (B. Vasic, ed.) CRC Press, 2004.
//available online (as of May 8, 2006) at 
//http://www.ece.arizona.edu/~ryan/New%20Folder/ryan-crc-ldpc-chap.pdf

//beliefPropagation() runs several iterations belief propagation until
//either the decoded bitstream agrees with the transmitted portion of 
//accumulated syndrome or convergence or the max number of iterations.
//Returns 1 if decoded bitstream agrees with 
//transmitted portion of accumulated syndrome.
bool DLDPCEnDec::beliefPropagation(int *ir, int *jc, int m, int n, int nzmax, 
                       const double *LLR_intrinsic, int *syndrome,
                       int *decoded)  // by szli, change syndrome decoded to int*
{
    int iteration, k, l, sameCount;
    
    sameCount = 0;
    for(k=0; k<n; k++)
        decoded[k] = 0;
    
    //initialize variable-to-check messages
    for(k=0; k<n; k++)
        for(l=jc[k]; l<jc[k+1]; l++)
            LLR_extrinsic[l] = LLR_intrinsic[k];
    
    for(iteration=0; iteration<MAX_BP_ITER; iteration++)
    {
        //Step 1: compute check-to-variable messages
        
        for(k=0; k<nzmax; k++)
        {
            check_LLR[k] = (double) ((LLR_extrinsic[k]<0) ? -1 : 1);
            check_LLR_mag[k] = ((LLR_extrinsic[k]<0) ? -LLR_extrinsic[k] : LLR_extrinsic[k]);
        }
        
        for(k=0; k<m; k++)
            rowTotal[k] = (double) ((syndrome[k]==1) ? -1 : 1);
        for(k=0; k<nzmax; k++)
            rowTotal[ir[k]] *= check_LLR[k];        
        for(k=0; k<nzmax; k++)
            check_LLR[k] = check_LLR[k] * rowTotal[ir[k]];
            //sign of check-to-variable messages
        
        for(k=0; k<nzmax; k++)
			check_LLR_mag[k] = -log( tanh( std::max(check_LLR_mag[k], 0.000000001)/2 ) );
        for(k=0; k<m; k++)
            rowTotal[k] = (double) 0;
        for(k=0; k<nzmax; k++)
            rowTotal[ir[k]] += check_LLR_mag[k];        
        for(k=0; k<nzmax; k++)
			check_LLR_mag[k] = -log( tanh( std::max(rowTotal[ir[k]] - check_LLR_mag[k], 0.000000001)/2 ) );
            //magnitude of check-to-variable messages
            
        for(k=0; k<nzmax; k++)
            check_LLR[k] = check_LLR[k] * check_LLR_mag[k];
            //check-to-variable messages
            
        //Step 2: compute variable-to-check messages
        
        for(k=0; k<n; k++)
        {
            LLR_overall[k] = LLR_intrinsic[k];
            for(l=jc[k]; l<jc[k+1]; l++)
                LLR_overall[k] += check_LLR[l];
        }
            
        for(k=0; k<n; k++)
            for(l=jc[k]; l<jc[k+1]; l++)
                LLR_extrinsic[l] = LLR_overall[k] - check_LLR[l];
                //variable-to-check messages
            
        //Step 3: test convergence and syndrome condition
        
        l = 0;
        for(k=0; k<n; k++)
            if(decoded[k] == ((LLR_overall[k]<0) ? 1 : 0))
                l++;
            else
                decoded[k] = ((LLR_overall[k]<0) ? 1 : 0);
        
        sameCount = ((l==n) ? sameCount+1 : 0); 
        
        if(sameCount==5)

		{
            return 0; //convergence (to wrong answer)
		}
        
        for(k=0; k<m; k++)
            rowTotal[k] = syndrome[k];
        for(k=0; k<n; k++)
            for(l=jc[k]; l<jc[k+1]; l++)
                rowTotal[ir[l]] += decoded[k];
                
        for(k=0; k<m; k++)
            if((static_cast<int>(rowTotal[k]+0.5) % 2) != 0) // was (int)(rowTotal[k]), modified Sep 15,2010
                break;
            else if(k==m-1)
			{
				
                return 1; //all syndrome checks satisfied
			}
           
    }
	
    return 0;
}












	









}