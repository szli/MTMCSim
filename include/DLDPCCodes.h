#ifndef _INCLUDED_DLDPCCODES_H
#define _INCLUDED_DLDPCCODES_H
#include "StdHeaders.h"

namespace MTMCSim{
using std::map;

/** Rate adaptive LDPC codes (for distributed source coding) management class.
*	Reads and process the LDPC codes parameter. The code is designed by David Varodayan, see more details:
*	http://ivms.stanford.edu/~varodayan/ldpca.html . And his paper. 
*
*	copysafe
*/
class DLDPCCodes{ 
friend class DLDPCEnDec;
public:
	DLDPCCodes();
	
	/** Constructor. 
	*	@param ladfn The file name of the .lad file that contains a series of rate-adatpive LDPC codes.
	*	@param [beg,end) The iterators in this range specifies the code IDs to be used. Usually it will be simple the returns
	*	from begin() and end() of a vector<int>
	*/
	template<typename InputIterator>
		DLDPCCodes(const string& ladfn, InputIterator beg, InputIterator end); // [beg,end) is the codes that will be used

	template<typename InputIterator>
		void setCodes(const string& ladfn, InputIterator beg, InputIterator end);
private:
	int codeLen, nzmax, period;
	vector<vector<int> > txSeqs, ir;
	vector<int> jc; ///< cumulative sum of column sum of H  (fixed across codes)
	int nCodesUsed; ///< number of codes used in current program (current object)
	vector<int> codeUsed;///< codeIDs, the codeID is the number of transmitted bits of a code. Sorted and no repeat.
	map<int,int> codeIDMapInd;///< codeID ==> Internal index (0, 1,..., nCodeUsed-1)
	vector<int> irWhole; ///< ir for the last code (rate 0), for encoding.
	int nCodes;
	
};

template<typename InputIterator>
	void DLDPCCodes::setCodes(const string& ladfn, InputIterator beg, InputIterator end)
	{
		ifstream ladF(ladfn.c_str());

		if (!ladF)
		throw runtime_error("Read File Error!" + ladfn);

		ladF >> nCodes >> codeLen >> nzmax >> period;

		if (!ladF)
			throw runtime_error("Error while reading LDPC file!");

		jc.resize(codeLen + 1);

		for (int i = 0; i< codeLen+1; i++)
			if (!(ladF >> jc[i]))
				throw runtime_error("Error while reading LDPC file!");

		codeUsed.resize(end-beg);

		std::copy(beg,end,codeUsed.begin()); // raw input

		sort(codeUsed.begin(), codeUsed.end());

		vector<int>::iterator it = std::unique(codeUsed.begin(), codeUsed.end());

		codeUsed.erase(it, codeUsed.end());

		// Now codeUsed is clean: sorted and contains unique IDs

		if ((codeUsed.front() < 2)||(codeUsed.back() > period))
			throw runtime_error("Wrong code ID!");

		nCodesUsed = static_cast<int>(codeUsed.size());

		for (int i = 0; i<nCodesUsed; i++)
			codeIDMapInd.insert(std::make_pair(codeUsed[i],i));

		txSeqs.resize(nCodesUsed);

		for (int i = 0; i<nCodesUsed; i++)
			txSeqs[i].resize(codeUsed[i]);

		ir.resize(nCodesUsed, vector<int>(nzmax));

		int reqCodeI = 0; 

		while (reqCodeI < nCodesUsed)
		{
			int codeID;
			if (!(ladF >> codeID))
				throw runtime_error("Error while reading LDPC file!");


			if (codeID != codeUsed[reqCodeI])
			{  // ignore this code (two rows)
				string str;
				ladF.ignore(numeric_limits<streamsize>::max(), '\n');
				getline(ladF,str);
				getline(ladF,str);
				if (!ladF)
					throw runtime_error("Error while reading LDPC file!");
				continue;
			}

			for (int j = 0; j<codeID; j++)
			ladF >> txSeqs[reqCodeI][j];

			for (int j = 0; j<nzmax; j++)
				ladF >> ir[reqCodeI][j];

			reqCodeI ++ ; // # codes that have been read
		}//end of while

		// let irWhole be the ir for the last code, rate 0 code
		if (codeUsed.back() == nCodes+1) // if the last code is used
			irWhole = ir.back();// irWhole is the ir for the last code
		else
		{
			// read the last code
			while(1)
			{
				int codeID;
				if (!(ladF >> codeID))
					throw runtime_error("Error while reading LDPC file!");

				string str;
				if (codeID != nCodes+1) // The last code
				{  // ignore this code (two rows)
					
					ladF.ignore(numeric_limits<streamsize>::max(), '\n');
					getline(ladF,str);
					getline(ladF,str);
					if (!ladF)
						throw runtime_error("Error while reading LDPC file!");
					continue;
				}
				
				ladF.ignore(numeric_limits<streamsize>::max(), '\n');
				getline(ladF,str); // ignore the txSeqs;
				irWhole.resize(nzmax);
				for (int j = 0; j<nzmax; j++)
					ladF >> irWhole[j];
				break;
			}
		}

		ladF.close();
	}




template<typename InputIterator>
	DLDPCCodes::DLDPCCodes(const string& ladfn, InputIterator beg, InputIterator end)
{
	setCodes(ladfn, beg, end);
}











	





		

		














}

#endif 