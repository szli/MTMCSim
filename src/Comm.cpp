#include "Comm.h"

namespace MTMCSim{
	namespace Comm{
	


		void BaseBPSK::modulate(const vector<int>& in, vector<double>& out)
		{
			unsigned int len = in.size();
			if (out.size() != len * bitPerSym)
				out.resize(len*bitPerSym);


			for (unsigned int i = 0; i<len; i++)
			{
				for (int j = 0; j<bitPerSym; j++)
					out[i*bitPerSym+j] = symBitMap[in[i]][j];
			}
		}

		void BaseBPSK::demGetReli(const vector<double>& in, vector<vector<double> >& reliMat, double sigma)
		{
			unsigned int nCol = reliMat[0].size();
			unsigned int symSize = 1<<bitPerSym;


			for (unsigned int col = 0; col<nCol; col++) // jth Column represents jth location
			{
				double colSum = 0;
				for (int i = 0; i<symSize; i++)
				{
					vector<int>& bits = symBitMap[i];
					double sum = 0;
					for (int j = 0; j<bitPerSym; j++)
					{
						sum += bits[j]==1? -(in[col*bitPerSym+j] - 1)*(in[col*bitPerSym+j] - 1) : 
										 -(in[col*bitPerSym+j] + 1)*(in[col*bitPerSym+j] + 1);
					}
					sum = sum / (2*sigma*sigma);

					reliMat[i][col] = exp(sum);
					colSum += reliMat[i][col];
				}

				for (int i = 0; i<symSize; i++)
					reliMat[i][col] /= colSum; // normalization
			}

		}

		BaseBPSK::BaseBPSK(int bitPerSym_):bitPerSym(bitPerSym_)
		{
				setBitPerSym(bitPerSym);
					
		}

		void BaseBPSK::setBitPerSym(int bitPerSym_)
		{
			bitPerSym = bitPerSym_;
			unsigned int symSize = 1<<bitPerSym;
				symBitMap.resize(symSize, vector<int>(bitPerSym,0));
				int temp = 0;
				for (unsigned int i = 0; i<symSize; i++)
				{
					temp = i;
					for (int j = 0; j<bitPerSym; j++)
					{
						symBitMap[i][j] = (temp % 2 == 1) ? 1 : -1;
						temp = temp >> 1;
					}
				}
		}










	
	
	}







}
