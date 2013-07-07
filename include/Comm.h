#ifndef _COMM_H_INCLUDED
#define _COMM_H_INCLUDED
#include "StdHeaders.h"
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

namespace MTMCSim{
	namespace Comm{
		
		/** BPSK Baseband Modulator and Demodulator.
		*/
		class BaseBPSK
		{
		public:
			void modulate(const vector<int>& in, vector<double>& out);

			/** Soft output. relimat is the reliability matrix.*/
			void demGetReli(const vector<double>& in, vector<vector<double> >& reliMat, double sigma);
			BaseBPSK(int bitPerSym_);
			BaseBPSK()
			{}

			void setBitPerSym(int bitPerSym_);

		private:
			int bitPerSym;
			vector<vector<int> > symBitMap;

		};

		/** Baseband real AWGN channel.
		*/
		class BaseAWGNReal
		{
		public:
			/** normRV is a vector of standard normal RVs. */
			void pass(const vector<double>& in, vector<double>& out,const  vector<double>& normRV)
			{
				unsigned int len = in.size();
				for (unsigned int i = 0; i<len; i++)
					out[i] = in[i] + sigma * normRV[i];
			}

			/** Assume the input power is 1, SNR is in dB.*/
			BaseAWGNReal(double SNR_, double codeRate_):SNR(SNR_), codeRate(codeRate_)
			{
				sigma = sqrt(1.0 / pow(10,SNR/10) / codeRate / 2);
			}

			double getSigma()
			{
				return sigma;
			}

		private:
			double SNR;
			double sigma, codeRate;
		

			

		};




	}


}



#endif
