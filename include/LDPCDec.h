#ifndef LDPC_DEC_H
#define LDPC_DEC_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <stdexcept>
//#include "LDPCDriver.h"
using namespace std;


#define NEARCW_START_ITER 45 //Keep tracking till the end
#define SMALL_ERR_SIZE 12

//For the phi function
#define SMALL_LLR 1e-10
#define LARGE_LLR 23.719


class LDPCDec
{
public:	
	int nRows;
	int nCols;
	int maxRowDeg;
	int maxColDeg;
	int no_edges;

	int *HRowProfile;
	int *HColProfile;
	int **HRows;

	int *cum_check_edges;
	int *cum_var_edges;

	int *check_to_var_perm;
	int *var_to_check_perm;
	
	double *edge_llrs_rows;
	double *edge_llrs_cols;

	double *channel_llrs;


	
	~LDPCDec(void);
	LDPCDec(int _nRows, int _nCols, int _maxRowDeg, int _maxColDeg, char* HRowProfileFile, char* HRowFile, char* HColProfileFile);
	void var_node_proc(int curr_node_index);
	void check_node_proc(int curr_node_index);
	void check_node_proc(int curr_node_index, int coset_bit);

	void compute_perm();
	int syndrome_decode(int *tx_cw, double* _channel_llrs, int maxIter, int *decodedCW);
	int decode(double* _channel_llrs, int maxIter, int *decodedCW);
	int decode_nearCW(double* _channel_llrs, int maxIter, int *decodedCW, char* nearCW_file, char* small_err_file);
	int decode_nearCW(double* _channel_llrs, int maxIter, int *decodedCW, char* nearCW_file);
	int coset_decode_nearCW(int *tx_cw,double* _channel_llrs, int maxIter, int *decodedCW, char* nearCW_file, char* small_err_file);
	int coset_decode_nearCW(int *tx_cw,double* channel_llrs,int maxIter,int* decodedCW,char* nearCW_file);
	int compute_syndrome(int *hardDecision);
	int compute_syndrome(int *hardDecision, int *coset);
	int final_decision_proc(int curr_node_index);
	double sgn(double x);
private:
	double phi(double x);
	
	void reset_llrs();
};
#endif
