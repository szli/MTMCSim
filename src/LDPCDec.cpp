#include "LDPCDec.h"
#include <cstdlib>

LDPCDec::LDPCDec(int _nRows, int _nCols, int _maxRowDeg, int _maxColDeg, char* HRowProfileFile, char* HRowFile, char* HColProfileFile)
{
	char *ch_num;
	int i;
	int k;
	int ct;
	int nRowsPlus1;
	int nColsPlus1;
	int curr_row;
	ifstream fp_row_profile;
	ifstream fp_col_profile;
	ifstream fp_rows;
	//ifstream fp;

	nRows = _nRows;
	nCols = _nCols;
	maxRowDeg = _maxRowDeg;
	maxColDeg = _maxColDeg;

	ch_num = new char[100];
	//Input the row profile file into the memory
	HRowProfile = new int[nRows];

	fp_row_profile.open(HRowProfileFile);
	

    if (!fp_row_profile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }

#ifdef DEBUG
	cout << "thru row profile file open" <<endl;
#endif
	ct = 0;
	while (fp_row_profile >> ch_num) {
		HRowProfile[ct] = atoi(ch_num);

		if (HRowProfile[ct] > maxRowDeg)
			throw runtime_error("Actual row degree greater than specified maxRowDeg!");

		ct++;
	}
	#ifdef DEBUG
	cout << "Number of row entries:" << ct;
#endif
	fp_row_profile.close();


	//Compute the cumulative sum of the number of rows
	nRowsPlus1 = nRows + 1;
	cum_check_edges = new int[nRowsPlus1];
	
	cum_check_edges[0] = 0;
	for (i=1;i<nRowsPlus1;i++) {
		cum_check_edges[i] = cum_check_edges[i-1] + HRowProfile[i-1];
	}


	//Input the column profile file into the memory
	HColProfile = new int[nCols];

	fp_col_profile.open(HColProfileFile);

    if (!fp_col_profile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
	#ifdef DEBUG
	cout << "thru col profile file open" <<endl;
#endif
	ct = 0;
	while (fp_col_profile >> ch_num) {
		HColProfile[ct] = atoi(ch_num);
		ct++;
	}
	#ifdef DEBUG
	cout << "Number of col entries:" << ct;
#endif
	fp_col_profile.close();

	//Compute the cumulative sum of the number of cols
	nColsPlus1 = nCols + 1;
	cum_var_edges = new int[nColsPlus1];
	
	cum_var_edges[0] = 0;
	for (i=1;i<nColsPlus1;i++) {
		cum_var_edges[i] = cum_var_edges[i-1] + HColProfile[i-1];
	}



	//Input the rows into the memory
	HRows = new int*[nRows];

	for (i=0;i<nRows;i++)
		HRows[i] = new int[maxRowDeg];

	fp_rows.open(HRowFile);

	if (!fp_rows) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
	#ifdef DEBUG
	cout << "thru row file open" <<endl;
#endif 
	ct = 0;
	curr_row = 0;
	k = 0;
	while (fp_rows >> ch_num) {
		if (ct == cum_check_edges[curr_row + 1]) {
			curr_row++;
			k = 0;
		}
		HRows[curr_row][k] = atoi(ch_num);
		k++;
		ct++;
	}

	fp_rows.close();

	//Instantiate edge permutation and de_permutation
	no_edges = cum_check_edges[nRowsPlus1 - 1];
	check_to_var_perm = new int[no_edges];
	var_to_check_perm = new int[no_edges];

	//Instantiate llr's
	edge_llrs_rows = new double[no_edges]; //LLR's from check to var nodes
	edge_llrs_cols = new double[no_edges]; //LLR's from var to check nodes


	//Compute the permutations
	compute_perm();

	//Instantiate the channel llrs
	channel_llrs = new double[nCols];


	//clean up

	delete ch_num;
}

void LDPCDec::compute_perm()
{
	int *count_cols;
	int ct;
	int i;
	int j;
	//int *sanity_check;

	count_cols = new int[nCols];
	//sanity_check = new int[no_edges];

	//for (i=0;i<no_edges;i++)
	//	sanity_check[i] = 0;

	for (i=0;i<nCols;i++)
		count_cols[i] = 0;
	//Computing the permutation from the check to the var nodes 
    ct = 0;
	for (i=0;i<nRows;i++) {
		//if (i==414)
		//	cout << endl;
		
		for (j=0;j<HRowProfile[i];j++) {
			check_to_var_perm[ct] = cum_var_edges[HRows[i][j]] + count_cols[HRows[i][j]];
			count_cols[HRows[i][j]]++;
			ct++;
		}
	}

	//Now compute the permutation from the var to the check nodes

	for (i=0;i<no_edges;i++) {
		for (j=0;j<no_edges;j++) {
			if (check_to_var_perm[j] == i) {
				break;
			}
		}
		var_to_check_perm[i] = j;
	}


	//for (i=0;i<no_edges;i++) {
	//	sanity_check[check_to_var_perm[i]]++;
	//}

	//for (i=0;i<no_edges;i++) {
	//	if (sanity_check[i] > 1) {
	//		cout << "Problem" << endl;
	//	}
	//}


	//clean up
	delete count_cols;
	//delete sanity_check;

}
void LDPCDec::var_node_proc(int curr_node_index) 
{
	int i;
	double s;

	s = 0.0;
	s += channel_llrs[curr_node_index];

	for (i=0;i<HColProfile[curr_node_index];i++) {
		s += edge_llrs_rows[var_to_check_perm[cum_var_edges[curr_node_index] + i]];
	}

	for (i=0;i<HColProfile[curr_node_index];i++) {
		edge_llrs_cols[cum_var_edges[curr_node_index] + i] = s - edge_llrs_rows[var_to_check_perm[cum_var_edges[curr_node_index] + i]];
	}
}

int LDPCDec::final_decision_proc(int curr_node_index)
{
	int i;
	double s;
	int val;

	s = channel_llrs[curr_node_index];
	for (i=0;i<HColProfile[curr_node_index];i++) {
		s += edge_llrs_rows[var_to_check_perm[cum_var_edges[curr_node_index] + i]];
	}
	//Important change 09/03/07. Using the mapping
	// c_i -> 1 - 2*c_i.
	val = (s < 0.0);
	//if (s < 0.0) {
	//	//cout << "(" << curr_node_index << "," << s << ")" << " ";
	//	cout << curr_node_index << " ";
	//}

	return val;
}

void LDPCDec::check_node_proc(int curr_node_index)
{
	//Implements the Sum-Product algorithm currently
	//need to put in more approximations at a later date
	int i;
	double s;
	double p;

	s = 0.0;
	p = 1.0;

	for (i=0;i<HRowProfile[curr_node_index];i++) {
		s += phi(fabs(edge_llrs_cols[check_to_var_perm[cum_check_edges[curr_node_index] + i]]));
		p *= sgn(edge_llrs_cols[check_to_var_perm[cum_check_edges[curr_node_index] + i]]);
	}

	for (i=0;i<HRowProfile[curr_node_index];i++) {
		edge_llrs_rows[cum_check_edges[curr_node_index] + i] = (p/sgn(edge_llrs_cols[check_to_var_perm[cum_check_edges[curr_node_index] + i]])) * phi(s - phi(fabs(edge_llrs_cols[check_to_var_perm[cum_check_edges[curr_node_index] + i]])));
	}
	
}


void LDPCDec::check_node_proc(int curr_node_index, int coset_bit)
{
	//Implements the Sum-Product algorithm currently
	//need to put in more approximations at a later date
	int i;
	double s;
	double p;

	s = 0.0;

	//if the coset_bit is a 1, then there is a extra sign flip
	if (coset_bit == 0) {
		p = 1.0;
	}
	else {
		p = -1.0;
	}

	for (i=0;i<HRowProfile[curr_node_index];i++) {
		s += phi(fabs(edge_llrs_cols[check_to_var_perm[cum_check_edges[curr_node_index] + i]]));
		p *= sgn(edge_llrs_cols[check_to_var_perm[cum_check_edges[curr_node_index] + i]]);
	}

	for (i=0;i<HRowProfile[curr_node_index];i++) {
		edge_llrs_rows[cum_check_edges[curr_node_index] + i] = (p/sgn(edge_llrs_cols[check_to_var_perm[cum_check_edges[curr_node_index] + i]])) * phi(s - phi(fabs(edge_llrs_cols[check_to_var_perm[cum_check_edges[curr_node_index] + i]])));
	}
	
}

double LDPCDec::phi(double x)
{
	double tmp;

	if (x < 0.0) {
		cout << "Error phi has a negative argument" << endl;
		exit(1);
	}
	if (x < SMALL_LLR) {
		return LARGE_LLR;
	}
	else if (x > LARGE_LLR) {
		return SMALL_LLR;
	}
	else {	
		tmp = (exp(x) + 1)/(exp(x) - 1);
		return log(tmp);
	}
}

double LDPCDec::sgn(double x)
{
	if (x >= 0.0) {
		return 1.0;
	}
	else {return -1.0;}
}

int LDPCDec::compute_syndrome(int *hardDecision)
{
	int i;
	int j;
	int tmp;
	int s;

	s = 0;
	tmp = 0;

	for (i=0;i<nRows;i++) {
		tmp = 0;
		for (j=0;j<HRowProfile[i];j++) {
			tmp += hardDecision[HRows[i][j]];
		}
		s += tmp%2;
	}

	return s;
}

int LDPCDec::compute_syndrome(int *hardDecision, int *coset)
{
	int i;
	int j;
	int tmp;
	int s;

	s = 0;
	tmp = 0;

	for (i=0;i<nRows;i++) {
		tmp = 0;
		for (j=0;j<HRowProfile[i];j++) {
			tmp += hardDecision[HRows[i][j]];
		}
		coset[i] = tmp%2;
		s += tmp%2;
	}

	return s;
}

void LDPCDec::reset_llrs()
{
	int i;

	for (i=0;i<no_edges;i++) {
		edge_llrs_rows[i] = 0.0;
		edge_llrs_cols[i] = 0.0;
	}
}
int LDPCDec::decode(double* _channel_llrs, int maxIter, int *decodedCW)
{
	int i;
	int j;
	int synd_count;
	int *hardDecision;
	
	hardDecision = new int[nCols];

	for (i=0;i<nCols;i++)
		channel_llrs[i] = _channel_llrs[i];

	//Compute the hard decision corresponding to the bits
	for (i=0;i<nCols;i++){
		//Important change 09/03/07. Using the mapping
		// c_i -> 1 - 2*c_i.
		hardDecision[i] = (channel_llrs[i] < 0.0);
		//hardDecision[i] = (channel_llrs[i] > 0.0);
	}

	//Compute the syndrome based on the hardDecision
	synd_count = compute_syndrome(hardDecision);

	if (synd_count == 0) {

		for (i=0;i<nCols;i++) {
			decodedCW[i] = hardDecision[i];
		}
		//cout << "Decoding success" << endl;
		return 1;
	}

	reset_llrs(); //may have to change if we implement turbo-eq


	for (i=0;i<maxIter;i++) {
		
 		for (j=0;j<nCols;j++) {
			//if (j==151) {
			//	cout << endl;
			//}
			var_node_proc(j); //processing the j-th var node
		}
		//At this point edge_llrs_cols should be updated

		for (j=0;j<nRows;j++) {
			//if ((j==26) || (j==29)) {
			//	cout << endl;
			//}
			check_node_proc(j);
		}
		//At this point edge_llrs_rows should be updated
		//cout << "Iteration number:" << i << endl;
		for (j=0;j<nCols;j++) {
			hardDecision[j] = final_decision_proc(j);
		}
		//cout << endl << endl << endl;
		synd_count = compute_syndrome(hardDecision);

		if (synd_count == 0) 
			break;
	}

	for (i=0;i<nCols;i++)
		decodedCW[i] = hardDecision[i];

	//clean up
	delete hardDecision;

	if (synd_count == 0) {
		//cout << "Decoding success" << endl;
		return 1;
	}
	else {
		//cout << "Decoding Failure" << endl;
		return -1;
	}

}

int LDPCDec::syndrome_decode(int *tx_cw, double* _channel_llrs, int maxIter, int *decodedCW)
{
	int i;
	int j;
	int synd_count;
	int *hardDecision;
	int tmp=0;
	int *tx_cw_coset;
	int *dec_coset;
	
	hardDecision = new int[nCols];
	tx_cw_coset = new int[nRows];
	dec_coset = new int[nRows];

	//First compute the coset of the transmitted codeword
	compute_syndrome(tx_cw,tx_cw_coset);

	for (i=0;i<nCols;i++)
		channel_llrs[i] = _channel_llrs[i];

	//Compute the hard decision corresponding to the bits
	for (i=0;i<nCols;i++){
		//Important change 09/03/07. Using the mapping
		// c_i -> 1 - 2*c_i.
		hardDecision[i] = (channel_llrs[i] < 0.0);
		//hardDecision[i] = (channel_llrs[i] > 0.0);
	}

	//Compute the syndrome based on the hardDecision
	synd_count = compute_syndrome(hardDecision,dec_coset);

	tmp = 0;
	for (i=0;i<nRows;i++) {
		if (dec_coset[i] != tx_cw_coset[i]) {
			tmp++;
		}
	}

	if (tmp == 0) {

		for (i=0;i<nCols;i++) {
			decodedCW[i] = hardDecision[i];
		}
		//cout << "Decoding success" << endl;
		return 1;
	}

	reset_llrs(); //may have to change if we implement turbo-eq


	for (i=0;i<maxIter;i++) {
		
 		for (j=0;j<nCols;j++) {
			//if (j==151) {
			//	cout << endl;
			//}
			var_node_proc(j); //processing the j-th var node
		}
		//At this point edge_llrs_cols should be updated

		for (j=0;j<nRows;j++) {
			//if ((j==26) || (j==29)) {
			//	cout << endl;
			//}
			check_node_proc(j,tx_cw_coset[j]);
		}
		//At this point edge_llrs_rows should be updated
		//cout << "Iteration number:" << i << endl;
		for (j=0;j<nCols;j++) {
			hardDecision[j] = final_decision_proc(j);
		}
		//cout << endl << endl << endl;
		synd_count = compute_syndrome(hardDecision,dec_coset);

		tmp = 0;
		for (j=0;j<nRows;j++) {
			if (dec_coset[j] != tx_cw_coset[j]) {
				tmp++;
			}
		}
		if (tmp == 0)
			break;
	}

	for (i=0;i<nCols;i++)
		decodedCW[i] = hardDecision[i];

	//clean up
	delete hardDecision;
	delete tx_cw_coset;
	delete dec_coset;

	if (tmp == 0) {
		//cout << "Decoding success" << endl;
		return 1;
	}
	else {
		//cout << "Decoding Failure" << endl;
		return -1;
	}

}

//Collecting the near codewords assuming that the all-zeros codeword was transmitted
int LDPCDec::decode_nearCW(double* _channel_llrs, int maxIter, int *decodedCW, char* nearCW_file, char* small_err_file)
{
	int i;
	int j;
	int tmp=0;
	int nearCW_flag=0;
	int nearViol=0;
	int synd_count;
	int *hardDecision;
	int *prev_hardDecision;
	
	hardDecision = new int[nCols];
	prev_hardDecision = new int[nCols];

	for (i=0;i<nCols;i++)
		channel_llrs[i] = _channel_llrs[i];

	//Compute the hard decision corresponding to the bits
	for (i=0;i<nCols;i++){
		//Important change 09/03/07. Using the mapping
		// c_i -> 1 - 2*c_i.
		hardDecision[i] = (channel_llrs[i] < 0.0);
		//hardDecision[i] = (channel_llrs[i] > 0.0);
	}

	//Compute the syndrome based on the hardDecision
	synd_count = compute_syndrome(hardDecision);

	if (synd_count == 0) {

		for (i=0;i<nCols;i++) {
			decodedCW[i] = hardDecision[i];
		}
		//cout << "Decoding success" << endl;
		return 1;
	}

	reset_llrs(); //may have to change if we implement turbo-eq

	for (i=0;i<maxIter;i++) {

		//before beginning the 45th iteration record the hardDecision at the end of 44th iteration
		if (i >= NEARCW_START_ITER) {
			//start tracking the hard decision
			for (j=0;j<nCols;j++) {
				prev_hardDecision[j] = hardDecision[j];
				//if (hardDecision[j] == 1) {
				//	cout << j << " ";
				//}				
			}
			//cout << endl;
		}
 		
		for (j=0;j<nCols;j++) {
			var_node_proc(j); //processing the j-th var node
		}
		//At this point edge_llrs_cols should be updated

		for (j=0;j<nRows;j++) {
			check_node_proc(j);
		}
		//At this point edge_llrs_rows should be updated

		for (j=0;j<nCols;j++) {			
			hardDecision[j] = final_decision_proc(j);
		}

		if (i >= NEARCW_START_ITER) {
			for (j=0;j<nCols;j++) {
				if (hardDecision[j] != prev_hardDecision[j]) {
					nearViol++;
					//cout << j << " ";
				}
			}
			//cout << endl;
		}
		synd_count = compute_syndrome(hardDecision);

		if (synd_count == 0)
			break;
	}

	for (i=0;i<nCols;i++)
		decodedCW[i] = hardDecision[i];
	
	//clean up
	delete hardDecision;
	delete prev_hardDecision;

	nearCW_flag = (nearViol == 0) && (synd_count != 0);

	//small error flag	
	for (i=0;i<nCols;i++) {
		if (decodedCW[i] != 0)
			tmp++;
	}

	if ((tmp <= SMALL_ERR_SIZE) && (tmp > 0)) {
		ofstream fp_small_err;

		fp_small_err.open(small_err_file,ios::app);
		for (i=0;i<nCols;i++) {
			if (decodedCW[i] != 0) {				
				fp_small_err << i << " ";
			}
		}
		fp_small_err << endl;
		fp_small_err.close();
	}

	if (nearCW_flag == 1) {
		//Write the near codeword onto the file
		ofstream fp;

		fp.open(nearCW_file,ios::app);

		for (i=0;i<nCols;i++) {
			if (decodedCW[i] != 0) {
				fp << i << " ";
			}
		}
		fp << endl;
		fp.close();
	}

	return nearCW_flag;

}

//Collecting the near codewords assuming that the all-zeros codeword was transmitted
int LDPCDec::decode_nearCW(double* _channel_llrs, int maxIter, int *decodedCW, char* nearCW_file)
{
	int i;
	int j;
	int tmp=0;
	int nearCW_flag=0;
	int nearViol=0;
	int synd_count;
	int *hardDecision;
	int *prev_hardDecision;
	
	hardDecision = new int[nCols];
	prev_hardDecision = new int[nCols];

	for (i=0;i<nCols;i++)
		channel_llrs[i] = _channel_llrs[i];

	//Compute the hard decision corresponding to the bits
	for (i=0;i<nCols;i++){
		//Important change 09/03/07. Using the mapping
		// c_i -> 1 - 2*c_i.
		hardDecision[i] = (channel_llrs[i] < 0.0);
		//hardDecision[i] = (channel_llrs[i] > 0.0);
	}

	//Compute the syndrome based on the hardDecision
	synd_count = compute_syndrome(hardDecision);

	if (synd_count == 0) {

		for (i=0;i<nCols;i++) {
			decodedCW[i] = hardDecision[i];
		}
		//cout << "Decoding success" << endl;
		return 1;
	}

	reset_llrs(); //may have to change if we implement turbo-eq

	for (i=0;i<maxIter;i++) {

		//before beginning the 45th iteration record the hardDecision at the end of 44th iteration
		if (i >= NEARCW_START_ITER) {
			//start tracking the hard decision
			for (j=0;j<nCols;j++) {
				prev_hardDecision[j] = hardDecision[j];
				//if (hardDecision[j] == 1) {
				//	cout << j << " ";
				//}				
			}
			//cout << endl;
		}
 		
		for (j=0;j<nCols;j++) {
			var_node_proc(j); //processing the j-th var node
		}
		//At this point edge_llrs_cols should be updated

		for (j=0;j<nRows;j++) {
			check_node_proc(j);
		}
		//At this point edge_llrs_rows should be updated

		for (j=0;j<nCols;j++) {			
			hardDecision[j] = final_decision_proc(j);
		}

		if (i >= NEARCW_START_ITER) {
			for (j=0;j<nCols;j++) {
				if (hardDecision[j] != prev_hardDecision[j]) {
					nearViol++;
					//cout << j << " ";
				}
			}
			//cout << endl;
		}
		synd_count = compute_syndrome(hardDecision);

		if (synd_count == 0)
			break;
	}

	for (i=0;i<nCols;i++)
		decodedCW[i] = hardDecision[i];
	
	//clean up
	delete hardDecision;
	delete prev_hardDecision;

	nearCW_flag = (nearViol == 0) && (synd_count != 0);

	if (nearCW_flag == 1) {
		//Write the near codeword onto the file
		ofstream fp;

		fp.open(nearCW_file,ios::app);

		for (i=0;i<nCols;i++) {
			if (decodedCW[i] != 0) {
				fp << i << " ";
			}
		}
		fp << endl;
		fp.close();
	}

	return nearCW_flag;

}



int LDPCDec::coset_decode_nearCW(int *tx_cw,double* _channel_llrs, int maxIter, int *decodedCW, char* nearCW_file, char* small_err_file)
{
	int i;
	int j;
	int tmp=0;
	int nearCW_flag=0;
	int nearViol=0;
	int synd_count;
	int *hardDecision;
	int *prev_hardDecision;
	int *tx_cw_coset;
	int *dec_coset;
	
	hardDecision = new int[nCols];
	prev_hardDecision = new int[nCols];
	tx_cw_coset = new int[nRows];
	dec_coset = new int[nRows];

	//First compute the coset of the transmitted codeword
	compute_syndrome(tx_cw,tx_cw_coset);

	for (i=0;i<nCols;i++)
		channel_llrs[i] = _channel_llrs[i];

	//Compute the hard decision corresponding to the bits
	for (i=0;i<nCols;i++){
		//Important change 09/03/07. Using the mapping
		// c_i -> 1 - 2*c_i.
		hardDecision[i] = (channel_llrs[i] < 0.0);
		//hardDecision[i] = (channel_llrs[i] > 0.0);
	}

	//Compute the syndrome based on the hardDecision
	synd_count = compute_syndrome(hardDecision,dec_coset);

	tmp = 0;
	for (i=0;i<nRows;i++) {
		if (dec_coset[i] != tx_cw_coset[i]) {
			tmp++;
		}
	}
	
	if (tmp == 0) {
		for (i=0;i<nCols;i++) {
			decodedCW[i] = hardDecision[i];
		}//Decoding succeeded without iters
		return 1;
	}

	reset_llrs(); //may have to change if we implement turbo-eq

	for (i=0;i<maxIter;i++) {

		//if (i == 12)
		//	cout << endl;
	
		//before beginning the 45th iteration record the hardDecision at the end of 44th iteration
		if (i >= NEARCW_START_ITER) {
			//start tracking the hard decision
			for (j=0;j<nCols;j++) {
				prev_hardDecision[j] = hardDecision[j];
				//if (hardDecision[j] != tx_cw[j]) {
				//	cout << j << " ";
				//}				
			}
			//cout << endl;
		}
 		
		for (j=0;j<nCols;j++) {
			var_node_proc(j); //processing the j-th var node
		}
		//At this point edge_llrs_cols should be updated

		for (j=0;j<nRows;j++) {
			//if (j==69)
			//	cout << endl;

			check_node_proc(j,tx_cw_coset[j]);
		}
		//At this point edge_llrs_rows should be updated

		for (j=0;j<nCols;j++) {
			hardDecision[j] = final_decision_proc(j);
		}

		if (i >= NEARCW_START_ITER) {
			for (j=0;j<nCols;j++) {
				if (hardDecision[j] != prev_hardDecision[j]) {
					nearViol++;
					//cout << j << " ";
				}
			}
			//cout << endl;
		}
		synd_count = compute_syndrome(hardDecision,dec_coset);

		tmp = 0;
		for (j=0;j<nRows;j++) {
			if (dec_coset[j] != tx_cw_coset[j]) {
				tmp++;
			}
		}
		if (tmp == 0)
			break;

	}

	for (i=0;i<nCols;i++)
		decodedCW[i] = hardDecision[i];
	
	//clean up
	delete hardDecision;
	delete prev_hardDecision;
	delete tx_cw_coset;
	delete dec_coset;

	nearCW_flag = (nearViol == 0) && (tmp != 0);

	//small error flag	
	tmp = 0;
	for (i=0;i<nCols;i++) {
		if (decodedCW[i] != tx_cw[i])
			tmp++;
	}

	if ((tmp <= SMALL_ERR_SIZE) && (tmp > 0)) {
		ofstream fp_small_err;

		fp_small_err.open(small_err_file,ios::app);
		for (i=0;i<nCols;i++) {
			if (decodedCW[i] != tx_cw[i]) {				
				fp_small_err << i << " ";
			}
		}
		fp_small_err << endl;
		fp_small_err.close();
	}

	if (nearCW_flag == 1) {
		//Write the near codeword onto the file
		ofstream fp;

		fp.open(nearCW_file,ios::app);

		for (i=0;i<nCols;i++) {
			if (decodedCW[i] != tx_cw[i]) {
				fp << i << " ";
			}
		}
		fp << endl;
		fp.close();
	}

	return nearCW_flag;

}

int LDPCDec::coset_decode_nearCW(int *tx_cw,double* _channel_llrs, int maxIter, int *decodedCW, char* nearCW_file)
{
	int i;
	int j;
	int tmp=0;
	int nearCW_flag=0;
	int nearViol=0;
	int synd_count;
	int *hardDecision;
	int *prev_hardDecision;
	int *tx_cw_coset;
	int *dec_coset;
	
	hardDecision = new int[nCols];
	prev_hardDecision = new int[nCols];
	tx_cw_coset = new int[nRows];
	dec_coset = new int[nRows];

	//First compute the coset of the transmitted codeword
	compute_syndrome(tx_cw,tx_cw_coset);

	for (i=0;i<nCols;i++)
		channel_llrs[i] = _channel_llrs[i];

	//Compute the hard decision corresponding to the bits
	for (i=0;i<nCols;i++){
		//Important change 09/03/07. Using the mapping
		// c_i -> 1 - 2*c_i.
		hardDecision[i] = (channel_llrs[i] < 0.0);
		//hardDecision[i] = (channel_llrs[i] > 0.0);
	}

	//Compute the syndrome based on the hardDecision
	synd_count = compute_syndrome(hardDecision,dec_coset);

	tmp = 0;
	for (i=0;i<nRows;i++) {
		if (dec_coset[i] != tx_cw_coset[i]) {
			tmp++;
		}
	}
	
	if (tmp == 0) {
		for (i=0;i<nCols;i++) {
			decodedCW[i] = hardDecision[i];
		}//Decoding succeeded without iters
		return 1;
	}

	reset_llrs(); //may have to change if we implement turbo-eq

	for (i=0;i<maxIter;i++) {

		//if (i == 12)
		//	cout << endl;
	
		//before beginning the 45th iteration record the hardDecision at the end of 44th iteration
		if (i >= NEARCW_START_ITER) {
			//start tracking the hard decision
			for (j=0;j<nCols;j++) {
				prev_hardDecision[j] = hardDecision[j];
				//if (hardDecision[j] != tx_cw[j]) {
				//	cout << j << " ";
				//}				
			}
			//cout << endl;
		}
 		
		for (j=0;j<nCols;j++) {
			var_node_proc(j); //processing the j-th var node
		}
		//At this point edge_llrs_cols should be updated

		for (j=0;j<nRows;j++) {
			//if (j==69)
			//	cout << endl;

			check_node_proc(j,tx_cw_coset[j]);
		}
		//At this point edge_llrs_rows should be updated

		for (j=0;j<nCols;j++) {
			hardDecision[j] = final_decision_proc(j);
		}

		if (i >= NEARCW_START_ITER) {
			for (j=0;j<nCols;j++) {
				if (hardDecision[j] != prev_hardDecision[j]) {
					nearViol++;
					//cout << j << " ";
				}
			}
			//cout << endl;
		}
		synd_count = compute_syndrome(hardDecision,dec_coset);

		tmp = 0;
		for (j=0;j<nRows;j++) {
			if (dec_coset[j] != tx_cw_coset[j]) {
				tmp++;
			}
		}
		if (tmp == 0)
			break;

	}

	for (i=0;i<nCols;i++)
		decodedCW[i] = hardDecision[i];
	
	//clean up
	delete hardDecision;
	delete prev_hardDecision;
	delete tx_cw_coset;
	delete dec_coset;

	nearCW_flag = (nearViol == 0) && (tmp != 0);

	if (nearCW_flag == 1) {
		//Write the near codeword onto the file
		ofstream fp;

		fp.open(nearCW_file,ios::app);

		for (i=0;i<nCols;i++) {
			if (decodedCW[i] != tx_cw[i]) {
				fp << i << " ";
			}
		}
		fp << endl;
		fp.close();
	}

	return nearCW_flag;

}


LDPCDec::~LDPCDec(void)
{
	int i;
	delete HRowProfile;
	delete HColProfile;

	for (i=0;i<nRows;i++) 
		delete[] HRows[i];

	delete HRows;

	delete cum_check_edges;
	delete cum_var_edges;

	delete check_to_var_perm;
	delete var_to_check_perm;
	
	delete edge_llrs_rows;
	delete edge_llrs_cols;
}
