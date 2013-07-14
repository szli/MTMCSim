#ifndef _STDHEADERS_H_INCLUDED
#define _STDHEADERS_H_INCLUDED


#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <limits>
#include <map>
#include <algorithm>
#include <iterator>

using std::vector;
using std::string;
using std::ostream;
using std::istream;
using std::cout;
using std::cin;
using std::ofstream;
using std::ifstream;
using std::runtime_error;
using std::istringstream;
using std::getline;
using std::numeric_limits;
using std::streamsize;
using std::sort;
using std::endl;
using std::stringstream;
using std::map;
using std::ostringstream;

typedef map<string,string> paraMap;
typedef map<string,string> compResMap;

//extern const double EPS;
namespace MTMCSim{
	const double EPS = 1e-15;// changed from 1e-8 to 1e-15 on Sep 28, 2010
	const int MAX_VALUE = 1000000000;
}


#define SIM_ERR(err) \
do {                 \
  std::ostringstream  err_msg; \
  err_msg << "Error in: " << endl <<  \
          "Function: " << __FUNCTION__ << endl << \
          "File: " << __FILE__ << endl <<  \
          "Line: " << __LINE__ << endl << \
          "Message: " << err;   	   \
  throw std::runtime_error(err_msg.str()); \
 } while(0) 

#define SIM_ASSERT(cond, msg) if (cond); else SIM_ERR(msg)


#endif 
