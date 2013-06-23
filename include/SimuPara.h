#ifndef _SIMUPARA_H_INCLUDED
#define _SIMUPARA_H_INCLUDED
#include "StdHeaders.h"
#include "DispUI.h"
#include <boost/shared_ptr.hpp>

using boost::shared_ptr;

namespace MTMCSim{

/** Base class for simulation paramaters. */
class SimuPara { 
	public:
	/** Load the parameters from a mapping, which is obtained from some form of user interface. The user interface type is irrelevent of the parameter
	*	object. 
	*	@param para a <string,string> (STL) mapping that contains parameters. 
	*	@param ui display user interface to display error message. 
	*/
	virtual void loadPara(const paraMap& para, DispUI& ui) = 0;

	virtual shared_ptr<SimuPara> clone() const = 0;

	string saveFN;///< The file  name of the result file. 

	virtual ~SimuPara()
	{}

};
}
#endif