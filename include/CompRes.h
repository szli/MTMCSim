#ifndef _COMPRES_H_INCLUDED
#define _COMPRES_H_INCLUDED
#include "StdHeaders.h"

namespace MTMCSim{

/** Base class for computed results.
*/
class CompRes { 

	/** Transfer the data (results) stored in the object to a mapping between the variable name and variable value (both string).
	*/
	virtual void setMapping(compResMap& crMap)const  = 0;
public:
	virtual ~CompRes()
	{}
};
}
#endif