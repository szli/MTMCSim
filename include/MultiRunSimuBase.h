#ifndef _MULTIRUNSIMUBASE_H_INCLUDED
#define _MULTIRUNSIMUBASE_H_INCLUDED
#include "MTMCSimBase.h"
#include "MultiRunnable.h"

namespace MTMCSim{
/** The base class for the simulators that are multi-runnable and support our multithread interface. 
*/
class MultiRunSimuBase: public MTMCSimBase, public MultiRunnable
{
public:
	MultiRunSimuBase(int mode_, SimuPara& simupara):MultiRunnable(mode_,simupara)
	{}

	virtual ~MultiRunSimuBase()
	{}
};
}

#endif