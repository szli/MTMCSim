#ifndef _MULTIRUNABLE_H_INCLUDED
#define _MULTIRUNABLE_H_INCLUDED

#include <boost/shared_ptr.hpp>
#include "SimuPara.h"
using boost::shared_ptr;

namespace MTMCSim{
/** The interface (abstract base class) for multi-runnable simulators. In such a simulator, the simulator can update the parameters 
	for the next simulation and multiple simulations can run without human control. */
class MultiRunnable{
public:
	/** Get the parameter for the next simulation. */
	virtual shared_ptr<SimuPara> getNextPara()
	{ return nextPara;}

	/** Check if terminate the whole simulation series. */
	virtual bool quit() const
	{ return isQuit;}

	/** Default constructor. \ Not in the multi run mode by default. 
	*/
	MultiRunnable(): mulMode(0), isQuit(false) // Not in multirun mode
	{}
	
	
	MultiRunnable(int mode_, SimuPara& simupara): mulMode(mode_),isQuit(false)
	{
		if (mulMode != 0)
			nextPara = simupara.clone();
	}

	virtual ~MultiRunnable()
	{}

protected:
	int mulMode;
	bool isQuit;
	shared_ptr<SimuPara> nextPara;
};


}


#endif