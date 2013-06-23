#ifndef _MTSIM_H_INCLUDED
#define _MTSIM_H_INCLUDED

#include "StdHeaders.h"
#include "MCThreadInfo.h"
#include "RNObj.h"
#include <boost/shared_ptr.hpp>

using boost::shared_ptr;

namespace MTMCSim{

/** The (abstract) base class for a multithread supported simulator. 
*/
class MTMCSimBase
{
public:

	/** This function will be called by the scheduler before the threads are created. \  Typically do resource allocation here.
	*/
	virtual void prepare(const MCThreadInfo&)
	{}

	/** The work to be done by a worker thread. This function should be written in a thread-safe manner, i.e., no conflict. no critical section. This typically
	* can be done by creating workspace for each thread, either as an automatic variable in the function, or as a data member in the simulator class. */
	virtual void oneThread(shared_ptr<RNObj>, const MCThreadInfo&) = 0;

	/** Returns the random number object, usually they are random sequences. Each object is for one thread. The object is created by "new" inside the 
	function and can be accessed througth the renturned shared_ptr. */
	virtual shared_ptr<RNObj> genRNObj() = 0;//generate randome seq for EACH THREAD

	
	/** Put the result from one thread to the final result object.
	*	It itself does not need to consider thread safty because a mutex will be locked before it is called. */
	virtual void dumpOneResult(const MCThreadInfo&) = 0;

	/** Return true if the terminate condition satisfies. Thread safty is partially guaranteed by the scheduler. The result object will 
	*	be locked but only the result object will be locked. So, typically it only reads the parameter object (which is always safe) and reads the result
	*	object (because the result object is locked by the scheduler)
	*/
	virtual bool checkTerm()  = 0;

	/** This function will be called after the RNG thread and worker threads are done. */
	virtual void postProcess()
	{}

	/** Single thread simulation. */
	void singleThreadSim()
	{
		MCThreadInfo threadInfo(1,0);

		prepare(threadInfo);

		while (!checkTerm())
		{
			oneThread(genRNObj(), threadInfo);
			dumpOneResult(threadInfo);
		}
		postProcess();
	}

	/** Calls compRes.setMapping(compResM_);.*/
	virtual void setCompResMap(compResMap& compResM_) const = 0;

	/** Virtual Destructor. */
	virtual ~MTMCSimBase()
	{}

};


}

#endif