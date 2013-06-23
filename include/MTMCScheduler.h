#ifndef _MTMCSCHEDULER_H_INCLUDED
#define _MTMCSCHEDULER_H_INCLUDED

#include "MCThreadInfo.h"
#include "RNObj.h"
#include <boost/shared_ptr.hpp>
using boost::shared_ptr;
#include "MTMCSimBase.h"
#include <boost/thread.hpp>
#include <queue>
#include "StdHeaders.h"

using namespace boost;
using std::queue;



namespace MTMCSim{

/** Multithread Monte Carlo Simulator Boss-Worker model scheduler. 
*	In general, the Monte Carlo Simulation computes E[Y], Y = f(X), and X is the random sequence driving the simulation. 
*	The scheduler runs in the boss worker model, i.e., one thread generates 
*	random sequences X and several threads computes f(X) and obtain the samples of Y. Between the boss and the worker, there 
*	is a queue as a buffer of the random sequences. The simulator should be derived classes of the abstract calss MTMCSimBase
*
*	This model is suitable when the computation of the function f()
*	is high. Usually this is the case in coding theory. In particular, in coding theory, X is the random sequences, f() contains encoding, 
*	channel effect if there is any, decoding. Y = 1 if the frame is in error after decoding so E[Y] = probably of error. This is a top level description of
*	Monte Carlo. The X and Y here and X, Y in the multistage LDPC DSC are different.(There, they are source sequences and they corresponds to the  X here.
*/
class MTMCScheduler{
public:

	/**
	* @param nRunThreads_ The number of running threads (workers)
	* @param RNQSize The size of the random number queue
	* @param simulator The (smart) pointer to the simulator object
	*/
	MTMCScheduler(unsigned int nRunThreads_, unsigned int RNQSize,shared_ptr<MTMCSimBase> simulator );

	/** The default constructor does nothing. Should call setSimulator, setNRunThreads and setRNQSize before calling start().
	*/
	MTMCScheduler():isTerminated(false)
	{}

	/** Generates the random number sequences and put into the queue. */
	void genRN();

	/** Each worker thread call this function to work. */
	void oneThreadSimu(MCThreadInfo threadInfo);

	/** Start simulation.*/
	void start();

	void setSimulator(shared_ptr<MTMCSimBase> simulator_)
	{
		MTMCSimulator = simulator_;
	}

	void setNRunThreads(unsigned int nRunThreads_)
	{
		nRunThreads = nRunThreads_;
		MCRunThreads.resize(nRunThreads_);
	}

	void setRNQSize(unsigned int RNQSize_)
	{
		RNQSize = RNQSize_;
	}

private:
	unsigned int nRunThreads; ///< running threads, working threads, not the totoal number of threads

	unsigned int RNQSize;

	bool isTerminated;

	shared_ptr<MTMCSimBase> MTMCSimulator; // The scheduler object does not own simulator object, if we want it to own, need a wrapper.

	vector<shared_ptr<thread> >  MCRunThreads;

	queue<shared_ptr<RNObj> > RNQueue;

	shared_ptr<thread> RNGThread;


	mutex queueMutex;

	mutex resultMutex;

	condition_variable CVQueueFull;
	
	condition_variable CVQueueEmpty;

	shared_mutex termMutex;

	mutex coutMutex;

	




};



}





#endif