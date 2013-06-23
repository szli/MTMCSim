#ifndef _MCTHREADINFO_H_INCLUDED
#define _MCTHREADINFO_H_INCLUDED

namespace MTMCSim{
/** Thread info. Like a POD (plain old datetype).
*/
struct MCThreadInfo
{
	int nThread; ///< Number of solver runninging threads, does not include RNG thread.
	int threadID; ///< Thread ID in MTMCSim project, not the boost::thread::id type. \ Range from (0,1,..., nThread-1)
	MCThreadInfo(int nThread_, int threadID_)
		: nThread(nThread_), threadID(threadID_)
	{}
};
}
#endif