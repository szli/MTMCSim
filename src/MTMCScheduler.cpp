#include "MTMCScheduler.h"

using boost::shared_ptr;
using boost::thread;
using boost::unique_lock;
using boost::shared_lock;
using boost::mutex;
using boost::shared_mutex;

namespace MTMCSim{

MTMCScheduler::MTMCScheduler(unsigned int nRunThreads_, unsigned int RNQSize_, shared_ptr<MTMCSimBase> simulator)
: nRunThreads(nRunThreads_), RNQSize(RNQSize_), MCRunThreads(nRunThreads_),
 MTMCSimulator(simulator),isTerminated(false)
{}

void MTMCScheduler::genRN()
{
    try
    {
        while (true)
        {
            shared_ptr<RNObj> pRN = MTMCSimulator->genRNObj();// should be be done with in the lock 
            { 
                unique_lock<mutex> lockQueueMut(queueMutex); // lock
                while (RNQueue.size() >= RNQSize) //queue is full
                {
                    CVQueueFull.wait(lockQueueMut);
                }

                RNQueue.push(pRN);
                
                if (RNQueue.size() == 1)
                    CVQueueEmpty.notify_one();
            }//release the lock

            {
                shared_lock<shared_mutex>  lockTerm(termMutex);
                if (isTerminated)
                {
                    break;
                }
            }//release the lock

        }// end while
    }//end try
    catch(const boost::thread_interrupted & ){ //Using interruption is one way to wake up the threads and terminate
            
            unique_lock<mutex> lockCoutMut(coutMutex);
            //cout << "genUniRN Thread:" << this_thread::get_id() << "Interrupted:" << endl;
            
        }

}



void MTMCScheduler::oneThreadSimu(MCThreadInfo threadInfo)
{
    shared_ptr<RNObj> currRN;

    try
    {

    while (true)
    {
        {
            unique_lock<mutex> lockQueueMut(queueMutex);

            while (RNQueue.empty())
            {
                CVQueueEmpty.wait(lockQueueMut);
            }

            currRN = RNQueue.front();
            RNQueue.pop();

            if (RNQueue.size() == RNQSize - 1)
                CVQueueFull.notify_one();
        }//release the lock


        MTMCSimulator->oneThread(currRN, threadInfo);

    //  {
    //      shared_lock<shared_mutex>  lockTerm(termMutex);
    //      if (isTerminated)// If the simulation is terminated, no need to dump the result
    //              break;
    //  }


        {
            unique_lock<mutex> lockResultMut(resultMutex);
            // If the simulation is terminated, no need to dump the result, it may cause error or inaccuracy
            // No need to lock shared_lock because only one thread will write isTerminated and it is this thread
            //E.g. Thread 1 dump result and set isTerminated = true, at the same time, Thread 2 is waiting for the resultMutex
            // after Thread 1 release the mutex, Thread 2 enters here and should get out immediately. 
            if (isTerminated)
                break;

            MTMCSimulator->dumpOneResult(threadInfo);
        

        if (MTMCSimulator->checkTerm())// if the terminal condition satisfies
        {
            unique_lock<shared_mutex> lockWriteTerm(termMutex);
            isTerminated = true;

            //interrupt all threads
            RNGThread->interrupt(); 
            for (unsigned int i = 0; i<nRunThreads; i++)
                MCRunThreads[i]->interrupt(); //Does not stop all other threads immediately because only if a thread is at interruption points, interupt takes effect. See boost thread doc. Note, lock is not interruption point.
            break;

        }
        }// release the result lock, it is released here because checkTerm() visits the member in result object, which may be changed
        // by other threads if not locked. 

        { // still need this, perhaps some threads will not go into 
        //"if (MTMCSimulator->checkTerm())" because it is false, but still needs exiting. Happens when some other thread
        // obtains resultMutex quickly and the termination condition satisfies there.

                    shared_lock<shared_mutex>  lockTerm(termMutex);
                    if (isTerminated)
                        break;
                }//release the lock

    }//end while
    }//end try
    catch(const boost::thread_interrupted & ){
                unique_lock<mutex> lockCoutMut(coutMutex);
            //cout << "solve Thread:" << this_thread::get_id() << "Interrupted:" << endl;
    }
}

void MTMCScheduler::start()
{
    isTerminated = false;
    MCThreadInfo initInfo(nRunThreads, -1);
    MTMCSimulator->prepare(initInfo);

    RNGThread.reset(new thread(boost::bind(std::mem_fun<void, MTMCScheduler>(&MTMCScheduler::genRN),this)));
    

    

    for (unsigned int i = 0; i<nRunThreads; i++)
    {
        MCThreadInfo threadInfo(nRunThreads,i);
        MCRunThreads[i].reset(new thread
            (boost::bind(std::mem_fun<void, MTMCScheduler,MCThreadInfo>(&MTMCScheduler::oneThreadSimu),
            this, threadInfo)));
    }

    RNGThread->join();
    for (unsigned int i = 0; i<nRunThreads;i++)
        MCRunThreads[i]->join();

    MTMCSimulator->postProcess();
}














}
