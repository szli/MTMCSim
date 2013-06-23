#ifndef _SIMUITEMCREATOR_H_INCLUDED
#define _SIMUITEMCREATOR_H_INCLUDED

#include "SimuFactory.h"
namespace MTMCSim{
/** class template to create objects that are derived classes of MTMCSimBase. T is the actual (derived) simulator type 
	and U is the actual (derived) simulation parameter class. 
	*/
template <class T, class U>
class SimuItemCreator{
public:
	static shared_ptr<MTMCSimBase> create(SimuPara& para, TextSaveUI& resFile_, TextSaveUI& tempFile_, DispUI& dispUI_)
	{
		return shared_ptr<T>(new T(dynamic_cast<U&>(para), resFile_,  tempFile_,dispUI_));
	}

	SimuItemCreator(string name)
	{
		(SimuFactory::instance()).regSimulator(name,SimuItemCreator<T,U>::create);
	}



};
}
#endif