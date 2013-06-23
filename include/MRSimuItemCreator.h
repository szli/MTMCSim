#ifndef _MRSIMUITEMCREATOR_H_INCLUDED
#define _MRSIMUITEMCREATOR_H_INCLUDED

#include "MRSimuFactory.h"
namespace MTMCSim{

/** class template to create objects that are derived classes of MultiRunSimuBase. T is the actual (derived) simulator type 
	and U is the actual (derived) simulation parameter class. 
	*/
template <class T, class U>
class MRSimuItemCreator{
public:
	static shared_ptr<MultiRunSimuBase> create(SimuPara& para, TextSaveUI& resFile_, TextSaveUI& tempFile_, DispUI& dispUI_)
	{
		return shared_ptr<T>(new T(dynamic_cast<U&>(para), resFile_,  tempFile_,dispUI_));
	}

	MRSimuItemCreator(string name)
	{
		(MRSimuFactory::instance()).regSimulator(name,MRSimuItemCreator<T,U>::create);
	}



};
}
#endif