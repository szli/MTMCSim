#ifndef _SIMUPARAITEMCREATOR_H_INCLUDED
#define _SIMUPARAITEMCREATOR_H_INCLUDED

#include "SimuParaFactory.h"
namespace MTMCSim{
/** class template to create objects that are derived classes of SimuPara. T is the actual (derived) simulation 
	 parameter class. 
	*/
template <class T>
class SimuParaItemCreator{
public:
	
	static shared_ptr<SimuPara> create()
	{
		return shared_ptr<T>(new T());
	}


	SimuParaItemCreator(string name)
	{
		(SimuParaFactory::instance()).regSimuPara(name,SimuParaItemCreator<T>::create);
	}

};
}
#endif