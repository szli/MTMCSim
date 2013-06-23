#ifndef _DOSDISPUI_H_INCLUDED
#define _DOSDISPUI_H_INCLUDED
#include "StdHeaders.h"
#include "DispUI.h"

namespace MTMCSim{
class DOSDispUI: public DispUI
{
public:
	virtual void dispStr(const string& str) 
	{
		cout << str;
		cout.flush();
	}

};
}
#endif
