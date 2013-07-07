#ifndef _DISPUI_H_INCLUDED
#define _DISPUI_H_INCLUDED
#include "StdHeaders.h"
#include "OutputUI.h"

namespace MTMCSim{
/** For displaying infomation on screen.
*	The derived class could be DOS prompt or GUI.
*/
class DispUI:public OutputUI
{
public:
	virtual void dispStr(const string& str) = 0;

	virtual ~DispUI()
	{}
};

}
#endif

