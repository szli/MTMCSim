#ifndef _SAVEUI_H_INCLUDED
#define _SAVEUI_H_INCLUDED
#include "OutputUI.h"
#include "StdHeaders.h"

namespace MTMCSim{

/** An output UI that saves to a text file. */
class TextSaveUI: public OutputUI
{
public: 
	TextSaveUI(string saveFN_):saveStream(saveFN_.c_str(), ofstream::app)
	{}
		

	void saveStr(const string& str) 
	{
		saveStream << str;
		saveStream.flush();
	}
		

private:
	ofstream saveStream;

};
}
#endif