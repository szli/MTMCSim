#ifndef _TEXTFILEINPUT_H_INCLUDED
#define _TEXTFILEINPUT_H_INCLUDED
#include "StdHeaders.h"

namespace MTMCSim{
/** Process input from a text file and convert to a parameter map (map<string, string>). */
class TextFileInput{
public:
	TextFileInput(std::string inFN_):inFN(inFN_) {}

	paraMap getParaMap();
private:
	string inFN;
};


}

#endif