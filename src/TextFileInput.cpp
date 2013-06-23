#include "TextFileInput.h"

namespace MTMCSim{

paraMap TextFileInput::getParaMap(){
		
	ifstream inFile(inFN.c_str());
	if (!inFile)
		throw runtime_error("Read File Error!" + inFN);
	paraMap inputPara;

	string oneline;

	getline(inFile, oneline);
	while (inFile)
	{
		
		int eqInd = oneline.find('=');
		if (eqInd == string::npos) // if '=' is not found 
		{ 
			getline(inFile, oneline);
			continue;
		}

		string key = oneline.substr(0, eqInd);
		string value = oneline.substr(eqInd+1, oneline.find_first_of("%\r") - (eqInd + 1)); // changed to find_first_of from find_last_of on Sep 21, 2010. 

		// remove spaces/tabs at the begining and end of key/value
		if (!key.empty())
		{
		while ((*(key.begin()) == ' ')||(*(key.begin()) == '\t'))
			value.erase(key.begin());
		while ((*(key.end()-1) == ' ')||(*(key.end()-1) == '\t'))
			key.erase(key.end()-1);
		}

		if (!value.empty())
		{

		while ((*(value.begin()) == ' ')||(*(value.begin()) == '\t'))
			value.erase(value.begin());
		while ((*(value.end()-1) == ' ')||(*(value.end()-1) == '\t'))
			value.erase(value.end()-1);
		while ((*(value.end()-1) == '\r'))
			value.erase(value.end()-1);

		}

		inputPara[key] = value;

		getline(inFile, oneline);
	}// end while

	return inputPara;
}
}