#include "Util.h"
#include <boost/lexical_cast.hpp>

namespace MTMCSim{
    namespace Util{
    bool getIntPara(const paraMap& para, const string& name, int& data)
    {
        paraMap::const_iterator it = para.find(name);
        if (it != para.end())
        {
            data = boost::lexical_cast<int>(it->second);
            return true;
        }

        return false;
    }

    bool getDouPara(const paraMap& para, const string& name, double& data)
    {
        paraMap::const_iterator it = para.find(name);
        if (it != para.end())
        {
            data = boost::lexical_cast<double>(it->second);
            return true;
        }
        return false;
    }

    bool getStrPara(const paraMap& para, const string& name, string& data)
    {
        paraMap::const_iterator it = para.find(name);
        if (it != para.end())
        {
            data = it->second;
            return true;
        }       
        return false;
    }
    
string getFormattedStr( const paraMap& inputPara, const compResMap& compRes, string formatFN)
{
    string oneline;
    string DispStr;

    ifstream formatStream(formatFN.c_str());

    if (!formatStream)
        throw runtime_error("Read File Error!" + formatFN);

    getline(formatStream, oneline);

    while (formatStream)
    {
        int symInd = 0, secInd = 0;
        int start = DispStr.length();
        int next = 0;

        DispStr += oneline;
        DispStr += "\n";

        symInd = DispStr.find_first_of("$@",0);

        while (symInd != string::npos)
        {
            
            secInd = DispStr.find_first_of("$@", symInd+1);

            string key = DispStr.substr(symInd+1,secInd - symInd - 1);

            bool flag = false;

            
            paraMap::const_iterator it = inputPara.find(key);
            
            if ((DispStr[symInd] == '$')&& 
                (it != inputPara.end()))

            {
                DispStr.replace(symInd, secInd - symInd + 1, it->second);
                next = symInd + it->second.length();
                flag = true;
            }

            compResMap::const_iterator it2 = compRes.find(key);

            if ((DispStr[symInd] == '@')&&
                (it2!=compRes.end()))
            {
                DispStr.replace(symInd,  secInd - symInd + 1, it2->second);
                next = symInd + it2->second.length();
                flag = true;
            }

            if (!flag) 
            {   cout << "No entry corresponding to key " << key << endl;
                break;
            }
            symInd = DispStr.find_first_of("$@",next); // even if next is out-of-bound, no exception will be thrown, returns -1
        }

        getline(formatStream, oneline);
    }// end of while

    return DispStr;
}//end of function

void getConf(const string& confFN, string& workingDir,string&relative_profileFN,string& tempFN,string& formatFN, int& queueSize, int& nThreads)
{
    //string workingDir,relative_profileFN, tempFN, formatFN;
    //int queueSize, nThreads;
    ifstream configFile(confFN.c_str());
    if (!configFile)
        throw runtime_error("Read File Error!" + confFN);
    string oneline;
    getline(configFile, oneline);
    workingDir = oneline.substr(oneline.find('=')+1,oneline.find_last_of("%\r") - (oneline.find('=') + 1));
    getline(configFile, oneline);
    relative_profileFN = oneline.substr(oneline.find('=')+1,oneline.find_last_of("%\r") - (oneline.find('=') + 1));
    getline(configFile, oneline);
    tempFN = oneline.substr(oneline.find('=')+1,oneline.find_last_of("%\r") - (oneline.find('=') + 1));
    getline(configFile, oneline);
    formatFN = oneline.substr(oneline.find('=')+1,oneline.find_last_of("%\r") - (oneline.find('=') + 1));
    string intstr;
    getline(configFile,oneline);
    intstr = oneline.substr(oneline.find('=')+1,oneline.find_last_of("%\r") - (oneline.find('=') + 1));
    queueSize = atoi(intstr.c_str());
    getline(configFile,oneline);
    intstr = oneline.substr(oneline.find('=')+1,oneline.find_last_of("%\r") - (oneline.find('=') + 1));
    nThreads = atoi(intstr.c_str());
}

    }// end namespace util
}//end namespace MTMCSim



