/** \mainpage Multithreaded Monte Carlo Simulation platform. (MTMCSim)
*
*   \section author Author and contact information
*   Author: Shizheng Li

    Contact: szli.code AT gmail

    This project was mostly done when he was a PhD student at Dept. of Electrical and Computer Engineering, Iowa State University. 
    
    Most development was done between 2009 and 2010. Revisited in 2013 for a build system rework and making it open source.   
*     
*   \section intro Introduction
*   This project indends to have a generic platform for mulththreaded Monte Carlo simulations. As long as the simulation is written
*   according to the interface (abstract class MTMCSim::MTMCSimBase), the multithreading will be handled behind the scene. The multithreading
*   scheduler runs under boss-worker mode. The boss thread generates random number sequences and the worker threads performs
*   specified simulations. In current version, there is only one random number generator (RNG), thus one boss thread. We may have parallel RNG in the 
future. Therefore, it is best suited for simulations where the RN generation is simple compared to the algorithm run at each instance of the simulation. 
The project is designed in a purely object oriented manner in C++. This makes the project easy to maintain and extend. Almost all pointers 
in this project are boost::shared_ptr. All containers are from C++ STL. 
*    In this project, we call the data that gives one sample output to be a frame. For example, in communications theory, a frame is a data frame to be transmitted and received. In Monte Carlo simulation, we simulation 
*    many frames to get many samples and average over them. This terminalogy is used in communication and coding theory. In other fields, the terminalogy could be different.
*   
*   This project also provides a friendly user interface. All the inputs are handled through text files and the user can 
*   change the display and save format easily. It is not hard to provide a Python interface given the current design. 
*
This documentation is brief. Only displays public members of each classes and documentations for some of the members. 

\section structure Main Structure
\subsection core Core Infrastructure Classes

MTMCSim::MTMCScheduler: Multithreading scheduler class. Handles thread creation, synchronization, buffering, etc. 

MTMCSim::MTMCSimBase: Base class of simulators. All simulator should be derived from this class. 

MTMCSim::SimuPara: Base class of simulation parameters. All simulaton parameters should be derived from this class. See details in the \ref ui section.

Factory pattern (MTMCSim::SimuFactory, MTMCSim::SimuParaFactory) is used to produce simulator objects and simulation parameter objects.

The main() function is in RumSimu.cpp.

\subsection simulators Current Implemented Simulators (Mainly for coding theory)
MTMCSim::DSCKVSim: Koetter-Vardy algorithm for distributed source coding.

MTMCSim::KVChanSim: Koetter-Vardy algorithm for channel coding. This is the simplest use of this simulation program.

MTMCSim::MulStageSim: Multistage LDPC codes for distributed source coding. Use rate adaptive LDPC codes, designed by David Varodayan et.al. (http://ivms.stanford.edu/~varodayan/ldpca.html)

MTMCSim::MSALSim: Multistage LDPC codes for distributed source coding. Use optimized codes for AWGN channels, designed using PEG algorithm. Not rate adaptive code, but better performance

MTMCSim::MulSourSim: Multisource distributed source coding based LDPC codes for linear equation correlation.

The above simulations work in a "FER" mode. We keep simulating frames until we see 100 error frames (configurable), and compute FER. If multi run mode is implemented, we change the code rate and **run new simulations** in order to get close to a code rate that gives FER < 10^-3. The code rate is the "error free" code rate.

For rate adaptive codes, we can also run simulation in a "Rate" mode. In one instance of Monte Carlo simulation, we change the code rate until the decoding is successful. This can be done in practice since the codes are rate-adaptive, just transmit more parity bits. We simulate 100 (configurable) frames and average the code rate, this is the average required code rate to achieve "error free" transmission.

MTMCSim::MulStageSimRate and MTMCSim::MulSourSimRate implements the "Rate" mode.

A simulation that implements MTMCSim::MultiRunSimuBase can run in "multi run" mode. It allows you to change your parameter and run the simulation again, until you set a quit flag to be true. Currently it is used to change the code rate heuristically to find rate that gives FER<10^-3.

For a complete list of class hierarchy, please see <a href="hierarchy.html">Class Hierarchy</a> .

\subsection utils Utilility classes
Here are some classes and functionarities that may be reused by other research projects.

MTMCSim::pdfXY Joint probablity densigy (mass) function for two non-binary sources. It can be qary symmetric sources, or pdf (pmf more procisely) given by a file. It computes entropy, conditional entropy, mutual information, marginal pdf, etc.

MTMCSim::discrete_distribution Random number generator for 1D general discrete probability distributions. It follows boost::random interface convention. It takes in a vector of double that gives the pmf. The algorithm is "slow start and fast generate". At startup time, it creates two tables for generating RNs, it is a bit slow, but when actual generating discrete random RVs, it is very fast. (So it is different from the "interval look up"). The table can be cached and reused, further reduces the start up time after one distribution is used once.

MTMCSim::nbToBinPDFConv Converts a pmf on alphabet with size nbSize to log_2(nbSize) binary sources with a joint pdf. Also computes conditional entropies and conditional probabilities condition on K number of bits where K between 0 to log_2(nbSize).

MTMCSim::LDPCDec General LDPC decoder using belief-propgation algorithm written by my advisor Aditya Ramamoorthy. It supports ordinary LDPC decoding (find a code word) or decode to syndrome (finding a vector in a coset corresponding to the syndrome), which is useful for distributed source coding.


\section sys Build and Environment

The project is portable between Windows and Linux. Tested under Windows 7 32-bit and 64-bit, Visual Studio 2010 and Fedora 13 and Centos 6.3 with Gcc 4.4.
   
The project requires Boost library. No strict requirement on versions. Tested under 1.43/1.44/1.47 version. Boost::thread must be built on the system.  Boost shared_ptr is heavily used.

We use CMake as the build tool. This allow cross-platform builds. You can generate Visual Studio solution/project files, or make files, or nmake files, etc.

To obtain high performance, the hardware should offer good parallel compute capability. Typically the number of working threads equals to the number of CPU cores.

*   \section ui User Interface

The user interface is written in a pure OO manner and intended to potentially support  GUI or Python in the future.
\subsection input Input files
*   The input parameters are given in text files. Typically, there are two input files. One is the top level configuration file,
*   the file name is MTMCSimconf.txt by default, or is the input paramater of the exetuable file. A typical top configuration file looks like as follows.
*
*<em>
*   simulation name = MulStageSim   % Simulation name tells which simulator to load \n
*   workingDir=./ \n
*   profile file name=MuStSimuNew/QARY256_1.txt  % The file name of the simulation configuration..Relative path from working Dir. \n
*   temp file name=MuStSimuNew/MulStagetempFile.txt % Temporary resultion file name. Relative path from working Dir .\n
*   format file name=MuStSimuNew/MuStFormat.txt %Output format file name.Relative path from working Dir. \n
*   queueSize=8 % The size of the buffer queue between boss thread and working thread. \n
*   nThreads=4 %The number of working threads.\n
*   multiple run mode = 0   % 1 if automatically run simulatin using different parameters, the simulation should be MTMCSim::MultiRunnable.\n
*   no simulation = 0 % If this is 1, we do not run the simulation. Only displays the parameters and computed results. \n
*</em>
*   In general, the input files in this project have the form like above. Each effecive line has the format Parameter name = Parameter value,
*   and the space or tab before and after the "=" will be ignored. Any line without "=" will be comment line. In an effective line, any string after % 
*   will be comments. 
*
*   Another input file gives detailed configuration of a specific simulation. Its file name is specified by "profile file name". 
*   Here is an example.
*<em>
*   QARY symmetric correlation. KV for Distributed source coding. 

alphabet size = 256 \n
code length = 255\n
info length = 65\n
pdf para = #QARY 0.4\n
seed  = 50\n
lambda = 100.99\n

max error frame = 100 % Max # of error frame to be simulated, termination condition. \n
display freq = 10 % Display temporary results every "display freq" frames. \n
save freq = 500 % Save temporary results every "save freq" frames.\n
save file = DSCKVSim/QARY256Res_Debug. % Result file name. txt\n
multiple run mode = 1\n

</em>

The names of the parameters and meanings of the parameters are determined by specific simulator class. 

\subsection out Output files
The output file format is given by a format file. The results are categorized into Simulation Results and Computed Results. Computed Results
are computed from the parameters, e.g., the source entropy, the channel capacity, etc. The format file specifies the format to display input
parameters and the computed results. Here is one line of the format file. 

Alphabet Size $alphabet size$, PDF Type $pdf para$. 
HX = @ HX @, HY = @ HY @, H(Y|X) = @ HYX @.

The parameter names are in between the $ symbol and the computed result names are in between the @ symbol. In the result file, the strings containing
$ and @ will be replaced by the actual parameters or computed results. The user can change the display format and result file format easily. In the current
version, the display and result file formats are the same. But it is easy to modified the code so that they are different. 

\subsection under_the_hood Parameter handling under the hood
This part gives some implementation details about the user interface. This is useful for people who want to write the simulator by themselves. 

The parameters for a specific simulator are saved in a nested class of the simulator, which is derived from MTMCSim::SimuPara. The parameters are stored in public data members of the specific parameter class. The input text file is read by MTMCSim::TextFileInput 
class and converted the an intermediate map (std::map<string,string>) that maps the parameter name to the parameter value. Then, a specifed version of virtual function MTMCSim::SimuPara::loadPara() 
should be called (usually called in a polymorphism manner) and the parameters will be converted to values according to their types and save to the data members of a specified parameter class.
Each specified parameter class should define its own MTMCSim::SimuPara::loadPara() realization.

The computed results for a specific simulator are saved in a nested class of the simulator, which is derived from MTMCSim::CompRes. The MTMCSim::CompRes class defines a 
virtual function MTMCSim::CompRes::setMapping() to convert the results in data members to std::map<string,string>. MTMCSim::Util::getFormattedStr() will convert the parameters and 
computed results into string using given format.

*   \section howto How to write a simulation
*   To write a specified simulation, a simulator class should be derived from  MTMCSim::MTMCSimBase class. See the doc of MTMCSim::MTMCSimBase
*   for more information. 

A good example of such a simulator  with detailed comments is MTMCSim::KVChanSim and MTMCSim::MulStageSim.

*   
*   
*/