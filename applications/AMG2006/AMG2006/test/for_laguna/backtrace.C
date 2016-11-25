/*
 * backtrace.C
 *
 *  Created on: Dec 12, 2011
 *      Author: Ignacio Laguna
 *     Contact: ilaguna@purdue.edu
 */

#if HAVE_CONFIG_H
#include "statetracker-config.h"
#endif

#include "backtrace.h"
#include "debugging.h"
#include "mpi_state.h"
#include "utilities.h"
#include "io_utilities.h"

#include <execinfo.h>
#include <unistd.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

string Backtrace::add2lineProgram = "/usr/bin/addr2line";
string Backtrace::objdumpProgram = "/usr/bin/objdump";

string Backtrace::getBacktrace()
{
	string ret = getGNUBacktrace();
	return ret;
}

string Backtrace::getGNUBacktrace()
{
	int maxStackSize = 100;
	void *array[maxStackSize];
	size_t size;
	char **names;
 	int i;

	size = backtrace(array, maxStackSize);
        for (i=0; i < size-1; ++i) {
   		array[i] -= 1;
        } 
	names = backtrace_symbols(array, size);

	// Get string of stack-frames
	string frames("");
	char p[500];
	if (size != 0) {
		//for (size_t i = 0; i < size; ++i) {
		for (size_t i = 0; i < size-1; ++i) { // don't take the last frame
			/* Print return addresses */
			//p[0] = '\0';
                        //sprintf(p, "%p", array[i]);
			//frames = frames + "|" + p;

			/* Print function name, offset, and return addr. */
			frames = frames + string("|") + names[i];

			/* Print function name and offset only */
			//vector<string> tokens;
			//Tokenize(names[i], tokens, " ");
			//frames = frames + string("|") + tokens[0];
		}
	} else {
#if STATE_TRACKER_DEBUG
    		handleError("in Backtrace::getGNUBacktrace(): "
    			"insufficient memory for the strings of backtrace");
#endif
	}

	return frames;
}

string Backtrace::convertAddrToLine(const string &addr, const string &program)
{
	string ret("");
	string command = add2lineProgram + " -e " + program + " " + addr;
	vector<string> output = executeShellCommand(command);
        return output[0];
}

bool Backtrace::canConvertAddrToLine()
{
	bool ret = true;
	// check that addr2line utility exists
	if (access(add2lineProgram.c_str(), F_OK) != 0) {
		ret = false;
		string msg = "addr2line program is not accessible in " +
				add2lineProgram;
		handleError(msg.c_str());
	}
	return ret;
}

bool Backtrace::canFindFileAndFunctionFromObject()
{
	bool ret = true;
	// check that addr2line utility exists
	if (access(objdumpProgram.c_str(), F_OK) != 0) {
		ret = false;
		string msg = "objdump program is not accessible in " +
				objdumpProgram;
		handleError(msg.c_str());
	}
	return ret;
}

void separateInformation(
		const string &functOffset,
		string &module,
		string &func,
		string &offset
		)
{
	vector<string> tokens;
	Tokenize(functOffset, tokens, " ");
	string tmp = tokens[0]; //takes the function+offset part

	tokens.clear();
	Tokenize(tmp, tokens, "(");
	module = tokens[0];

	tmp = tokens[1];
	tokens.clear();
	Tokenize(tmp, tokens, "+");
	func = tokens[0];
	offset = tokens[1].substr(0, tokens[1].size()-1);
}

/*
 * Input corresponds to function offset as given by GNU backtrace API:
 * /path/to/module(FUNCTION+OFFSET) [address]
 */
FileAndFunction Backtrace::findFileAndFunctionFromObject(const string &functOffset)
{
   FileAndFunction ret;
	string module, func, offset;
	separateInformation(functOffset, module, func, offset);

	// Get base address for function
	string baseAddr("");
	string command = objdumpProgram + " -t " + module + " | /bin/grep -e \" " + func + "$\"";
	vector<string> output = executeShellCommand(command);
	if(output.size() == 0)
	{
		cout << "Error in command:=>" << command << "=>return empty vector" << endl; 
		return ret;
	}
	baseAddr += output[0];

	// Get only the address from 'objdump'
	vector<string> tokens;
	Tokenize(baseAddr, tokens, " ");
	baseAddr = tokens[0];
	if(tokens.size() == 0)
	{
		 cout << "Error in baseAddr:=>" << baseAddr << "=>return empty vector" << endl;
		return ret;	
	}

	// Convert base and offset from string to integers
	unsigned long baseValue = convertStringHexToInteger(baseAddr);
	unsigned long offsetValue = convertStringHexToInteger(offset);

	// Convert back to string and get actual address
	char buff[128];
	sprintf(buff, "%x", (baseValue + offsetValue));
	string addr(buff);
	string line = convertAddrToLine(addr, module);

	ret.fileNameAndLine = line;
	ret.functionName = func;
	if (module.find("libstracker.") != string::npos)
		ret.fromTool = true;
	else
		ret.fromTool = false;
	
	return ret;
}
