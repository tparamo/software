/*
 * Log.h
 *
 *  Created on: Jan 13, 2011
 *      Author: tp334
 */

#ifndef LOG_H_
#define LOG_H_

#include <string>
#include "stdio.h"

using namespace std;

class Log {
	public:
		string path;
		struct tm * timeinfo;

	public:
		Log();
		Log(string p);
		void setPath(string path);
		string getPath();
		void error(string message);
		void info(string message);
		void warning(string message);
		string getTime();
		virtual ~Log();
};


#endif /* LOG_H_ */
