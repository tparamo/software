/*
 * Log.cpp
 *
 *  Created on: Jan 13, 2011
 *      Author: tp334
 */

#include "Log.h"
#include "stdio.h"
#include <iostream>
#include <time.h>
#include "stdlib.h"

Log::Log() {
	// TODO Auto-generated constructor stub
}

Log::Log(string p) {
	this->setPath(p);
}

void Log::setPath(string p){
	path = p + "/Log.txt";
	//Erase previous content
	FILE *logfile = fopen(path.c_str(),"w");
	fclose(logfile);
}

string Log::getPath(){
	return path;
}

void Log::error(string message){
	FILE *logfile = fopen(path.c_str(),"a");
	message = getTime()+"  ERROR    : "+ message +"\n";
	fputs(message.c_str(),logfile);
	fclose(logfile);
}

void Log::info(string message){
	message = getTime()+"  INFO     : "+ message +"\n";
	FILE *logfile = fopen(path.c_str(),"a");
	fputs(message.c_str(),logfile);
	fclose(logfile);
}

void Log::warning(string message){
	FILE *logfile = fopen(path.c_str(),"a");
	message = getTime()+"  WARNING  : "+ message +"\n";
	fputs(message.c_str(),logfile);
	fclose(logfile);
}

string Log::getTime(){
	time_t rawtime;

	time (&rawtime);
	timeinfo = localtime(&rawtime);
	string ts = (string)asctime(timeinfo);
	int ch = ts.find_last_of("\n");
	ts.erase(ch, ch+1);

	return ts;
}

Log::~Log() {
}
