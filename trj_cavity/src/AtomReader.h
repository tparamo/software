/*
 * AtomReader.h
 *
 *  Created on: Nov 17, 2010
 *      Author: tp334
 */

#ifndef ATOMREADER_H_
#define ATOMREADER_H_

#include <iostream>
#include <vector>

#include "Atom.h"
#include "Grid.h"
#include "Coordinates.h"

#define FACTOR 10

using namespace std;


class AtomReader {

	protected:
		string fpath;

	public:
		AtomReader();
		AtomReader(string fpath);
		FILE* createFile(string fpath, string fname, string mode);
		string getFPath();
        //Grid readPDB(FILE* pdbfile, float gridspacing, bool gromos);
		void calculateLimitsGrid(float *limits, Atom a);
		Grid initializeGrid(Grid grid,vector<Atom> atoms);
		int* distanceCells(float spacing, float radius, int center, int limit);
		//void readTrajectory(char *inputfile, char *topology, float gridspacing, FILE* cavfile, const char *trjfile, FILE* volumefile, FILE* statfile, int mode, int tstat, bool gromos,map<string,map<string,double> > forceFieldRadii, FILE* watersfile, FILE* water_statistics);
		void start(int argc,char *argv[]);
		virtual ~AtomReader();
};

#endif /* ATOMREADER_H_ */
