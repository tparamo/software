/*
 * AtomWriter.h
 *
 *  Created on: Dec 1, 2010
 *      Author: tp334
 */

#include "Coordinates.h"
#include "Grid.h"
#include <vector>
#include <string>
#include <map>

using namespace std;

#ifndef ATOMWRITER_H_
#define ATOMWRITER_H_

class AtomWriter {
	public:
		AtomWriter();
		void writePDB(vector<Coordinates> coordinates, FILE* file);
		void writePDB(vector<vector<Coordinates> > cavities, FILE* file);
		void writeStatistics(Grid statistics, float gridspacing, int nframes, FILE* file);
		void writeVolumeProbability(float volume[]);
		void writeVolumeFrame(int gridpoints, float spacing, float time, FILE* file);
		//void writeVolumeFrame(int gridpoints, vector<vector<Coordinates> > current, vector<vector<Coordinates> >& original, float spacing, int frame, FILE* file);
		vector<int> preprocessVolumeCavities(int gridpoints,vector<vector<Coordinates> > current, vector<vector<Coordinates> >& original, float spacing);
		void writeVolumeAllFrames(map<float, vector<int> > volume_cavites, float spacing, FILE* file);
		void writeExactVolumeFrame(vector<Coordinates> cavity, Grid& grid, float spacing, float time, FILE* file);
		void writeSolventFrame(int sol, float spacing, float time, FILE* file);
		void writeBottleneckAreaFrame(float area, float time, FILE* file);
		void writeSectionAreaFrame(map<int, pair<vector<float>, float> > section, int nframes, FILE* file);
		virtual ~AtomWriter();
};

#endif /* ATOMWRITER_H_ */
