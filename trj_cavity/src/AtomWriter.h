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
		void writeVolumeFrame(vector<vector<Coordinates> > current, vector<vector<Coordinates> > original, float spacing, float time, FILE* file);
		void writeExactVolumeFrame(vector<Coordinates> cavity, Grid& grid, float spacing, float time, FILE* file);
		void writeWatersFrame(int waters, float spacing, float time, FILE* file);
		void writeSectionAreaFrame(float area, float time, FILE* file);
		virtual ~AtomWriter();
};

#endif /* ATOMWRITER_H_ */
