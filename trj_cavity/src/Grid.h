/*
 * Grid.h
 *
 *  Created on: Nov 19, 2010
 *      Author: tp334
 */

#ifndef GRID_H_
#define GRID_H_

#include <vector>
#include <set>
#include "Atom.h"
#include "Coordinates.h"

using namespace std;

class Grid {
	protected:
		int ***grid;
		int width;
		int height;
		int depth;
		float originX;
		float originY;
		float originZ;
		float spacing;
		int ca_limits[6];
		float tStartStatistics;
		int dimensionsSearch;
		int cutoff;
		vector<Atom> atoms;

	public:
		Grid();
		Grid(float spacingGrid,float x_max, float x_min, float y_max, float y_min, float z_max, float z_min, float cutoff);
		//do i need this anymore?
		Grid(float spacingGrid, float x_min, float y_min, float z_min, int width, int height, int depth);
		void initializeGrid();
		float calculateDistance(float a, float b);
		float calculateNumberCells(float a, float b);
		int calculatePositionInGrid(float coordinate, float origin);
		int* getDistanceCells(float radius, int center, int limit);
		vector<vector<Coordinates> > getCavitySeedPoint(Coordinates coordinates, Grid statistics, int &max_cavity_index, int min, bool rectify);
		vector<vector<Coordinates> > getCavitiesNoSeedPoint(Grid statistics, int &max_cavity_index, int &max_cavity_position, int min, bool rectify);
		vector<vector<Coordinates> > getCavitiesSeedPoint(Coordinates seed, Grid statistics, int &max_cavity_index, int &max_cavity_position, int min, bool rectify);
		vector<vector<Coordinates> > getCavities(Grid statistics, int &max_cavity_index, int &max_cavity_position, bool rectify, int x_min, int x_max, int y_min, int y_max, int z_min, int z_max, int min);
		vector<Coordinates> getCavityEvolution(float cutoff,vector<Coordinates> last_cavity, int ***statistics, bool rectify);
		vector<Coordinates> getCavity(int x, int y , int z, bool rectify, Grid statistics, int index);
		void getNeighbourCandidate(int x, int y, int z, vector<int*>& neighbours, bool rectify, Grid statistics, int index);
		vector<int*> getNeighbours(int* seed, bool rectify, Grid statistics, int index);
		Coordinates calculateCoordinateInGrid(int* position);
		bool isInsideCavity(int x, int y, int z);
		int surroundedBy(int x, int y, int z, int element, int positions);
		bool isInsideCavityNew(int x, int y, int z, int pos, int index);
		int surroundedByX(int x, int y, int z, int element, int pos);
		int surroundedByY(int x, int y, int z, int element, int pos);
		int surroundedByZ(int x, int y, int z, int element, int pos);
		int surroundedByXPositive(int x, int y, int z, int element, int pos);
		int surroundedByXNegative(int x, int y, int z, int element, int pos);
		int surroundedByYPositive(int x, int y, int z, int element, int pos);
		int surroundedByYNegative(int x, int y, int z, int element, int pos);
		int surroundedByZPositive(int x, int y, int z, int element, int pos);
		int surroundedByZNegative(int x, int y, int z, int element, int pos);
		vector<int*> getCandidatesToSeedPoint(int x, int y , int z, float cutoff);
		int*** getGrid();
		void deleteGrid();
		int getWidth();
		int getHeight();
		int getDepth();
		float getOriginX();
		float getOriginY();
		float getOriginZ();
		float getSpacing();
		void setCAplhaLimits(float x_max, float x_min, float y_max, float y_min, float z_max, float z_min);
		void setTStartStatisitics(float tstat);
		float getTStartStatisitics();
		void setDimensionsSearch(int dim);
		int getDimensionsSearch();
		void setAxisDependentInfo(int axis, int &ilimit, int &jlimit, int &klimit, float &origin, int &a, int &b);
		map<int,float> calculateBottleneckArea(int index,int axis, vector<Coordinates>& tunnel, float &bottleneck, float axis_value ,bool sector, bool radius);
		map<int, float> getSectorTunnel(int index, int axis, int a, int b, vector<Coordinates>& tunnel, float &bottleneck, bool radius, float axis_value);
		float getAreaAxis(int* seed, int index, int axis, int a, int b);
		virtual ~Grid();
		vector<Atom> getAtoms();
		void setAtoms(vector<Atom> atoms);
		//float getAreaAxis(vector<Coordinates> cavity, int axis, float value);
		//float refineVolume(float x, float y, float z);
		float refineArea(float xf, float yf, float zf, int axis);
};

#endif /* GRID_H_ */
