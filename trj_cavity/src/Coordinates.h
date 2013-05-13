/*
 * Coordinates.h
 *
 *  Created on: Nov 30, 2010
 *      Author: tp334
 */

#ifndef Coordinates_H_
#define Coordinates_H_

#include "math.h"

using namespace std;

class Coordinates {
	protected:
		float x;
		float y;
		float z;
		int gridCoor[3];

	public:
		Coordinates();
		Coordinates(float x, float y, float z);
		void setX(float x);
		void setY(float y);
		void setZ(float z);
		void setGridCoordinate(int *coor);
		float getX();
		float getY();
		float getZ();
		int* getGridCoordinate();
		float getCoordinate(int axis);
		~Coordinates();
};

class compare_float {
	private:
      int spacing;
   public:
      compare_float(float s_) : spacing(s_) {};
      bool operator()(const float x,const float y) {
    	  return (fabs(x-y)<spacing);
      }
};

#endif /* Coordinates_H_ */
