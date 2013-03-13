/*
 * Coordinates.h
 *
 *  Created on: Nov 30, 2010
 *      Author: tp334
 */

#ifndef Coordinates_H_
#define Coordinates_H_

using namespace std;

class Coordinates {
	protected:
		float x;
		float y;
		float z;

	public:
		Coordinates();
		Coordinates(float x, float y, float z);
		void setX(float x);
		void setY(float y);
		void setZ(float z);
		float getX();
		float getY();
		float getZ();
		~Coordinates();
};

#endif /* Coordinates_H_ */
