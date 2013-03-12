/*
 * Cavity.h
 *
 *  Created on: Jan 15, 2013
 *      Author: tp334
 */

#ifndef CAVITY_H_
#define CAVITY_H_

#include <vector>
#include "Coordinates.h"

using namespace std;

class Cavity {
	protected:
		int id;
		vector<Coordinates> coordinates;
		double bottleneckRadius;
		bool tunnel;

	public:
		Cavity();
		Cavity(int id);
		double getBottleneckRadius();
		vector<Coordinates> getCoordinates();
		int getId();
		bool getTunnel();
		void setBottleneckRadius(double bottleneckRadius);
		void setCoordinates(vector<Coordinates> coordinates);
		void setId(int id);
		void setTunnel(bool tunnel);
		~Cavity();
};

#endif /* CAVITY_H_ */
