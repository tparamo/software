/*
 * Coordinates.cpp
 *
 *  Created on: Nov 30, 2010
 *      Author: tp334
 */

#include "Coordinates.h"

Coordinates::Coordinates() {
	// TODO Auto-generated constructor stub

}

Coordinates::Coordinates(float x, float y, float z) {
	this->x = x;
	this->y = y;
	this->z = z;
}

void Coordinates::setX(float x) {
	this->x = x;
}

void Coordinates::setY(float y) {
	this->y = y;
}

void Coordinates::setZ(float z) {
	this->z = z;
}


float Coordinates::getX() {
	return x;
}

float Coordinates::getY() {
	return y;
}

float Coordinates::getZ() {
	return z;
}

Coordinates::~Coordinates() {
	// TODO Auto-generated destructor stub
}
