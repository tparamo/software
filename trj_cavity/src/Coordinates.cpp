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

void Coordinates::setGridCoordinate(int *coor){
	this->gridCoor = coor;
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

int* Coordinates::getGridCoordinate(){
	return this->gridCoor;
}

float Coordinates::getCoordinate(int axis){
	switch (axis){
		case 0:
			return this->x; break;
		case 1:
			return this->y; break;
		case 2:
			return this->z; break;
		default:
			return 0;
	}
}

Coordinates::~Coordinates() {
	//delete this->gridCoor;
}
