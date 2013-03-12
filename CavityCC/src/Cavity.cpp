/*
 * Cavity.cpp
 *
 *  Created on: Jan 15, 2013
 *      Author: tp334
 */

#include "Cavity.h"

Cavity::Cavity() {
	// TODO Auto-generated constructor stub

}

Cavity::Cavity(int id) {
	this->id = id;
}

Cavity::~Cavity() {
	// TODO Auto-generated destructor stub
}

double Cavity::getBottleneckRadius()
{
    return bottleneckRadius;
}

vector<Coordinates> Cavity::getCoordinates()
{
    return coordinates;
}

int Cavity::getId()
{
    return id;
}

bool Cavity::getTunnel()
{
    return tunnel;
}

void Cavity::setBottleneckRadius(double bottleneckRadius)
{
    this->bottleneckRadius = bottleneckRadius;
}

void Cavity::setCoordinates(vector<Coordinates> coordinates)
{
    this->coordinates = coordinates;
}

void Cavity::setId(int id)
{
    this->id = id;
}

void Cavity::setTunnel(bool tunnel)
{
    this->tunnel = tunnel;
}


