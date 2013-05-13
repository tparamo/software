/*
 * Atom.cpp
 *
 *  Created on: Nov 17, 2010
 *      Author: tp334
 */

#include "Atom.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include "stdio.h"
#include <sstream>

using namespace std;

Atom::Atom(){

}

Atom::Atom(char name, Coordinates coordinates) {
	this->setCoordinates(coordinates);
	this->setName(name);
	this->setDefaultVdWRadius(name);
}

Atom::Atom(char name, Coordinates coordinates, float radius) {
	this->setCoordinates(coordinates);
	this->setName(name);
	this->setDefaultVdWRadius(name);
	this->vdwRadius=radius;
}

void Atom::setDefaultVdWRadius(char name){
	if(name == 'H') vdwRadius = 0.40; //wiki raius -> vdwRadius = 1.20;
	else if(name == 'N') vdwRadius = 1.55; // gromacs default radius-> vdwRadius = 1.10;
	else if(name == 'C') vdwRadius = 1.70; // gromacs default radius-> vdwRadius = 1.50;
	else if(name == 'O') vdwRadius = 1.52; // gromacs default radius-> vdwRadius = 1.05;
	else if(name == 'S') vdwRadius = 1.80; // gromacs default radius->vdwRadius = 1.60;
	else if(name == 'P') vdwRadius = 1.80; // gromacs default radius->vdwRadius = 1.20;
	else if(name == 'K') vdwRadius = 2.75; // K ions
	else if(name == '1') vdwRadius = 2.27; // NA ions
	else if(name == '2') vdwRadius = 1.75; // CL ions
	else vdwRadius = 1.20;
}

float Atom::getVdwRadius(){
	return vdwRadius;
}

void Atom::setName(char name){
	this->name = name;
}

string Atom::getName(){
	return name;
}

void Atom::setCoordinates(Coordinates coordinates){
	this->coordinates.setX(coordinates.getX());
	this->coordinates.setY(coordinates.getY());
	this->coordinates.setZ(coordinates.getZ());
}

Coordinates Atom::getCoordinates(){
	return coordinates;
}

void Atom::setResId(int res){
	this->resId = res;
}

int Atom::getResId(){
	return resId;
}

bool Atom::getCa()
{
    return ca;
}

void Atom::setCa(bool ca)
{
    this->ca = ca;
}

Atom::~Atom() {
	// TODO Auto-generated destructor stub
}
