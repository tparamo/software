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
#include "AtomRadii.h"

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
	switch(name){
		case 'H' : vdwRadius = H; break;
		case 'N' : vdwRadius = N; break;
		case 'C' : vdwRadius = C; break;
		case 'O' : vdwRadius = O; break;
		case 'S' : vdwRadius = S; break;
		case 'P' : vdwRadius = P; break;
		case 'K' : vdwRadius = K; break;
		case '1' : vdwRadius = Na; break;
		case '2' : vdwRadius = Cl; break;
		default: vdwRadius = 1.20;
	}
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
