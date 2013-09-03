/*
 * Atom.cpp
 *
 *  Created on: Nov 17, 2010
 *      Author: tp334
 */


#include "../include/Atom.h"

#include "stdio.h"
#include <cstring>
#include "cstdlib"


using namespace std;

float Atom::H;
float Atom::N;
float Atom::C;
float Atom::S;
float Atom::P;
float Atom::O;
float Atom::K;
float Atom::Na;
float Atom::Cl;

Atom::Atom(){
}

void Atom::initialize(char *path){

	char *pch = strrchr(path,'/');
	strncpy (pch,"/AtomRadii.txt",14);

	FILE* file = fopen(path,"r");
	if(!file){
		printf("File AtomRadii.txt not found!\n");
	    exit(1);
	}else{
		char atom[2];
		float radius;
		char line[10];
		while (fgets(line, sizeof(line), file)) {
			sscanf(line,"%s %f", atom, &radius);
			if(strstr(atom,"H")) {H=radius;}
			if(strstr(atom,"C")) {C=radius;}
			if(strstr(atom,"N")) {N=radius;}
			if(strstr(atom,"S")) {S=radius;}
			if(strstr(atom,"P")) {P=radius;}
			if(strstr(atom,"O")) {O=radius;}
			if(strstr(atom,"K")) {K=radius;}
			if(strstr(atom,"Na")) {Na=radius;}
			if(strstr(atom,"Cl")) {Cl=radius;}
		}
		fclose(file);
	}
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
