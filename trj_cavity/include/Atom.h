/*
 * Atom.h
 *
 *  Created on: Nov 17, 2010
 *      Author: tp334
 */

#ifndef ATOM_H_
#define ATOM_H_

#include <string>
#include "Coordinates.h"
#include <map>

using namespace std;

class Atom {

	static float H;
	static float N;
	static float C;
	static float S;
	static float P;
	static float O;
	static float K;
	static float Na;
	static float Cl;

	protected:
		Coordinates coordinates;
		string name;
		float vdwRadius; //Van der Waals radius in angstrons
		int resId;
		bool ca;

	public:
		Atom();
		Atom(char name,Coordinates coordinates);
		Atom(char name,Coordinates coordinates, float radius);
		static void initialize(char *path);
		//Atom(char name,Coordinates coordinates, int resId);
		//Atom(char name,Coordinates coordinates, string atom, string residue, map<string,map<string, double> > radii);
		void setDefaultVdWRadius(char name);
		//void setForceFieldVdWRadius(char name, string atom, string residue, map<string,map<string, double> > radii);
		float getVdwRadius();
		void setCoordinates(Coordinates coordinates);
		Coordinates getCoordinates();
		void setName(char name);
		string getName();
		void setResId(int resId);
		int getResId();
		virtual ~Atom();
		bool getCa();
		void setCa(bool ca);
};

#endif /* ATOM_H_ */
