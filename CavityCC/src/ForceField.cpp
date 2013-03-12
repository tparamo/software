/*
 * ForceField.cpp
 *
 *  Created on: Mar 22, 2012
 *      Author: tp334
 */

#include "ForceField.h"
#include "Log.h"
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include "stdio.h"

using namespace std;

ForceField::ForceField(){

}

ForceField::ForceField(string forcefield) {
	this->setForceFieldPath(forcefield);
	this->setForceFieldRadii();
}

void ForceField::setForceFieldPath(string forcefield){
	string lib = getenv("LD_LIBRARY_PATH");
	string gromacs = lib.substr(0, lib.find_last_of('/'));
	string ff = gromacs + "/" + "share/top/" + forcefield + ".ff";
	Log log;
	log.info(ff);
	this->forceFieldPath =  ff;
}

string ForceField::getForceFieldPath(){
       return this->forceFieldPath;
}

void ForceField::setForceFieldRadii(){
	Log log;
	stringstream ss;

	string aa = this->forceFieldPath + "/aminoacids.rtp";

	FILE* aa_file = fopen(aa.c_str(),"r");

    map<string,double> radii = this->calculateRadii();
	map<string, double> atomsmap;

	char line[100];

	if(aa_file!=NULL){
		bool atoms = false;
		string residue = "";
		while (fgets(line, sizeof(line), aa_file)) {
			 char str[2];
			 strncpy(str,line,2);
			 string aux = "[ ";
			 if(strstr(str,aux.c_str())!=NULL && strstr(line,"[ bondedtypes ]")==NULL){
				 string l = line;
				 residue = l.substr(2,3);
				 continue;
			 }else{
				 if(atoms==false){
					 char str2[10];
					 strncpy(str2,line,10);
					 string aux2 = " [ atoms ]";
					 if(strstr(str2,aux2.c_str())!=NULL){
						 atoms = true;
						 continue;
					 }
				 }else{
					 char str3[10];
					 strncpy(str3,line,10);
					 string aux3 = " [ bonds ]";
					 if(strstr(str3,aux3.c_str())!=NULL or str3[0]=='\n'){
						 forceFieldRadii[residue] = atomsmap;
						 atomsmap.clear();
						 atoms = false;
						 continue;
					 }else{
						 char name[5], type[5];
						 float charge;
						 int chargegroup;
						 sscanf(line,"%s %s %f %i", name, type, &charge, &chargegroup);
						 string n = name;
						 string t = type;
						 atomsmap[n]=radii[t];

						 ss<<n<<" "<<residue<<" "<<t<<" "<<atomsmap[n]<<" "<<endl;
					 }
				 }

			 }
		}

		log.info("Force field van der waals radius:");
		log.info(ss.str());

		fclose(aa_file);
	}else{
		throw 1;
	}
}

map<string, map<string, double> > ForceField::getForceFieldRadii(){
        return this->forceFieldRadii;
}


map<string, double> ForceField::calculateRadii(){
	string ffnonbonded = this->forceFieldPath + "/ffnonbonded.itp";

	map<string, double> radii;

	FILE* nb_file = fopen(ffnonbonded.c_str(),"r");

	Log log;
	stringstream ss;
	ss<<"Calculating van der waals radii according to forcefield "<<this->forceFieldPath<<endl;

	if(nb_file!=NULL){
		char line[100];
		bool atoms = false;
		bool sigma = false;

		while (fgets(line, sizeof(line), nb_file)) {
			if(strstr(line,"[ atomtypes ]")){
				atoms = true;
				continue;
			}else{
				if(atoms==true){
					if(strstr(line,"[ nonbond_params ]")){
						atoms = false;
						break;
					}else{
						if(strstr(line,"sigma")){
							sigma = true;
						}
						char str[1];
						strncpy(str,line,1);
						string aux = ";";
						if(strstr(str,aux.c_str())!=NULL){
							continue;
						}else{
							char name[5], ptype[1];
							int atnum;
							float mass, charge, c6, c12;
							sscanf(line,"%s %i %f %f %s %f %f", name,&atnum,&mass,&charge,ptype,&c6,&c12);

							double radio = 0.0;
							double exponent =1.0/6.0;

							if(sigma){
								radio = (pow(2.0, exponent)*c6)/2;
							}else{
								double base = 2.0*c12/c6;
								radio = pow(base,exponent)/2.0;
							}

							radii[name] = radio*10.0; //convert to angstroms

							ss<<name<<" "<<radii[name]<<endl;
						}
					}
				}
			}
		}
		fclose(nb_file);
	}else{
		throw 2;
	}

	//log.info(ss.str());

	return radii;
}

ForceField::~ForceField() {
	// TODO Auto-generated destructor stub
}
