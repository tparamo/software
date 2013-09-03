/*
 * ForceField.cpp
 *
 *  Created on: Mar 22, 2012
 *      Author: tp334
 */

#include "../include/ForceField.h"

#include <iostream>
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
	this->forceFieldPath =  forcefield;
}

string ForceField::getForceFieldPath(){
       return this->forceFieldPath;
}

void ForceField::setForceFieldRadii(){

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
						 if(forceFieldRadii.find(residue)==forceFieldRadii.end()){
							 forceFieldRadii[residue] = atomsmap;
						 }
						 atomsmap.clear();
						 atoms = false;
						 continue;
					 }else{
						 char name[10], type[10];
						 float charge;
						 int chargegroup;
						 sscanf(line,"%s %s %f %i", name, type, &charge, &chargegroup);
						 string n = name;
						 string t = type;
						 if(atomsmap.find(n)==atomsmap.end()){
							 atomsmap[n]=radii[t];
						 }

						 //cout<<n<<" "<<residue<<" "<<t<<" "<<atomsmap[n]<<" "<<endl;
					 }
				 }

			 }
		}

		fclose(aa_file);
	}else{
		printf("File %s not found", aa.c_str());
		exit(1);
	}
}

map<string, map<string, double> > ForceField::getForceFieldRadii(){
        return this->forceFieldRadii;
}


map<string, double> ForceField::calculateRadii(){
	string ffnonbonded = this->forceFieldPath + "/ffnonbonded.itp";

	map<string, double> radii;

	FILE* nb_file = fopen(ffnonbonded.c_str(),"r");

	cout<<"Calculating van der waals radii according to forcefield "<<this->forceFieldPath<<endl;

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
							char name[10], ptype[5];
							int atnum;
							float mass, charge, c6, c12;
							sscanf(line,"%s %i %f %f %s %f %f", name,&atnum,&mass,&charge,ptype,&c6,&c12);

							if(radii.find(name)==radii.end()){
								if(sigma){
									//radii[name] = (0.5*c6)*10;
									radii[name] = (0.5*(pow(2.0,(1.0/6.0))*c6))*10;
								}else{
									float sig6 = c12/c6;
									radii[name]= (0.5*pow(sig6,1.0/6.0))*10;
								}
							}
						}
					}
				}
			}
		}
		fclose(nb_file);
	}else{
		printf("File %s not found", ffnonbonded.c_str());
		exit(1);
	}

	return radii;
}

ForceField::~ForceField() {
	// TODO Auto-generated destructor stub
}
