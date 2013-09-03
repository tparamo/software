/*
 * ForceField.h
 *
 *  Created on: Mar 22, 2012
 *      Author: tp334
 */

#include <string>
#include <map>

#ifndef FORCEFIELD_H_
#define FORCEFIELD_H_

using namespace std;

class ForceField {
	protected:
           map<string,map<string,double> > forceFieldRadii;
	       string forceFieldPath;

	public:
		ForceField();
        ForceField(string forcefield);
		map<string,map<string,double> > getForceFieldRadii();
		void setForceFieldRadii();
		map<string, double> calculateRadii();
		void setForceFieldPath(string forcefield);
        string getForceFieldPath();
		virtual ~ForceField();

};

#endif /* FORCEFIELD_H_ */
