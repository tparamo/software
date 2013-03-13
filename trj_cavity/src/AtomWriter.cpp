/*
 * AtomWriter.cpp
 *
 *  Created on: Dec 1, 2010
 *      Author: tp334
 */

#include "AtomWriter.h"
#include "math.h"
#include "stdio.h"
#include <string>
#include <sstream>
#include <iostream>

AtomWriter::AtomWriter() {
	// TODO Auto-generated constructor stub
}

void AtomWriter::writePDB(vector<Coordinates> cavities, FILE* file){
	for(unsigned int i=0;i<cavities.size();i++){
		fprintf(file, "ATOM %6i PNT CAV 9999     %8.3f%8.3f%8.3f  1.00 0.00\n",
				i+1,cavities[i].getX(),cavities[i].getY(),cavities[i].getZ());
	}
	fprintf(file, "END");
	fclose(file);
}

void AtomWriter::writePDB(vector<vector<Coordinates> > cavities, FILE* file){
	int cont = 1;
	for(unsigned int c=0;c<cavities.size();c++){
		vector<Coordinates> cavity = cavities[c];
		if(cavity.size()>1){
			char chain[3];
			//sprintf(chain,"%3d", c+1);
			for(unsigned int i=0;i<cavity.size();i++){
				fprintf(file, "ATOM%7i PNT  CAV A%4i    %8.3f%8.3f%8.3f  1.00 0.00\n", cont, c+1, cavity[i].getX(),cavity[i].getY(),cavity[i].getZ());
				cont ++;
			}
		}
	}
	fprintf(file, "END");
	fflush(file);
}

void AtomWriter::writeStatistics(Grid statistics, float gridspacing, int nframes, FILE* file){
	stringstream ss;

	float volume[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	int atom = 1;

	for(int i=0; i<statistics.getWidth();i++){
		for(int j=0; j<statistics.getHeight();j++){
			for(int k=0; k<statistics.getDepth();k++){
				if(statistics.getGrid()[i][j][k]!=0){
					int* position = new int[3];
					position[0] = i;
					position[1] = j;
					position[2] = k;
					Coordinates c = statistics.calculateCoordinateInGrid(position);

					float pc = (float)statistics.getGrid()[i][j][k]/(float)(nframes);

					/*if(pc>=0.1) volume[0] = volume[0] + (gridspacing*gridspacing*gridspacing);
					if(pc>=0.2) volume[1] = volume[1] + (gridspacing*gridspacing*gridspacing);
					if(pc>=0.3) volume[2] = volume[2] + (gridspacing*gridspacing*gridspacing);
					if(pc>=0.4) volume[3] = volume[3] + (gridspacing*gridspacing*gridspacing);
					if(pc>=0.5) volume[4] = volume[4] + (gridspacing*gridspacing*gridspacing);
					if(pc>=0.6) volume[5] = volume[5] + (gridspacing*gridspacing*gridspacing);
					if(pc>=0.7) volume[6] = volume[6] + (gridspacing*gridspacing*gridspacing);
					if(pc>=0.8) volume[7] = volume[7] + (gridspacing*gridspacing*gridspacing);
					if(pc>=0.9) volume[8] = volume[8] + (gridspacing*gridspacing*gridspacing);*/


					fprintf(file, "ATOM %6i PNT CAV  9999     %8.3f%8.3f%8.3f  1.00%6.2f\n",
						atom,c.getX(),c.getY(),c.getZ(), pc);

					atom ++;
				}
			}
		}
	}

	//writeVolumeProbability(volume);

	fprintf(file, "END");
	fflush(file);
	fclose(file);
}

void AtomWriter::writeVolumeProbability(float volume[]){
	/*FILE *file = fopen (fpath,"w");

	for(int i=0; i<9; i++){
		fprintf(file,"%8.3f,%8.3f\n", volume[i], (i+1)*0.1);
	}

	fclose(file);*/


	for(int i=0; i<9; i++){
		stringstream ss;
		ss<<"Volume with probability >="<<(i+1)*0.1<<" :"<<volume[i]<<" A3";
	}
}

void AtomWriter::writeVolumeFrame(int gridpoints, float spacing, float time, FILE* file){
	fprintf(file,"%8.3f, %8.3f\n",  time, spacing*spacing*spacing*gridpoints);
	fflush(file);
}

void AtomWriter::writeExactVolumeFrame(vector<Coordinates> cavity, Grid& grid, float spacing, float time, FILE* file){

	float acc = 0;

	for(unsigned int i=0; i<cavity.size();i++){
		//acc = acc + grid.refineVolume(cavity[i].getX(), cavity[i].getY(), cavity[i].getZ());
	}

	fprintf(file,"%8.3f, %8.3f\n",  time, acc);
	fflush(file);
}

void AtomWriter::writeVolumeFrame(vector<vector<Coordinates> > current, vector<vector<Coordinates> > original, float spacing, float time, FILE* file){

	float result[original.size()];
	float voxel = spacing*spacing*spacing;

	for(unsigned o=0; o<original.size(); o++){
		result[o] = 0.0;
	}

	bool exists;
	for(unsigned o=0; o<original.size(); o++){
		if(time==0.0){
			result[o] = original[o].size()*voxel;
		}else{
			for(unsigned c=0; c<current.size(); c++){
				exists = false;
				for(unsigned i=0; i<current[c].size(); i++){
					for(unsigned j=0; j<original[o ].size(); j++){
						if(fabs(original[o][i].getX()-current[c][j].getX())<spacing and fabs(original[o][i].getY()-current[c][j].getY())<spacing and fabs(original[o][i].getZ()-current[c][j].getZ())<spacing){
							result[o] = current[c].size()*voxel;
							exists = true;
							break;
						}
					}
					if(exists) break;
				}
			}
		}
	}

	fprintf(file,"%8.3f", time);
	for(unsigned o=0; o<original.size(); o++){
		fprintf(file," %8.3f",result[o]);
	}
	fprintf(file,"\n");

	fflush(file);
}


void AtomWriter::writeWatersFrame(int waters, float spacing, float time, FILE* file){
	fprintf(file,"%8.3f, %8i\n",  time, waters);
	fflush(file);
}

void AtomWriter::writeSectionAreaFrame(float area, float time, FILE* file){
	fprintf(file,"%8.3f, %8.3f\n",  time, area);
	fflush(file);
}

AtomWriter::~AtomWriter() {
	// TODO Auto-generated destructor stub
}
