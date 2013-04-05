/*
 * AtomWriter.cpp
 *
 *  Created on: Dec 1, 2010
 *      Author: tp334
 */

#include "AtomWriter.h"
#include "stdio.h"
#include <string>
#include <sstream>
#include <iostream>
#include "math.h"

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
		if(!cavity.empty() && cavity.size()>1){
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

	//float volume[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

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

vector<int> AtomWriter::preprocessVolumeCavities(int gridpoints, vector<vector<Coordinates> > current, vector<vector<Coordinates> >& original, float spacing){

	vector<int> result (original.size(), 0);
	Coordinates coor1, coor2;
	//vector<vector<Coordinates> > ordered (original.size(), cavity);

	bool same_cavity = false;
	unsigned int c1, c2, i, j;

	if(original.size()>0){
		for(c1=0; c1<original.size();c1++){
			//for each cavity
			for(i=0; i<original[c1].size();i++){
				//for each coordinate check if it is represented in any of the current cavities
				coor1 = original[c1][i];
				c2 = 0;

				while(current.size()>0 && c2<current.size()){
					same_cavity = false;
					for(j=0; j<current[c2].size();j++){
						coor2 = current[c2][j];
						if(fabs(coor1.getX()-coor2.getX())<spacing && fabs(coor1.getY()-coor2.getY())<spacing && fabs(coor1.getZ()-coor2.getZ())<spacing){
							//This current coordinate is in an already existing cavity
							if(!same_cavity){
								same_cavity = true;
								result[c1] = current[c2].size();
								//ordered[c1] = current[c2]);
							}
						}else{
							if(same_cavity){
								original[c1].push_back(coor2);
							}
						}
					}

					if(same_cavity){
						current.erase(current.begin()+c2);
						break;
					}else{
						c2 ++;
					}
				}

				if(same_cavity) break;
			}
		}

		if(current.size()>0){
			for(unsigned int c=0; c<current.size(); c++){
				original.push_back(current[c]);
				result.push_back(current[c].size());
				//ordered.push_back(current[c]);
			}
		}

	}else{
		original = current;
		for(unsigned int c1=0; c1<original.size();c1++){
			result.push_back(original[c1].size());
		}
	}

	//I add the total size
	result.insert(result.begin(), gridpoints);

	return result;
}

void AtomWriter::writeVolumeAllFrames(map<float, vector<int> > volume_cavities, float spacing, FILE* file){

	float voxel_vol = spacing*spacing*spacing;
	vector<int> cavities;
	int size, max = 0;

	map<float, vector<int> >::iterator it;

	for(it=volume_cavities.begin(); it!=volume_cavities.end();++it){
		size = it->second.size();
		if(size>max) max = size;
	}

	fprintf(file, "@ s0 legend \"Total\"\n");
	for(int i=1; i<max; i++){
		fprintf(file, "@ s%d legend \"Cavity %d\"\n",i,i);
	}

	for(it=volume_cavities.begin(); it!=volume_cavities.end();++it){
		cavities = it->second;
		fprintf(file,"%10.3f ", it->first, voxel_vol*cavities[0]);
		for(int i=0; i<max; i++){
			if(i<cavities.size()){
				fprintf(file,"%8.3f ", voxel_vol*cavities[i]);
			}else{
				fprintf(file,"%8.3f ", 0.0);
			}
		}
		fprintf(file,"\n");
	}
	fflush(file);
}


void AtomWriter::writeSolventFrame(int sol, float spacing, float time, FILE* file){
	fprintf(file,"%8.3f, %8i\n",  time, sol);
	fflush(file);
}

void AtomWriter::writeBottleneckAreaFrame(float area, float time, FILE* file){
	fprintf(file,"%8.3f, %8.3f\n",  time, area);
	fflush(file);
}

void AtomWriter::writeSectionAreaFrame(map<int,pair<vector<float>, float> > section, int nframes, FILE* file){
	map<int,pair<vector<float>,float> >::iterator it;
	float mean, sd, n = 0;

	for(it=section.begin(); it!=section.end(); ++it){
		//Calculation of the population standard deviation
		n = it->second.first.size();
		mean = (it->second.second)/(float)n;
		sd = 0;
		for(int i=0; i<n; i++){
			sd = sd + (it->second.first[i]-mean)*(it->second.first[i]-mean);
		}
		sd = sqrt(sd/n);

		fprintf(file,"%d %8.3f %8.3\n",  it->first, mean, sd);
	}

	fflush(file);
}

AtomWriter::~AtomWriter() {
	// TODO Auto-generated destructor stub
}
