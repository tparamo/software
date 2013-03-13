//============================================================================
// Name        : CavityCC.cpp
// Author      : T Paramo
// Version     :
// Copyright   : 
// Description :
//============================================================================

#include "stdio.h"
#include <iostream>
#include <vector>
#include "src/AtomReader.h"
#include "src/AtomWriter.h"
#include <time.h>
#include <sstream>
#include "string.h"


using namespace std;

int main(int argc, char *argv[]) {

	/*zchar* topology;
	char* trajectory;
	bool trajectory_m = false;
	string fpath = "";

	bool gromos = false;;
	float spacing;
	int mode = 0; //mode==0 stands for biggest cavity, mode==1 stands for all cavities
	int tstat = 0;
	float  minvol=1;

	float seed_x = 0.0;
	float seed_y = 0.0;
	float seed_z = 0.0;

	for(int i = 1; i < argc; i++){
		if(argv[i][0]=='-'){
			string sw = argv[i];
			if(sw=="-s") {
				topology = argv[i+1];
				string inputfile = argv[i+1];
				fpath = inputfile.substr(0, inputfile.find_last_of('/'));
				log.setPath(fpath);
				log.info("Start. Input file: "+inputfile);
			}
			else if(sw=="-f"){
				trajectory = argv[i+1];
				trajectory_m = true;
			}
			else if(sw=="-mode"){sscanf(argv[i+1],"%i",&mode);}
			else if(sw=="-sp") {sscanf(argv[i+1],"%f",&spacing);}
			else if(sw=="-gromos"){gromos = true;}
			else if(sw=="-mv"){sscanf(argv[i+1],"%f", &minvol);}
			else if(sw=="-tstat"){sscanf(argv[i+1],"%i", &tstat);}
			else if(sw=="-coords"){;
				sscanf(argv[i+1],"%f", &seed_x);
				sscanf(argv[i+2],"%f", &seed_x);
				sscanf(argv[i+3],"%f", &seed_x);
			}
		}
	}

	FILE *file = fopen (topology,"r");

	if(file == NULL){
		log.error("No file!");
	}else{
		Grid statistics;
		AtomWriter writer;
		AtomReader reader(fpath);

		if(!trajectory_m){
			Grid grid;
			Grid statistics;
			vector<vector<Coordinates> > v;

			grid = reader.readPDB(file, spacing, gromos);

			if (mode==0){
				v = grid.getCavityNoSeedPoint(statistics, false);
			}else{
				v = grid.getCavitiesNoSeedPoint(statistics, false);
			}

			if(v.size()>0){
				writer.writePDB(v,reader.getCavityFile());
			}

			grid.deleteGrid();

		}else{
			string cav_traj = reader.getFPath()+"/"+"trajectory.xtc";
			FILE* cav = reader.getCavityFile();
			map<string,map<string,double> > forceFieldRadii;
			if(gromos){
			  ForceField ff("charmm27");
			  forceFieldRadii = ff.getForceFieldRadii();
			}
			reader.readTrajectory(trajectory, topology, spacing, cav, cav_traj.c_str() , reader.getVolumeFile(), reader.getStatisticsFile() , mode, tstat, gromos,forceFieldRadii, reader.getWatersFile(), reader.getWaterStatisticsFile());
		}
	}

	log.info("End");*/

	AtomReader reader;
	reader.start(argc,argv);

	return 0;
}
