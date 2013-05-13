/*
 * AromReader.cpp
 *
 *  Created on: Nov 17, 2010
 *      Author: tp334
 */

#include "AtomReader.h"
#include "AtomWriter.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include <string>
#include "ctype.h"
#include <stdlib.h>

//Gromacs libraries
#include "xtcio.h"
#include <sstream>
#include "trajana.h"

#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "copyrite.h"
#include "tpxio.h"
#include "names.h"
#include "gmx_random.h"
#include "gmx_ana.h"
#include "filenm.h"
#include "macros.h"
#include "pbc.h"
#include "vec.h"
#include "xvgr.h"

using namespace std;

AtomReader::AtomReader() {
	// TODO Auto-generated constructor stub

}

AtomReader::~AtomReader() {
}

struct t_file_pointers
{
	//string fpath; /* where the files might be saved*/
	FILE *fvolume; /* volume file*/
	FILE *fcavity; /* cavity file*/
	FILE *ftunnel; /* tunnel file */
	FILE *fsolcount; /* water count file*/
	FILE *fbottleneck; /* tunnel bottleneck area*/
	FILE *log; /* log file */
	t_trxstatus *ftrajectory; /* cavity trajectory file*/
	t_trxstatus *ftunnel_trajectory; /* tunnel trajectory file*/
};

struct t_cavity_params
{
	atom_id *index;   /* the index for the atom numbers */
	int isize;    /* the size of each group */
	float *radii; /* ff radii */
	atom_id *sol_index;   /* the index for the solvent  */
	int sol_isize;    /* the size of each group */
	float *sol_radii; /* ff radii for solvent */
	Grid statistics; /* statistics*/
	Grid water_statistics; /* water statistics */
	int min_size; /* cavity min size for trajectory*/
	int tunnel_size; /* tunnel min size*/
	int frame_stat;
	int frame;
	vector<vector<Coordinates> > v; /* Set of reference cavities*/
	map<float, vector<int> > volume_cavities; /* Set of cavities for each frame */
	map<int,pair<vector<float>,float> > acc_sector; /* tunnel section statistics */
};


struct t_cavity_options
{
	int mode; /* Mode of execution*/
	float min; /* Minimun size of cavities */
	int dim; /* Number of dimensions to cover by the protein */
	int skip; /* Frames to skip when processing trajectory */
	float spacing; /* Grid spacing */
	float t_stat; /* Time to strat calculating statistics */
	bool rectify; /* Rectify output?*/
	bool ff_radius; /* use force field radius*/
	bool sol; /* water monitoring */
	bool ca; /* Restrict to ca*/
	bool tunnel; /* water monitoring */
	rvec seed; /* seed point */
	int axis; /* Tunnel orientation */
	bool use_seed;
	float area; /* Axis value to inspect area*/
	float cutoff; /* Axis value to inspect area*/

	t_file_pointers files; /* file pointers */

    t_cavity_params params; /* params needed for cavity calculation */
};

struct gmx_ana_traj_t
{
    unsigned long             flags;
    int                       nrefgrps;
    int                       nanagrps;
    int                       frflags;
    gmx_bool                  bRmPBC;
    gmx_bool                  bPBC;
    char                     *trjfile;
    char                     *topfile;
    char                     *topfile_notnull;
    char                     *ndxfile;
    char                     *selfile;
    char                     *selection;
    t_topology               *top;
    gmx_bool                  bTop;
    rvec                     *xtop;
    matrix                    boxtop;
    int                       ePBC;
    t_trxframe               *fr;
    t_trxstatus              *status;
    int                       nframes;
    int                       ngrps;
    gmx_ana_selection_t     **sel;
    char                    **grpnames;
    gmx_ana_poscalc_coll_t   *pcc;
    gmx_ana_selcollection_t  *sc;
    output_env_t              oenv;
};

AtomReader::AtomReader(string fpath) {

	/*TODO: is it neccesary? */
	this->fpath = fpath;
}

void AtomReader::calculateLimitsGrid(float *limits, Atom atom){

	 if(atom.getCoordinates().getX() + atom.getVdwRadius()>limits[0]){
		 limits[0] = atom.getCoordinates().getX() + atom.getVdwRadius();
	 }else{
		 if(atom.getCoordinates().getX() - atom.getVdwRadius()<limits[3]) limits[3] = atom.getCoordinates().getX() - atom.getVdwRadius();
	 }

	 if(atom.getCoordinates().getY() + atom.getVdwRadius()>limits[1]){
		 limits[1] = atom.getCoordinates().getY() + atom.getVdwRadius();
	 }else{
		 if(atom.getCoordinates().getY() - atom.getVdwRadius()<limits[4]) limits[4] = atom.getCoordinates().getY() - atom.getVdwRadius();
	 }

	 if(atom.getCoordinates().getZ() + atom.getVdwRadius()>limits[2]){
		 limits[2] = atom.getCoordinates().getZ() + atom.getVdwRadius();
	 }else{
		 if(atom.getCoordinates().getZ() - atom.getVdwRadius()<limits[5]) limits[5] = atom.getCoordinates().getZ() - atom.getVdwRadius();
	 }
}

Grid AtomReader::initializeGrid(Grid grid, vector<Atom> atoms){

	float spacing = grid.getSpacing();
	int cell_center_x, cell_center_y,cell_center_z;
	int *x_cells,*y_cells,*z_cells;
	Atom atom;

	int x, y, z;

	for(unsigned int i=0;i<atoms.size();i++){

		atom = atoms[i];

		cell_center_x = grid.calculatePositionInGrid(atom.getCoordinates().getX(), grid.getOriginX());
		cell_center_y = grid.calculatePositionInGrid(atom.getCoordinates().getY(), grid.getOriginY());
		cell_center_z = grid.calculatePositionInGrid(atom.getCoordinates().getZ(), grid.getOriginZ());

		//Apply the atom radius

		if((atom.getVdwRadius())*2>spacing){
			x_cells = grid.getDistanceCells(atom.getVdwRadius(), cell_center_x, grid.getWidth());
			y_cells = grid.getDistanceCells(atom.getVdwRadius(), cell_center_y, grid.getHeight());
			z_cells = grid.getDistanceCells(atom.getVdwRadius(), cell_center_z, grid.getDepth());
		}else{
			x_cells = new int[2]; x_cells[0] = cell_center_x; x_cells[1] = cell_center_x;
			y_cells = new int[2]; y_cells[0] = cell_center_y; y_cells[1] = cell_center_y;
			z_cells = new int[2]; z_cells[0] = cell_center_z; z_cells[1] = cell_center_z;
		}


		for(x=x_cells[0]; x<=x_cells[1];x++){
			for(y=y_cells[0]; y<=y_cells[1];y++){
				for(z=z_cells[0]; z<=z_cells[1];z++){
					grid.getGrid()[x][y][z] = i+1;
				}
			}
		}

		delete[] x_cells;
		delete[] y_cells;
		delete[] z_cells;
	}

	return grid;
}

/*int* AtomReader::distanceCells(float spacing, float radius, int center, int limit){
	int* result = new int[2];

	float p_distance = radius/spacing;
	int p_positions = (int)p_distance;
	if(center + p_positions <limit-1){
		if((p_distance - p_positions)>0.5){
			p_positions ++;
		}
	}else{
		while(center + p_positions >limit-1 and p_positions>0){
			p_positions --;
		}
	}

	float n_distance = radius/spacing;
	int n_positions = (int)n_distance;
	if(center - n_positions > 0){
		if((n_distance - n_positions)>0.5){
			n_positions ++;
		}
	}else{
		while(center - n_positions < 0 and n_positions>0){
			n_positions --;
		}
	}

	result[0] = center - n_positions;
	result[1] = center + p_positions;

	return result;
}*/

//GROMACS area

void writeTrajectory(t_trxstatus *trjfile, vector<vector<Coordinates> > v, t_trxframe *fr, int size);

void writeTrajectory(t_trxstatus *trjfile, vector<vector<Coordinates> > v, t_trxframe *fr, int size){
	t_trxframe frout;

	frout=*fr;
	frout.x=NULL;
	snew(frout.x, size);
	int cont = 0;
	for(int i=0; i<(int)v.size(); i++){
		vector<Coordinates> c = v[i];
		for(int j=0; j<(int)c.size(); j++){
			frout.x[cont][XX] = c[j].getX()/FACTOR;
			frout.x[cont][YY] = c[j].getY()/FACTOR;
			frout.x[cont][ZZ] = c[j].getZ()/FACTOR;

			cont ++;
		}
	}

	frout.natoms = size;

	write_trxframe(trjfile, &frout, NULL);

	free(frout.x);
}

int getSolventCount(vector<Atom> waters, Grid grid, Grid& water_statistics, float frame);

int getSolventCount(vector<Atom> waters, Grid grid, vector<vector<Coordinates> > v, Grid& water_statistics, float time){

	int solvent_count = 0;
	Atom sol;
	vector<Coordinates> c;
	Coordinates coor_sol, coor_stat;
	set<int> cavity_indexes;

	int cell_center_x, cell_center_y, cell_center_z = 0;

	float spacing = grid.getSpacing();

	for(unsigned int i=0; i<v.size(); i++){
		cavity_indexes.insert(grid.getGrid()[v[i][0].getGridCoordinate()[0]][v[i][0].getGridCoordinate()[1]][v[i][0].getGridCoordinate()[2]]);
	}

	for(unsigned int k=0; k<waters.size(); k++){

		sol = waters[k];
		coor_sol = sol.getCoordinates();

		cell_center_x = grid.calculatePositionInGrid(coor_sol.getX(), grid.getOriginX());
		cell_center_y = grid.calculatePositionInGrid(coor_sol.getY(), grid.getOriginY());
		cell_center_z = grid.calculatePositionInGrid(coor_sol.getZ(), grid.getOriginZ());

		if(cell_center_x>=0 && cell_center_x<grid.getWidth() && cell_center_y>=0 && cell_center_y<grid.getHeight() && cell_center_z>=0 && cell_center_z<grid.getDepth()){

			if(cavity_indexes.size()>0 && cavity_indexes.find(grid.getGrid()[cell_center_x][cell_center_y][cell_center_z])!=cavity_indexes.end()){

				solvent_count ++;

				if(water_statistics.getTStartStatisitics()<=time){

					int *x_cells, *y_cells, *z_cells, *coor;
					int x, y, z = 0;

					//Apply the atom radius

					if((sol.getVdwRadius())*2>spacing){
						x_cells = grid.getDistanceCells(sol.getVdwRadius(), cell_center_x, grid.getWidth());
						y_cells = grid.getDistanceCells(sol.getVdwRadius(), cell_center_y, grid.getHeight());
						z_cells = grid.getDistanceCells(sol.getVdwRadius(), cell_center_z, grid.getDepth());
					}else{
						x_cells = new int[2]; x_cells[0] = cell_center_x; x_cells[1] = cell_center_x;
						y_cells = new int[2]; y_cells[0] = cell_center_y; y_cells[1] = cell_center_y;
						z_cells = new int[2]; z_cells[0] = cell_center_z; z_cells[1] = cell_center_z;
					}

					for(x=x_cells[0]; x<=x_cells[1];x++){
						for(y=y_cells[0]; y<=y_cells[1];y++){
							for(z=z_cells[0]; z<=z_cells[1];z++){

								int *coor = new int[3]; coor[0] =x; coor[1]=y; coor[2]=z;
								coor_stat = grid.calculateCoordinateInGrid(coor);

								int x = water_statistics.calculatePositionInGrid(coor_stat.getX(),water_statistics.getOriginX());
								int y = water_statistics.calculatePositionInGrid(coor_stat.getY(),water_statistics.getOriginY());
								int z = water_statistics.calculatePositionInGrid(coor_stat.getZ(),water_statistics.getOriginZ());

								if(water_statistics.getWidth()!=0 && (x<water_statistics.getWidth()) && (y>0) && (y<water_statistics.getHeight()) && (z>0) && (z<water_statistics.getDepth())){
									water_statistics.getGrid()[x][y][z] = water_statistics.getGrid()[x][y][z]+1;
								}

								delete [] coor;
							}
						}
					}
					delete[] x_cells;
					delete[] y_cells;
					delete[] z_cells;
				}
			}
		}
	}

	return solvent_count;
}

int updateStatistics(vector<vector<Coordinates> > v, Grid& statistics);

int updateStatistics(vector<vector<Coordinates> > v, Grid& statistics){

	for(unsigned int i=0; i<v.size(); i++){
		vector<Coordinates> c = v[i];
		for(unsigned int j=0; j<c.size(); j++){
			int x = statistics.calculatePositionInGrid(c[j].getX(),statistics.getOriginX());
			int y = statistics.calculatePositionInGrid(c[j].getY(),statistics.getOriginY());
			int z = statistics.calculatePositionInGrid(c[j].getZ(),statistics.getOriginZ());

			if(statistics.getWidth()!=0 && (x<statistics.getWidth()) && (y>0) && (y<statistics.getHeight()) && (z>0) && (z<statistics.getDepth())){
				statistics.getGrid()[x][y][z] = statistics.getGrid()[x][y][z]+1;
			}
		}
	}
	return 0;
}

void updateTunnelSection(map<float,float> section, map<float,pair<vector<float>, float> >& acc);

void updateTunnelSection(map<int,float> section, map<int,pair<vector<float>, float> >& acc){

	//TODO: standard deviation

	map<int,float>::iterator it;

	if(acc.size()==0){
		for(it=section.begin(); it!=section.end(); ++it){
			vector<float> sec (1, section[it->first]);
			acc[it->first] = make_pair(sec, section[it->first]);
		}
	}else{
		for(it=section.begin(); it!=section.end(); ++it){
			if(acc.find(it->first)!=acc.end()){
				acc[it->first].first.push_back(section[it->first]);
				acc[it->first].second = acc[it->first].second + section[it->first];
				//acc[it->first] = make_pair(acc[it->first].first + section[it->first], acc[it->first].second + 1);
			}else{
				vector<float> sec (1, section[it->first]);
				acc[it->first] = make_pair(sec, section[it->first]);
				//acc[it->first] = make_pair(section[it->first], 1);
			}
		}
	}
}

/* This function forces the trajectory to have the same number of atoms than the topology cavity (frame=0) by adding dummy atoms or removing central
 * atoms, depending on the difference. Note that the volume was calculated before, so it will only affect the representation. */

int cavityToTrajectory(t_cavity_options *opt, t_trxframe *fr, vector<vector<Coordinates> >& v, int size, bool tunnel);

int cavityToTrajectory(t_cavity_options *opt, t_trxframe *fr, vector<vector<Coordinates> >& v, int size, bool tunnel){

	int i;
	AtomWriter writer;

	if(opt->params.frame == 0) {

		if(!tunnel){
			opt->params.min_size = size + 1000;
		}else{
			opt->params.tunnel_size = size + 1000;
		}

		for(i=0; i<1000; i++){
			//Dummy coordinates! To ensure we have an initial cavity in case it just appears along the simulation
			Coordinates dummy(v[0][0].getX(), v[0][0].getY(), v[0][0].getZ());
			v[0].push_back(dummy);
		}

		if(!tunnel){
			writer.writePDB(v, opt->files.fcavity);
			fclose(opt->files.fcavity);
		}else{
			writer.writePDB(v, opt->files.ftunnel);
			fclose(opt->files.ftunnel);
		}
	}else{
		int dif = 0;
		if(!tunnel){
			dif = opt->params.min_size - size;
		}else{
			dif = opt->params.tunnel_size - size;
		}

		if(dif!=0){
			if(dif>0){
				for(int i=0; i<dif; i++){
					//Dummy coordinates!
					Coordinates dummy(v[0][0].getX(), v[0][0].getY(), v[0][0].getZ());
					v[0].push_back(dummy);
				}
			}else{
				//Erase the ones in the center of the cavity
				int dif = size - opt->params.min_size;

				if(v.size()==1){
					v[0].erase(v[0].begin(),v[0].begin()+dif);
				}else{
					while(dif>0){
						for(unsigned int i=0; i<v.size();i++){
							if(v[i].size()>2){
								v[i].erase(v[i].begin(),v[i].begin()+1);
								dif--;
								if(dif==0) break;
							}
						}
					}
				}
				if(!tunnel){
					fprintf(opt->files.log, "I have erased %d cavity points so cavity matches cavity topology. Careful!\n", dif);
				}else{
					fprintf(opt->files.log, "I have erased %d tunnel points so tunnel matches tunnel topology. Careful!\n", dif);
				}
				fflush(opt->files.log);
			}
		}
	}

	return 0;
}

int getVolume(Grid grid, t_cavity_options *opt, t_trxframe *fr, vector<Atom> waters);

int getVolume(Grid grid, t_cavity_options *opt, t_trxframe *fr, vector<Atom> waters){

	vector<vector<Coordinates> > v;
	string file;
	AtomWriter writer;
	int index;

	int max_cavity = 0;

	if(opt->use_seed==true){
		Coordinates seed(opt->seed[0], opt->seed[1], opt->seed[2]);
		if(opt->mode==0){
			v = grid.getCavitySeedPoint(seed, opt->params.statistics, index, (int)opt->min, opt->rectify);
			if(v.size()==0 && opt->params.v.size()>0){
				for(unsigned i = 0; i<opt->params.v[0].size();i++){
					v = grid.getCavitySeedPoint(opt->params.v[0][i], opt->params.statistics, index, (int)opt->min, opt->rectify);
					if(v.size()>0) break;
				}
			}
		}else{
			v = grid.getCavitiesSeedPoint(seed, opt->params.statistics, index, max_cavity, (int)opt->min, opt->rectify);
			if(v.size()==0 && opt->params.v.size()>0){
				for(unsigned i = 0; i<opt->params.v[0].size();i++){
					v = grid.getCavitiesSeedPoint(opt->params.v[0][i], opt->params.statistics, index, max_cavity, (int)opt->min, opt->rectify);
					if(v.size()>0) break;
				}
			}
		}
		if(v.size()==0) fprintf(opt->files.log, "No cavity found for coordinates %8.3f %8.3f %8.3f\n", opt->seed[0],opt->seed[1],opt->seed[2]);
	}else{
		vector<vector<Coordinates> > aux = grid.getCavitiesNoSeedPoint(opt->params.statistics, index, max_cavity, (int)opt->min, opt->rectify);
		if(opt->mode==0){
			if(aux.size()>0){
				v.push_back(aux[max_cavity]);
			}
		}else{
			v = aux;
		}
	}

	int size = 0;
	for(unsigned int i=0; i<v.size();i++){
		size = size + v[i].size();
	}

	fprintf(opt->files.log, "Number of cavities %d. Total size of the cavity: %d\n", (int)v.size(), size);

	if(opt->mode == 0){
		writer.writeVolumeFrame(size,opt->spacing, fr->time, opt->files.fvolume);
	}else{
		vector<int> volume_cavities = writer.preprocessVolumeCavities(size, v, opt->params.v, opt->spacing);
		opt->params.volume_cavities.insert(make_pair(fr->time, volume_cavities));
	}

	if(opt->sol){
		int water_count = getSolventCount(waters, grid, v, opt->params.water_statistics, fr->time);
		if(opt->files.fsolcount!=NULL) writer.writeSolventFrame(water_count,opt->spacing, fr->time, opt->files.fsolcount);
	}

	if(fr->time >= opt->params.statistics.getTStartStatisitics()){
		fprintf(opt->files.log, "Calculating statistics... \n");
		updateStatistics(v, opt->params.statistics);
		opt->params.frame_stat = opt->params.frame_stat + 1;
	}

	if(opt->params.frame==0 || opt->files.ftrajectory!=NULL){
		if(size>0){
			cavityToTrajectory(opt, fr, v, size, false);
		}else{
			int* aux = new int[3];
			aux[0] = 0; aux[1] = 0; aux[2] = 0;
			Coordinates dummy = grid.calculateCoordinateInGrid(aux);
			vector<Coordinates> a;
			for(int i=0; i<opt->params.min_size;i ++){
				a.push_back(dummy);
			}
			v.push_back(a);
			fprintf(opt->files.log, "No cavity!\n");
		}
	}

	fflush(opt->files.fvolume);

	if(opt->files.ftrajectory!=NULL){
		writeTrajectory(opt->files.ftrajectory, v, fr, opt->params.min_size);
	}

	// Tunnel

	if(opt->tunnel){
		vector<Coordinates> tunnel;
		float bottleneck = 0.0;

		if(fr->time>opt->t_stat){
			map<int, float> sector = grid.calculateBottleneckArea(index, opt->axis, tunnel, bottleneck, opt->area, true);
			updateTunnelSection(sector, opt->params.acc_sector);
		}else{
			grid.calculateBottleneckArea(index, opt->axis, tunnel, bottleneck, opt->area, false);
		}


		if(opt->files.fbottleneck!=NULL) writer.writeBottleneckAreaFrame(bottleneck, fr->time, opt->files.fbottleneck);
		fprintf(opt->files.log, "Tunnel found. Bottleneck section area %8.3f\n", bottleneck);

		if((opt->files.ftunnel!=NULL && opt->params.frame==0) || opt->files.ftunnel_trajectory!=NULL){
			vector<vector<Coordinates> > tunnels;
			if(tunnel.size()>0){
				tunnels.push_back(tunnel);
				cavityToTrajectory(opt, fr, tunnels, tunnel.size(), true);
			}else{
				int* aux = new int[3];
				aux[0] = 0; aux[1] = 0; aux[2] = 0;
				Coordinates dummy = grid.calculateCoordinateInGrid(aux);
				vector<Coordinates> a;
				for(int i=0; i<opt->params.tunnel_size;i ++){
					tunnel.push_back(dummy);
				}
				tunnels.push_back(a);
				fprintf(opt->files.log, "No tunnel!\n");
			}
			if(opt->files.ftunnel_trajectory!=NULL){
				writeTrajectory(opt->files.ftunnel_trajectory, tunnels, fr, opt->params.tunnel_size);
			}
		}
	}

	fflush(opt->files.log);

	return size;
}

/* This function maps the atoms in the grid and calls for the calculation of the cavity */

int analyze_frame(t_topology * top, t_trxframe * fr, t_pbc * pbc, int nr, gmx_ana_selection_t *sel[], void *data);

int analyze_frame(t_topology * top, t_trxframe * fr, t_pbc * pbc, int nr, gmx_ana_selection_t *sel[], void *data){

	t_cavity_options *opt = (t_cavity_options *)data;
    float limits[6];
	float ca_limits[6];
	int a, i;
	AtomReader reader;
	AtomWriter writer;

	vector<Atom> selection;
	vector<Atom> waters;

	if((opt->params.frame)%opt->skip == 0){

		string atom_n, file;

		a = opt->params.index[0];
		limits[0] = fr->x[a][XX]*FACTOR; limits[1] = fr->x[a][YY]*FACTOR; limits[2] = fr->x[a][ZZ]*FACTOR; limits[3] = fr->x[a][XX]*FACTOR; limits[4] = fr->x[a][YY]*FACTOR; limits[5] = fr->x[a][ZZ]*FACTOR;
		if(opt->ca==true){
			ca_limits[0] = fr->x[a][XX]*FACTOR; ca_limits[1] = fr->x[a][YY]*FACTOR; ca_limits[2] = fr->x[a][ZZ]*FACTOR; ca_limits[3] = fr->x[a][XX]*FACTOR; ca_limits[4] = fr->x[a][YY]*FACTOR; ca_limits[5] = fr->x[a][ZZ]*FACTOR;
		}

		for(i=0;(i<opt->params.isize);i++) {
			Atom atom;
			a = opt->params.index[i];
			atom_n = (string)*(top->atoms.atomname[a]);

			if(!opt->ff_radius /*|| (opt->ff_radius && opt->radii[i]==0)*/){
				Atom aux(top->atoms.atomname[a][0][0],Coordinates(fr->x[a][XX]*FACTOR,fr->x[a][YY]*FACTOR,fr->x[a][ZZ]*FACTOR));
				atom = aux;
			}else{
				Atom aux(top->atoms.atomname[a][0][0],Coordinates(fr->x[a][XX]*FACTOR,fr->x[a][YY]*FACTOR,fr->x[a][ZZ]*FACTOR), opt->params.radii[i]);
				atom = aux;
			}

			//TODO: delete C alpha stuff?
			/*If required, calculate limits of the backbone area*/
			if(opt->ca==true){
				if(atom_n.compare("CA")==0){
					atom.setCa(true);
					reader.calculateLimitsGrid(ca_limits, atom);
				}
			}

			selection.push_back(atom);
			reader.calculateLimitsGrid(limits, atom);
		}

		/* If required, calculate the water grid */
		if(opt->sol==true){
			for(i=0;(i<opt->params.sol_isize);i++) {
				a = opt->params.sol_index[i];
				if((fr->x[a][XX]*FACTOR<=limits[0]) and (fr->x[a][XX]*FACTOR >= limits[3] and (fr->x[a][YY]*FACTOR <= limits[1]) and (fr->x[a][YY]*FACTOR >= limits[4]) and (fr->x[a][ZZ]*FACTOR <= limits[2]) and (fr->x[a][ZZ]*FACTOR >= limits[5])) ){
					Atom atom;
					char atom_n;
					if(top->atoms.atomname[a][0][0]=='N' and top->atoms.atomname[a][0][1]=='A'){
						atom_n = '1';
					}else{
						if(top->atoms.atomname[a][0][0]=='C' and top->atoms.atomname[a][0][1]=='L'){
							atom_n = '2';
						}else{
							atom_n =top->atoms.atomname[a][0][0];
						}
					}

					if(!opt->ff_radius /*|| (opt->ff_radius && opt->sol_radii[i]==0)*/){
						Atom aux(atom_n,Coordinates(fr->x[a][XX]*FACTOR,fr->x[a][YY]*FACTOR,fr->x[a][ZZ]*FACTOR));
						aux.setResId(top->atoms.atom[a].resind);
						atom = aux;
					}else{
						Atom aux(atom_n,Coordinates(fr->x[a][XX]*FACTOR,fr->x[a][YY]*FACTOR,fr->x[a][ZZ]*FACTOR), opt->params.sol_radii[i]);
						atom = aux;
					}
					waters.push_back(atom);
				}
			}
		}

		Grid grid(opt->spacing, limits[0], limits[3], limits[1], limits[4], limits[2], limits[5], opt->cutoff);
		grid.setDimensionsSearch(opt->dim);
		grid.setAtoms(selection);

		if(opt->params.frame == 0){
			Grid stat(opt->spacing, limits[0], limits[3], limits[1], limits[4], limits[2], limits[5], opt->cutoff);
			stat.setTStartStatisitics(opt->t_stat);
			opt->params.statistics = stat;
			if(opt->sol==true){
				Grid wat(opt->spacing, limits[0], limits[3], limits[1], limits[4], limits[2], limits[5], opt->cutoff);
				wat.setTStartStatisitics(opt->t_stat);
				opt->params.water_statistics = wat;
			}
		}

		reader.initializeGrid(grid, selection);

		if(opt->ca==false){
			ca_limits[0] = limits[0]; ca_limits[1] = limits[1];ca_limits[2] = limits[2];ca_limits[3] = limits[3];ca_limits[4] = limits[4];ca_limits[5] = limits[5];
		}

		grid.setCAplhaLimits(ca_limits[0], ca_limits[3], ca_limits[1], ca_limits[4], ca_limits[2], ca_limits[5]);

		time_t rawtime;
		struct tm * timeinfo;

		time (&rawtime);
		timeinfo = localtime(&rawtime);
		string ts = (string)asctime(timeinfo);
		int ch = ts.find_last_of("\n");
		ts.erase(ch, ch+1);

		timeval tp;
		gettimeofday(&tp, 0);

		fprintf(opt->files.log, "\n%s. Frame %d. Time %.2f. Milliseconds %d .\n", ts.c_str(), opt->params.frame, fr->time, (int)tp.tv_usec/1000);
		fprintf(opt->files.log, "Number of atoms protein %d and solvent %d\n", opt->params.isize, opt->params.sol_isize);
		fprintf(opt->files.log, "Limits of the grid: %.2f %.2f %.2f %.2f %.2f %.2f", limits[0], limits[3], limits[1], limits[4], limits[2], limits[5]);
		fprintf(opt->files.log, ".Grid size N = %d \n", grid.getWidth()*grid.getHeight()*grid.getDepth());

		fflush(opt->files.log);

		getVolume(grid, opt, fr, waters);

		grid.deleteGrid();
	}

	opt->params.frame = opt->params.frame + 1;

	return 0;
}

/* I calculate the radii form the parameters of the force field (contained in tpr file)*/

void calculateRadius(t_topology * top, atom_id * selection, int selection_size, float * radii);

void calculateRadius(t_topology * top, atom_id * selection, int selection_size, float * radii){
	int i, itype;
	double c12, c6, sig6;
	int ntype = top->idef.atnr;

	for(i=0; i<selection_size; i++) {
		printf("Atom %d name %s",i, top->atoms.atomname[selection[i]][0]);
		itype = top->atoms.atom[selection[i]].type;
		c12   = top->idef.iparams[itype*ntype+itype].lj.c12;
		c6    = top->idef.iparams[itype*ntype+itype].lj.c6;
		if ((c6 != 0) && (c12 != 0)) {
			sig6 = c12/c6;
			radii[i]= (0.5*pow(sig6,1.0/6.0))*10;  /* Factor of 10 for nm -> Angstroms */
			//printf(" radius %f\n",radii[i]);
		}else{
			//printf("no LJ parameters\n");
			radii[i]=0.0;
		}
	}
}

void AtomReader::start(int argc,char *argv[]){
	t_cavity_options opt;
	gmx_ana_traj_t *trj;
	output_env_t oenv;
	const char *tmp_fnm;
	AtomWriter writer;

	opt.files.log = ffopen("log.txt", "w");
	for(int i=0;i<argc;i++) fprintf(opt.files.log, "%s ", argv[i]);
	fprintf(opt.files.log, "\n");

	char **grpname; /* the name of each group */

	const char *desc[] = {
		"This is an analysis program that implements cavity volume calculation",
	};

	const char *mode[4]={NULL, "all","max", NULL};
	const char *orientation[5]={NULL, "z","x","y", NULL};

	t_filenm fnm[] = {
			{ efTRX, NULL, NULL, ffOPTRD },
			{ efTPS, NULL, NULL, ffREAD },
			{ efNDX, NULL, NULL, ffOPTRD },
			{ efPDB, "-o", "cavity", ffWRITE },
			{ efXVG, "-ov", "volume", ffWRITE },
			{ efTRO, "-ot", "cav_traj",ffOPTWR },
			{ efPDB, "-ostat", "stat", ffOPTWR },
			{ efXVG, "-os", "sol_count", ffOPTWR },
			{ efPDB, "-osstat","sol_stat",ffOPTWR },
			{ efPDB, "-ob", "tunnel", ffOPTWR },
			{ efTRO, "-obt", "tun_traj",ffOPTWR },
			{ efXVG, "-oba", "min_area", ffOPTWR },
			{ efXVG, "-obs", "section_area", ffOPTWR },
	};

	t_pargs pa[] = {
			{ "-skip", NULL, etINT, { &opt.skip }, "Only write every nr-th frame" },
			{ "-mode", FALSE , etENUM, { mode }, "Cavities to output: " },
			{ "-min", NULL, etREAL, { &opt.min }, "Min. size of cavities (A^3)" },
			{ "-dim", NULL, etINT, { &opt.dim }, "Number of dimension to be surrounded by protein" },
			{ "-spacing", NULL, etREAL, {&opt.spacing},  "Grid spacing"},
			{ "-tstat", NULL, etINT, {&opt.t_stat},  "Time (ps) to start the calculation of statistics"},
			{ "-rectify", NULL, etBOOL, {&opt.rectify},  "Rectify output based on statistics"},
			{ "-ff_radius", NULL, etBOOL, {&opt.ff_radius},  "Use force field radii (only using tpr topology)"},
			{ "-axis", FALSE , etENUM, { orientation },  "Orientation of the protein for tunnel calculation"},
			{ "-ca", NULL, etBOOL, {&opt.ca},  "Only cavities enclosed in the backbone"},
			{ "-seed", FALSE, etRVEC,{&opt.seed}, "Coordinates of the seed point"},
			{ "-area", FALSE, etREAL ,{&opt.area}, "Axis value to be inspected"},
			{ "-cutoff", FALSE, etREAL ,{&opt.cutoff}, "Cut-off (test)"}
	};

	opt.dim = 5;
	opt.mode = 1;
	opt.skip = 1;
	opt.min  = 50;
	opt.spacing = 1.4;
	opt.t_stat = 0.0;
	opt.seed[0] = 0.0;
	opt.seed[1] = 0.0;
	opt.seed[2] = 0.0;
	opt.axis = 2;
	opt.area = NULL;
	opt.cutoff = 0;

	opt.params.min_size = 0;
	opt.params.tunnel_size = 0;
	opt.params.frame_stat = 0;
	opt.params.frame = 0;
	opt.rectify = false;
	opt.ff_radius = false;
	opt.sol = false;
	opt.ca = false;
	opt.tunnel = false;
	opt.use_seed = false;

	opt.files.fcavity = NULL;
	opt.files.ftrajectory = NULL;
	opt.files.fvolume = NULL;
	opt.files.fsolcount = NULL;
	opt.files.ftunnel = NULL;
	opt.files.fbottleneck = NULL;
	opt.files.ftunnel_trajectory = NULL;

	#define NFILE asize(fnm)
	#define NPA asize(pa)
	#define FLAGS (TRX_READ_X | ANA_USE_TOPX)

	gmx_ana_traj_create(&trj, FLAGS);
	gmx_ana_set_nrefgrps(trj, 1);

	parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_BE_NICE, NFILE ,fnm,asize(pa),pa,asize(desc),desc,0,NULL, &oenv);
	trj->oenv = oenv;

	/* Process our own options.
	 * Make copies of file names for easier memory management. */
	tmp_fnm            = ftp2fn_null(efTRX, NFILE, fnm);
	trj->trjfile         = tmp_fnm ? strdup(tmp_fnm) : NULL;
	tmp_fnm            = ftp2fn_null(efTPS, NFILE, fnm);
	trj->topfile         = tmp_fnm ? strdup(tmp_fnm) : NULL;
	trj->topfile_notnull = strdup(ftp2fn(efTPS, NFILE, fnm));
	tmp_fnm            = ftp2fn_null(efNDX, NFILE, fnm);
	trj->ndxfile         = tmp_fnm ? strdup(tmp_fnm) : NULL;

	if(nenum(mode)>1) opt.mode = 0;
	int axis = nenum(orientation);
	if(axis>2) opt.axis = axis-2;

	/* If no trajectory file is given, we need to set some flags to be able
	 * to prepare a frame from the loaded topology information. Also, check
	 * that a topology is provided. */

	if (!trj->trjfile)
	{
		if (!trj->topfile)
		{
			gmx_input("No trajectory or topology provided, nothing to do!");
		}
		trj->flags |= ANA_REQUIRE_TOP;
		trj->flags |= ANA_USE_TOPX;
	}else{
		if(opt2fn_null("-ot",NFILE,fnm)) opt.files.ftrajectory = open_trx(opt2fn("-ot", NFILE, fnm),"w");
		if(opt2fn_null("-obt",NFILE,fnm)) opt.files.ftunnel_trajectory = open_trx(opt2fn("-obt", NFILE, fnm),"w");
	}

	/* Read topology */

	snew(trj->top, 1);
	char title[STRLEN];
	trj->bTop = read_tps_conf(trj->topfile_notnull, title, trj->top, &trj->ePBC, &trj->xtop, NULL, trj->boxtop, FALSE);

	/* Read index file */

	if(ftp2bSet(efNDX,NFILE,fnm)){
		printf("Select group of atoms to calculate the cavities within:\n");
		snew(grpname,1);
		get_index(&trj->top->atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&opt.params.isize,&opt.params.index,grpname);
	}else{
		int i;
		opt.params.isize = trj->top->atoms.nr;
		snew(opt.params.index,opt.params.isize);
		for(i=0; i<opt.params.isize; i++){
			opt.params.index[i] = i;
		}
	}

	/* If solvent calculation required, identify solvent atoms*/

	if (opt2fn_null("-os",NFILE,fnm) || opt2fn_null("-osstat",NFILE,fnm)){
		opt.sol=true;

		if(ftp2bSet(efNDX,NFILE,fnm)){
			printf("Select solvent:\n");
			snew(grpname,1);
			get_index(&trj->top->atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&opt.params.sol_isize,&opt.params.sol_index,grpname);
		}else{
			int cont=0;
			opt.params.sol_isize = trj->top->atoms.nr;
			snew(opt.params.index,opt.params.sol_isize);
			for(int i=0; i<trj->top->atoms.nr; i++) {
				if(strcmp(trj->top->atoms.resinfo[trj->top->atoms.atom[i].resind].name[0], "SOL") == 0){
					opt.params.sol_index[cont]=i;
					cont++;
				}
			}
			opt.params.sol_isize = cont;
		}
	}

	if (opt2parg_bSet("-seed",NPA,pa)){
		opt.use_seed = true;
	}

	if(opt2fn_null("-ob",NFILE,fnm) || opt2fn_null("-obt",NFILE,fnm) || opt2fn_null("-oba",NFILE,fnm) || opt2fn_null("-obs",NFILE,fnm)){
		opt.tunnel = true;
	}

	/* Get force field radii */
	if(opt.ff_radius==true){
		if(ftp2bSet(efTPR,NFILE,fnm)){
			snew(opt.params.radii,opt.params.isize);
			calculateRadius(trj->top, opt.params.index, opt.params.isize, opt.params.radii);
			if(opt.sol==true){
				snew(opt.params.sol_radii,opt.params.sol_isize);
				calculateRadius(trj->top, opt.params.sol_index, opt.params.sol_isize, opt.params.sol_radii);
			}
		}else{
			opt.ff_radius = false;
			printf("Option -ff_radius requires a .tpr file for the -s option\n" );
		}
	}

	/* Min number of grids of cavity */
	opt.min = (int)((opt.min)/(opt.spacing*opt.spacing*opt.spacing) + 0.5);

	/* Open output files: default files */

	opt.files.fvolume = xvgropen(opt2fn("-ov", NFILE, fnm), "Volume evolution", "Time [ps]", "Volume [angstroms\\S3\\N]", oenv);
	opt.files.fcavity =ffopen(opt2fn("-o", NFILE, fnm),"w");

	/* Open output files: optional files */

	if(opt2fn_null("-os",NFILE,fnm)) opt.files.fsolcount =xvgropen(opt2fn("-os", NFILE, fnm), "Solvent Atoms", "Time [ps]", "Solvent atoms", oenv);
	if(opt2fn_null("-oba",NFILE,fnm)) opt.files.fbottleneck =xvgropen(opt2fn("-oba", NFILE, fnm), "Min. tunnel section area", "Time [ps]", "Area [angstroms\\S2\\N]", oenv);
	if(opt2fn_null("-ob",NFILE,fnm) || opt2fn_null("-obt",NFILE,fnm)) opt.files.ftunnel =ffopen(opt2fn("-ob", NFILE, fnm),"w");

	/* Do the actual analysis*/

	gmx_ana_do(trj, 0, &analyze_frame, &opt);

	/* summary results and close files */

	if(opt.mode == 1){
		writer.writeVolumeAllFrames(opt.params.volume_cavities, opt.spacing, opt.files.fvolume);
	}

	if (trj->trjfile){
		if(opt2fn_null("-ostat",NFILE,fnm)){
			fprintf(opt.files.log, "Statistics calculated for the last %d frames\n", opt.params.frame_stat);
			fflush(opt.files.log);
			FILE* fstatistics = ffopen(opt2fn("-ostat", NFILE, fnm),"w");
			writer.writeStatistics(opt.params.statistics, opt.spacing, opt.params.frame_stat ,fstatistics);
			if(opt2fn_null("-ot",NFILE,fnm)) close_trx(opt.files.ftrajectory);
		}
		if(opt2fn_null("-obt",NFILE,fnm)) close_trx(opt.files.ftunnel_trajectory);
	}

	if(opt.sol == true){
		if (trj->trjfile && opt2fn_null("-osstat",NFILE,fnm)){
			FILE* fwater_statistics = ffopen(opt2fn("-osstat", NFILE, fnm),"w");
			writer.writeStatistics(opt.params.water_statistics, opt.spacing, opt.params.frame_stat ,fwater_statistics);
		}
		if(opt2fn_null("-os",NFILE,fnm)) fclose(opt.files.fsolcount);
	}

	fclose(opt.files.fvolume);
	fclose(opt.files.log);

	if(opt.tunnel == true){
		if(opt2fn_null("-obs",NFILE,fnm)){
			const char *leyend;
			if(opt.axis==0){
				leyend="X [angstroms]";
			}else{
				if(opt.axis==1){
					leyend="Y [angstroms]";
				}else{
					leyend="Z [angstroms]";
				}
			}

			FILE* fsection =xvgropen(opt2fn("-obs", NFILE, fnm), "Tunnel section area", leyend, "Area [angstroms\\S2\\N]", oenv);
			writer.writeSectionAreaFrame(opt.params.acc_sector, opt.params.frame_stat, fsection);
			fclose(fsection);
		}
		if(opt2fn_null("-oba",NFILE,fnm)) fclose(opt.files.fbottleneck);
	}

	/* Free memory */

	gmx_ana_traj_free(trj);

	free(opt.params.index);
	free(opt.params.sol_index);
	opt.params.statistics.deleteGrid();
	opt.params.water_statistics.deleteGrid();

	printf("end");
}




