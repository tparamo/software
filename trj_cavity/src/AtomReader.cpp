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


struct t_cavity_options
{
	int mode; /* Mode of execution*/
	float min; /* Minimun size of cavities */
	int dim; /* Number of dimensions to cover by the protein */
	int skip; /* Frames to skip when processing trajectory */
	float spacing; /* Grid spacing */
	int t_stat; /* Time to strat calculating statistics */
	bool rectify; /* Rectify output?*/
	bool ff_radius; /* use force field radius*/
	bool water; /* water monitoring */
	bool ca; /* Restrict to ca*/
	bool tunnel; /* water monitoring */
	rvec seed; /* seed point */
	int axis; /* Tunnel orientation */

    //string fpath; /* where the files might be saved*/
	FILE *fvolume; /* volume file*/
	FILE *fcavity; /* cavity file*/
	FILE *ftunnel; /* tunnel file */
	FILE *fwatercount; /* water count file*/
	FILE *fbottleneck; /* water count file*/
	FILE *log; /* log file */
	t_trxstatus *ftrajectory; /* cavity trajectory file*/
	t_trxstatus *ftunnel_trajectory; /* tunnel trajectory file*/

	atom_id *index;   /* the index for the atom numbers */
	int isize;    /* the size of each group */
	float *radii; /* ff radii */
	atom_id *sol_index;   /* the index for the solvent  */
	int sol_isize;    /* the size of each group */
	float *sol_radii; /* ff radii for solvent */
    Grid statistics; /* statistics*/
    Grid water_statistics; /* water statistics */
    int min_size; /* cavity min size*/
    int tunnel_size; /* tunnel min size*/
    int frame_stat;
    int frame;
    bool use_seed;
    vector<vector<Coordinates> > v; /* First set of cavities*/
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

	int x, y, z;

	for(unsigned int i=0;i<atoms.size();i++){

		cell_center_x = grid.calculatePositionInGrid(atoms[i].getCoordinates().getX(), grid.getOriginX());
		cell_center_y = grid.calculatePositionInGrid(atoms[i].getCoordinates().getY(), grid.getOriginY());
		cell_center_z = grid.calculatePositionInGrid(atoms[i].getCoordinates().getZ(), grid.getOriginZ());

		//Apply the atom radius

		if((atoms[i].getVdwRadius())*2>spacing){
			x_cells = distanceCells(spacing, atoms[i].getVdwRadius(), cell_center_x, grid.getWidth());
			y_cells = distanceCells(spacing, atoms[i].getVdwRadius(), cell_center_y, grid.getHeight());
			z_cells = distanceCells(spacing, atoms[i].getVdwRadius(), cell_center_z, grid.getDepth());
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

int* AtomReader::distanceCells(float spacing, float radius, int center, int limit){
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
}

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
}

int getWaterCount(vector<Atom> waters, vector<vector<Coordinates> > v, Grid& water_statistics, int frame);

int getWaterCount(vector<Atom> waters, vector<vector<Coordinates> > v, Grid& water_statistics, int time){
	//this part is for water monitoring ->make it fancier!!
	int water_count = 0;
	bool water_flag = false;
	int water_molecule = 0;

	for(unsigned int k=0; k<waters.size(); k++){
		if(water_molecule != waters[k].getResId()){
			for(unsigned int i=0; i<v.size(); i++){
				vector<Coordinates> c = v[i];
				for(unsigned int j=0; j<c.size(); j++){
					if((waters[k].getCoordinates().getX()<=c[j].getX()+water_statistics.getSpacing()/2) and (waters[k].getCoordinates().getX()>c[j].getX()-water_statistics.getSpacing()/2)) {
						if((waters[k].getCoordinates().getY()<=c[j].getY()+water_statistics.getSpacing()/2) and (waters[k].getCoordinates().getY()>c[j].getY()-water_statistics.getSpacing()/2)) {
							if((waters[k].getCoordinates().getZ()<=c[j].getZ()+water_statistics.getSpacing()/2) and (waters[k].getCoordinates().getZ()>c[j].getZ()-water_statistics.getSpacing()/2)) {
								water_count ++;

								if(water_statistics.getTStartStatisitics()<=time){
									int x = water_statistics.calculatePositionInGrid(c[j].getX(),water_statistics.getOriginX());
									int y = water_statistics.calculatePositionInGrid(c[j].getY(),water_statistics.getOriginY());
									int z = water_statistics.calculatePositionInGrid(c[j].getZ(),water_statistics.getOriginZ());

									if(water_statistics.getWidth()!=0 && (x<water_statistics.getWidth()) && (y>0) && (y<water_statistics.getHeight()) && (z>0) && (z<water_statistics.getDepth())){
										water_statistics.getGrid()[x][y][z] = water_statistics.getGrid()[x][y][z]+1;
									}
								}
								water_flag = true;
								water_molecule = waters[k].getResId();
								break;
							}
						}
					}
					if(water_flag==true){
						water_flag = false;
						break;
					}
				}
			}
		}
	}

	return water_count;
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

int cavityToTrajectory(t_cavity_options *opt, t_trxframe *fr, vector<vector<Coordinates> >& v, int size, bool tunnel);

int cavityToTrajectory(t_cavity_options *opt, t_trxframe *fr, vector<vector<Coordinates> >& v, int size, bool tunnel){
	int i;
	AtomWriter writer;

	if(opt->frame == 0) {
		if(!tunnel){
			opt->min_size = size + 1000;
		}else{
			opt->tunnel_size = size + 1000;
		}
		for(i=0; i<1000; i++){
			//Dummy coordinates! To ensure we have an initial cavity in case it just appears along the simulation
			Coordinates dummy = v[0][0];
			v[0].push_back(dummy);
		}

		if(!tunnel){
			writer.writePDB(v, opt->fcavity);
			fclose(opt->fcavity);
		}else{
			writer.writePDB(v, opt->ftunnel);
			fclose(opt->ftunnel);
		}
	}else{
		int dif = 0;
		if(!tunnel){
			dif = opt->min_size - size;
		}else{
			dif = opt->tunnel_size - size;
		}

		if(dif!=0){
			if(dif>0){
				for(int i=0; i<dif; i++){
					//Dummy coordinates!
					Coordinates dummy = v[0][0];
					v[0].push_back(dummy);
				}
			}else{
				//Erase the ones in the center of the cavity
				int dif = size - opt->min_size;

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
					fprintf(opt->log, "I have erased %d cavity points so cavity matches cavity topology. Careful!\n", dif);
				}else{
					fprintf(opt->log, "I have erased %d tunnel points so tunnel matches tunnel topology. Careful!\n", dif);
				}
				fflush(opt->log);
			}
		}

		/*if(size!=opt->min_size){
			if(size<opt->min_size){
				int dif = opt->min_size - size;

				for(int i=0; i<dif; i++){
					//Dummy coordinates!
					Coordinates dummy = v[0][0];
					v[0].push_back(dummy);
				}
			}else{
				//Erase the ones in the center of the cavity
				int dif = size - opt->min_size;

				if(v.size()==1){
					v[0].erase(v[0].begin(),v[0].begin()+dif);
				}else{
					while(dif>0){
						for(int i=0; i<v.size();i++){
							if(v[i].size()>2){
								v[i].erase(v[i].begin(),v[i].begin()+1);
								dif--;
								if(dif==0) break;
							}
						}
					}
				}
				fprintf(opt->log, "I have erased %d cavity points so cavity matches topology. Careful!\n", dif);
				fflush(opt->log);
			}
		}*/
	}

	return 0;
}

int getVolume(Grid grid, t_cavity_options *opt, t_trxframe *fr, vector<Atom> waters);

int getVolume(Grid grid, t_cavity_options *opt, t_trxframe *fr, vector<Atom> waters){

	vector<vector<Coordinates> > v;
	string file;
	AtomWriter writer;
	int index;

	if(opt->use_seed==true){
		Coordinates seed(opt->seed[0], opt->seed[1], opt->seed[2]);
		v = grid.getCavitySeedPoint(seed, opt->statistics, index, opt->rectify);
		/*if(fr->time>0.0 && v.size()==0){
			for(unsigned i = 0; i<opt->v[0].size();i++){
				v = grid.getCavitySeedPoint(opt->v[0][i], opt->statistics, opt->rectify);
				if(v.size()>0) break;
			}
		}*/
		if(v.size()==0) fprintf(opt->log, "No cavity found for coordinates %8.3f %8.3f %8.3f\n", opt->seed[0],opt->seed[1],opt->seed[2]);
	}else{
		int max_cavity = 0;
		vector<vector<Coordinates> > aux = grid.getCavitiesNoSeedPoint(opt->statistics, index, max_cavity, opt->rectify);
		if(opt->mode==0){
			if(aux.size()>0){
				v.push_back(aux[max_cavity]);
			}
		}else{
			v = aux;
		}
	}

	if(opt->frame==0) opt->v = v;

	int size = 0;
	for(unsigned int i=0; i<v.size();i++){
		size = size + v[i].size();
	}

	fprintf(opt->log, "Number of cavities %d. Total size of the cavity: %d\n", (int)v.size(), size);

	if(opt->mode == 0){
		writer.writeVolumeFrame(size,opt->spacing, fr->time, opt->fvolume);
	}else{
		writer.writeVolumeFrame(v, opt->v, opt->spacing, fr->time, opt->fvolume);
	}

	if(opt->water == true){
		int water_count = getWaterCount(waters, v, opt->water_statistics, (int)fr->time);
		writer.writeWatersFrame(water_count,opt->spacing,(int)fr->time, opt->fwatercount);
	}

	if((int)fr->time >= opt->statistics.getTStartStatisitics()){
		fprintf(opt->log, "Calculating statistics... \n");
		updateStatistics(v, opt->statistics);
		opt->frame_stat = opt->frame_stat + 1;
	}

	if(size>0){
		cavityToTrajectory(opt, fr, v, size, false);
	}else{
		int* aux = new int[3];
		aux[0] = 0; aux[1] = 0; aux[2] = 0;
		Coordinates dummy = grid.calculateCoordinateInGrid(aux);
		vector<Coordinates> a;
		for(int i=0; i<opt->min_size;i ++){
			a.push_back(dummy);
		}
		v.push_back(a);
		fprintf(opt->log, "No cavity!\n");
	}

	fflush(opt->fvolume);

	if(opt->ftrajectory){
		writeTrajectory(opt->ftrajectory, v, fr, opt->min_size);
	}

	// Tunnel

	if(opt->tunnel){
		vector<Coordinates> tunnel;
		float bottleneck = grid.calculateBottleneckArea(index, tunnel);
		writer.writeSectionAreaFrame(bottleneck, (int)fr->time, opt->fbottleneck);
		fprintf(opt->log, "Tunnel found. Section area %8.3f\n", bottleneck);

		vector<vector<Coordinates> > tunnels;
		if(tunnel.size()>0){
			tunnels.push_back(tunnel);
			cavityToTrajectory(opt, fr, tunnels, tunnel.size(), true);
		}else{
			int* aux = new int[3];
			aux[0] = 0; aux[1] = 0; aux[2] = 0;
			Coordinates dummy = grid.calculateCoordinateInGrid(aux);
			vector<Coordinates> a;
			for(int i=0; i<opt->tunnel_size;i ++){
				tunnel.push_back(dummy);
			}
			tunnels.push_back(a);
			fprintf(opt->log, "No tunnel!\n");
		}
		if(opt->ftunnel_trajectory){
			writeTrajectory(opt->ftunnel_trajectory, tunnels, fr, opt->tunnel_size);
		}
	}

	fflush(opt->log);

	return size;
}

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

	if((opt->frame)%opt->skip == 0){

		string atom_n, file;

		a = opt->index[0];
		limits[0] = fr->x[a][XX]*FACTOR; limits[1] = fr->x[a][YY]*FACTOR; limits[2] = fr->x[a][ZZ]*FACTOR; limits[3] = fr->x[a][XX]*FACTOR; limits[4] = fr->x[a][YY]*FACTOR; limits[5] = fr->x[a][ZZ]*FACTOR;
		if(opt->ca==true){
			ca_limits[0] = fr->x[a][XX]*FACTOR; ca_limits[1] = fr->x[a][YY]*FACTOR; ca_limits[2] = fr->x[a][ZZ]*FACTOR; ca_limits[3] = fr->x[a][XX]*FACTOR; ca_limits[4] = fr->x[a][YY]*FACTOR; ca_limits[5] = fr->x[a][ZZ]*FACTOR;
		}

		fprintf(opt->log, "Number of atoms protein %d and solvent %d\n", opt->isize, opt->sol_isize);
		fflush(opt->log);

		for(i=0;(i<opt->isize);i++) {
			Atom atom;
			a = opt->index[i];
			atom_n = (string)*(top->atoms.atomname[a]);

			if(!opt->ff_radius /*|| (opt->ff_radius && opt->radii[i]==0)*/){
				Atom aux(top->atoms.atomname[a][0][0],Coordinates(fr->x[a][XX]*FACTOR,fr->x[a][YY]*FACTOR,fr->x[a][ZZ]*FACTOR));
				atom = aux;
			}else{
				Atom aux(top->atoms.atomname[a][0][0],Coordinates(fr->x[a][XX]*FACTOR,fr->x[a][YY]*FACTOR,fr->x[a][ZZ]*FACTOR), opt->radii[i]);
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
		if(opt->water==true){
			for(i=0;(i<opt->sol_isize);i++) {
				a = opt->sol_index[i];
				if((fr->x[a][XX]*FACTOR<=limits[0]) and (fr->x[a][XX]*FACTOR >= limits[3] and (fr->x[a][YY]*FACTOR <= limits[1]) and (fr->x[a][YY]*FACTOR >= limits[4]) and (fr->x[a][ZZ]*FACTOR <= limits[2]) and (fr->x[a][ZZ]*FACTOR >= limits[5])) ){
					Atom atom;
					if(!opt->ff_radius /*|| (opt->ff_radius && opt->sol_radii[i]==0)*/){
						Atom aux(top->atoms.atomname[a][0][0],Coordinates(fr->x[a][XX]*FACTOR,fr->x[a][YY]*FACTOR,fr->x[a][ZZ]*FACTOR),top->atoms.atom[a].resind);
						atom = aux;
					}else{
						Atom aux(top->atoms.atomname[a][0][0],Coordinates(fr->x[a][XX]*FACTOR,fr->x[a][YY]*FACTOR,fr->x[a][ZZ]*FACTOR), opt->sol_radii[i]);
						atom = aux;
					}
					waters.push_back(atom);
				}
			}
		}

		time_t rawtime;
		struct tm * timeinfo;

		time (&rawtime);
		timeinfo = localtime(&rawtime);
		string ts = (string)asctime(timeinfo);
		int ch = ts.find_last_of("\n");
		ts.erase(ch, ch+1);

		Grid grid(opt->spacing, limits[0], limits[3], limits[1], limits[4], limits[2], limits[5]);
		grid.setDimensionsSearch(opt->dim);
		grid.setAtoms(selection);

		if(opt->frame == 0){
			Grid stat(opt->spacing, limits[0], limits[3], limits[1], limits[4], limits[2], limits[5]);
			stat.setTStartStatisitics(opt->t_stat);
			opt->statistics = stat;
			if(opt->water==true){
				Grid wat(opt->spacing, limits[0], limits[3], limits[1], limits[4], limits[2], limits[5]);
				wat.setTStartStatisitics(opt->t_stat);
				opt->water_statistics = wat;
			}
		}

		reader.initializeGrid(grid, selection);

		if(opt->ca==false){
			ca_limits[0] = limits[0]; ca_limits[1] = limits[1];ca_limits[2] = limits[2];ca_limits[3] = limits[3];ca_limits[4] = limits[4];ca_limits[5] = limits[5];
		}

		grid.setCAplhaLimits(ca_limits[0], ca_limits[3], ca_limits[1], ca_limits[4], ca_limits[2], ca_limits[5]);

		fprintf(opt->log, "\n%s: Frame %d. Time %.2f. Atoms found: %d .\n", ts.c_str(), opt->frame, fr->time, fr->natoms);
		fprintf(opt->log, "Limits of the grid: %.2f %.2f %.2f %.2f %.2f %.2f", limits[0], limits[3], limits[1], limits[4], limits[2], limits[5]);
		fprintf(opt->log, ". Limits of the c-alpha grid: %.2f %.2f %.2f %.2f %.2f %.2f \n", ca_limits[0], ca_limits[3], ca_limits[1], ca_limits[4], ca_limits[2], ca_limits[5]);
		fflush(opt->log);

		getVolume(grid, opt, fr, waters);

		grid.deleteGrid();

	}

	opt->frame = opt->frame + 1;

	return 0;
}

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
			printf(" radius %f\n",radii[i]);
		}else{
			printf("no LJ parameters\n");
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

	char **grpname; /* the name of each group */

	const char *desc[] = {
		"This is an analysis program that implements cavity volume calculation",
	};

	const char *mode[4]={NULL, "all","max", NULL};
	const char *orientation[5]={NULL, "z","x","y", NULL};

	t_filenm fnm[] = {
			{ efTRX, NULL,   NULL,      ffOPTRD },
			{ efTPS, NULL,   NULL,      ffREAD },
			{ efNDX, NULL,   NULL,      ffOPTRD },
			{ efPDB, "-o",   "cavity",  ffWRITE },
			{ efPDB, "-otunnel",   "tunnel",  ffWRITE },
			{ efXVG, "-ov",  "volume",  ffWRITE },
			{ efXVG, "-ow",  "water",  ffOPTWR },
			{ efXVG, "-ob",  "bottleneck_area",  ffOPTWR },
			{ efTRO, "-ot",  "cav_traj",ffOPTWR },
			{ efPDB, "-ostat","statistics",ffOPTWR },
			{ efPDB, "-owstat","wstatistics",ffOPTWR },
			{ efTRO, "-ott",  "tun_traj",ffOPTWR },
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
			{ "-sol", NULL, etBOOL, {&opt.water},  "Perform solvent monitoring"},
			{ "-tunnel", NULL, etBOOL, {&opt.tunnel},  "Perform tunnel monitoring"},
			{ "-axis", FALSE , etENUM, { orientation },  "Orientation of the protein for tunnel calculation"},
			{ "-ca", NULL, etBOOL, {&opt.ca},  "Only cavities enclosed in the backbone"},
			{ "-seed", FALSE, etRVEC,{&opt.seed}, "Coordinates of the seed point"},
	};

	opt.dim = 5;
	opt.mode = 1;
	opt.skip = 1;
	opt.min  = 100;
	opt.spacing = 1.4;
	opt.t_stat = 0;
	opt.rectify = false;
	opt.ff_radius = false;
	opt.water = false;
	opt.min_size = 0;
	opt.frame_stat = 0;
	opt.frame = 0;
	opt.ca = false;
	opt.tunnel = false;
	opt.tunnel_size = 0;
	opt.use_seed = false;
	opt.seed[0] = 0.0;
	opt.seed[1] = 0.0;
	opt.seed[2] = 0.0;
	opt.axis = 2;


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
		opt.ftrajectory = open_trx(opt2fn("-ot", NFILE, fnm),"w");
		if(opt.tunnel) opt.ftunnel_trajectory = open_trx(opt2fn("-ott", NFILE, fnm),"w");
	}

	/* Read topology */

	snew(trj->top, 1);
	char title[STRLEN];
	trj->bTop = read_tps_conf(trj->topfile_notnull, title, trj->top, &trj->ePBC, &trj->xtop, NULL, trj->boxtop, FALSE);

	/* Read index file */

	if(ftp2bSet(efNDX,NFILE,fnm)){
		printf("Select group of atoms to calculate the cavities within:\n");
		snew(grpname,1);
		get_index(&trj->top->atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&opt.isize,&opt.index,grpname);
	}else{
		int i;
		opt.isize = trj->top->atoms.nr;
		snew(opt.index,opt.isize);
		for(i=0; i<opt.isize; i++){
			opt.index[i] = i;
		}
	}

	/* Identify solvent */

	if(opt.water==true){
		if(ftp2bSet(efNDX,NFILE,fnm)){
			printf("Select solvent:\n");
			snew(grpname,1);
			get_index(&trj->top->atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&opt.sol_isize,&opt.sol_index,grpname);
		}else{
			int cont=0;
			opt.sol_isize = trj->top->atoms.nr;
			snew(opt.index,opt.sol_isize);
			for(int i=0; i<trj->top->atoms.nr; i++) {
				if(strcmp(trj->top->atoms.resinfo[trj->top->atoms.atom[i].resind].name[0], "SOL") == 0){
					opt.sol_index[cont]=i;
					cont++;
				}
			}
			opt.sol_isize = cont;
		}
	}

	if (opt2parg_bSet("-seed",NPA,pa)){
		opt.use_seed = true;
	}

	/* Get force field radii */
	if(opt.ff_radius==true){
		if(ftp2bSet(efTPR,NFILE,fnm)){
			snew(opt.radii,opt.isize);
			calculateRadius(trj->top, opt.index, opt.isize, opt.radii);
			if(opt.water==true){
				snew(opt.sol_radii,opt.sol_isize);
				calculateRadius(trj->top, opt.sol_index, opt.sol_isize, opt.sol_radii);
			}
		}else{
			opt.ff_radius = false;
			printf("Option -ff_radius requires a .tpr file for the -s option\n" );
		}
	}

	/* Open output files */

	opt.fvolume = xvgropen(opt2fn("-ov", NFILE, fnm), "Volume evolution", "Time [ps]", "Volume [A^3]", oenv);
	opt.fcavity =ffopen(opt2fn("-o", NFILE, fnm),"w");
	opt.log = ffopen("log.txt", "w");
	if(opt.water == true){
		opt.fwatercount =xvgropen(opt2fn("-ow", NFILE, fnm), "Water Molecules Count", "Time [ps]", "Water molecules", oenv);
	}
	if(opt.tunnel == true){
		opt.ftunnel =ffopen(opt2fn("-otunnel", NFILE, fnm),"w");
		opt.fbottleneck =xvgropen(opt2fn("-ob", NFILE, fnm), "Min. tunnel section area", "Time [ps]", "Area [A^2]", oenv);
	}

	/* Do the actual analysis*/

	gmx_ana_do(trj, 0, &analyze_frame, &opt);

	/* summary results and close files */

	if (trj->trjfile){
		fprintf(opt.log, "Statistics calculated for the last %d frames\n", opt.frame_stat);
		fflush(opt.log);
		FILE* fstatistics = ffopen(opt2fn("-ostat", NFILE, fnm),"w");
		writer.writeStatistics(opt.statistics, opt.spacing, opt.frame_stat ,fstatistics);
		close_trx(opt.ftrajectory);
		if(opt.tunnel) close_trx(opt.ftunnel_trajectory);
	}

	if(opt.water == true){
		if (trj->trjfile){
			FILE* fwater_statistics = ffopen(opt2fn("-owstat", NFILE, fnm),"w");
			writer.writeStatistics(opt.water_statistics, opt.spacing, opt.frame_stat ,fwater_statistics);
		}
		fclose(opt.fwatercount);
	}

	fclose(opt.fvolume);
	fclose(opt.log);

	if(opt.tunnel == true){
		fclose(opt.fbottleneck);
	}

	/* Free memory */

	gmx_ana_traj_free(trj);

	free(opt.index);
	free(opt.sol_index);
	opt.statistics.deleteGrid();
	opt.water_statistics.deleteGrid();

	printf("end");
}




