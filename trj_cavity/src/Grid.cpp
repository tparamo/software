/*
 * Grid.cpp
 *
 *  Created on: Nov 19, 2010
 *      Author: tp334
 */

#include "../include/Grid.h"
#include "../include/AtomWriter.h"

#include <iostream>

#define Pi 3.1416

using namespace std;

Grid::Grid() {
	// TODO Auto-generated constructor stub

}

Grid::Grid(float spacingGrid,float x_max, float x_min, float y_max, float y_min, float z_max, float z_min, float cutoff) {
	spacing = spacingGrid;

	originX = x_min;
	originY = y_min;
	originZ = z_min;

	width = calculateNumberCells(x_max, x_min);
	height = calculateNumberCells(y_max, y_min);
	depth = calculateNumberCells(z_max, z_min);

	tStartStatistics = 0;

	if(cutoff!=0){
		this->cutoff = (int)((cutoff/spacing) + 0.5);
	}else{
		int max_aux = max(width, height);
		int max_final = max(max_aux, depth);
		this->cutoff = (max_final)/2;
	}

	tunnel = false;

	initializeGrid();
}

Grid::Grid(float spacingGrid, float x_min, float y_min, float z_min, int width, int height, int depth) {
	spacing = spacingGrid;

	originX = x_min;
	originY = y_min;
	originZ = z_min;

	this->width = width;
	this->height = height;
	this->depth = depth;

	tStartStatistics = 0;
	initializeGrid();
}

void Grid::initializeGrid(){
    grid = new int**[width];
	for (int x = 0; x < width; ++x) {
		grid[x] = new int*[height];

		for (int y = 0; y < height; ++y)
			grid[x][y] = new int[depth];
	 }

	for (int x = 0; x < width; ++x) {
		for (int y = 0; y < height; ++y) {
			for (int z = 0; z < depth; ++z) {
				grid[x][y][z] = 0;
			}
		}
	}
}

/*TODO: Does the C-alpha thing make sense anymore? */
void Grid::setCAplhaLimits(float x_max, float x_min, float y_max, float y_min, float z_max, float z_min){
	ca_limits[0] = calculatePositionInGrid(x_min + 1*spacing, originX);
	ca_limits[1] = calculatePositionInGrid(x_max - 1*spacing, originX);
	ca_limits[2] = calculatePositionInGrid(y_min + 1*spacing, originY);
	ca_limits[3] = calculatePositionInGrid(y_max - 1*spacing, originY);
	ca_limits[4] = calculatePositionInGrid(z_min + 1*spacing, originZ);
	ca_limits[5] = calculatePositionInGrid(z_max - 1*spacing, originZ);
}


float Grid::calculateDistance(float a, float b){
	float dif = pow(b-a,2);
	return sqrt(dif);
}

float Grid::calculateNumberCells(float a, float b){
	float distance = calculateDistance(a,b);
	int points = fabs(distance/spacing);
	if (fmod(distance,spacing)!= 0) points ++;
	return points;
}

int Grid::calculatePositionInGrid(float coordinate, float origin){
	float dist = calculateDistance(origin,coordinate)/spacing;
	int position = fabs(dist);
	if(dist-position>0.5) position ++;
	return position;
}

Coordinates Grid::calculateCoordinateInGrid(int* position){
	float x, y, z;
	x = (position[0]*spacing) + originX ;
	y = (position[1]*spacing) + originY;
	z = (position[2]*spacing) + originZ;

	Coordinates c(x,y,z);
	c.setGridCoordinate(position);

	return c;
}

int* Grid::getDistanceCells(float radius, int center, int limit){
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

vector<vector<Coordinates> > Grid::getCavitySeedPoint(Coordinates seed, Grid statistics, int &max_cavity_index, int min, bool rectify){
	vector<vector<Coordinates> > cavities;
	vector<Coordinates> cavity;

	int x = calculatePositionInGrid(seed.getX(), originX);
	int y = calculatePositionInGrid(seed.getY(), originY);
	int z = calculatePositionInGrid(seed.getZ(), originZ);

	if((x>0 && x<width) && (y>0 && y<height) && (z>0 && z<depth)){
		if(grid[x][y][z]==0){
			cavity = getCavity(x,y,z, rectify, statistics, -1);
			max_cavity_index=-1;
			if(cavity.size()>1){
				cavities.push_back(cavity);
			}else{
				grid[x][y][z]=0;
			}
		}
		if(cavities.size()<1){
			int max_cavity_position=0;
			vector<vector<Coordinates> > aux = getCavitiesSeedPoint(seed, statistics, max_cavity_index, max_cavity_position, min, rectify);
			if(aux.size()>0){
				cavities.push_back(aux[max_cavity_position]);
			}
		}
	}else{
		cout<<"Seed coordinate is out of the grid "<<x<<" "<<y<<" "<<z<<" "<<width<<" "<<height<<" "<<depth<<"\n";
	}

	return cavities;
}

vector<vector<Coordinates> > Grid::getCavitiesSeedPoint(Coordinates seed, Grid statistics, int &max_cavity_index, int &max_cavity_position, int min, bool rectify){
	int pos = 5/spacing;

	int x = calculatePositionInGrid(seed.getX(), originX);
	int y = calculatePositionInGrid(seed.getY(), originY);
	int z = calculatePositionInGrid(seed.getZ(), originZ);

	int start_x = x - pos;
	int end_x = x + pos;

	if(start_x<0) start_x=0;
	if(end_x>=width) end_x = width-1;

	int start_y = y - pos;
	int end_y = y + pos;

	if(start_y<0) start_y=0;
	if(end_y>=height) end_y = height-1;

	int start_z = z - pos;
	int end_z = z + pos;

	if(start_z<0) start_z=0;
	if(end_z>=depth) end_z = depth-1;

	return getCavities(statistics, max_cavity_index, max_cavity_position, rectify, start_x, end_x, start_y, end_y, start_z, end_z, min);
}

vector<vector<Coordinates> > Grid::getCavitiesNoSeedPoint(Grid statistics, int &max_cavity_index, int &max_cavity_position, int min, bool rectify){
	return getCavities(statistics, max_cavity_index, max_cavity_position, rectify, ca_limits[0], ca_limits[1], ca_limits[2], ca_limits[3], ca_limits[4], ca_limits[5], min);
}

vector<vector<Coordinates> > Grid::getCavities(Grid statistics, int &max_cavity_index, int &max_cavity_position, bool rectify, int x_min, int x_max, int y_min, int y_max, int z_min, int z_max, int min){
	vector<Coordinates> cavity;
	vector<vector<Coordinates> > cavities;

	int cont =1;
	int size, max = 0;

	for(int i=x_min; i<x_max;i++){
		for(int j=y_min; j<y_max;j++){
			for(int k=z_min; k<z_max;k++){
				if(grid[i][j][k]==0){
					cavity = getCavity(i,j,k, false, statistics, -(cont));
					size = cavity.size();
					if(size>1 && size>=min){
						cavities.push_back(cavity);
						if(size>max){
							max = size;
							max_cavity_index = -(cont);
							max_cavity_position = cavities.size()-1;
						}
					}
					cont ++;
				}
			}
		}
	}

	return cavities;
}

/*vector<vector<Coordinates> > Grid::getCavitiesSeedPoint(Coordinates seed, Grid statistics, int &max_cavity_index, int &max_cavity_position, bool rectify){
	vector<Coordinates> cavity;
	vector<vector<Coordinates> > cavities;

	int x = calculatePositionInGrid(seed.getX(), originX);
	int y = calculatePositionInGrid(seed.getY(), originY);
	int z = calculatePositionInGrid(seed.getZ(), originZ);

	int pos = 5/spacing;

	int start_x = x - pos;
	int end_x = x + pos;

	if(start_x<0) start_x=0;
	if(end_x>=width) end_x = width-1;

	int start_y = y - pos;
	int end_y = y + pos;

	if(start_y<0) start_y=0;
	if(end_y>=height) end_y = height-1;

	int start_z = z - pos;
	int end_z = z + pos;

	if(start_z<0) start_z=0;
	if(end_z>=depth) end_z = depth-1;

	return getCavities(statistics, max_cavity_index, max_cavity_position, rectify, start_x, end_x, start_y, end_y, start_z, end_z);

	return cavities;
}*/

vector<Coordinates> Grid::getCavity(int x, int y , int z, bool rectify, Grid statistics, int index){
	vector<Coordinates> cavity;

	vector<int*> candidates;
	vector<int*> neighbours;
	vector<int*> aux;

	int* candidate = new int[3];
	candidate[0] = x; candidate[1] = y; candidate[2] = z;
	candidates.push_back(candidate);
	grid[x][y][z]=index;

	while(!candidates.empty()){
		while(!neighbours.empty()) {delete [] neighbours.back(); neighbours.pop_back();}
		for(unsigned int i=0;i<candidates.size();i++){
			aux.clear();
			candidate =  candidates[i];
			aux = getNeighbours(candidate,rectify, statistics, index);
			neighbours.insert(neighbours.end(), aux.begin(), aux.end());
			cavity.push_back(calculateCoordinateInGrid(candidate));

			Coordinates coord = calculateCoordinateInGrid(candidate);
			int x = calculatePositionInGrid(coord.getX(), originX);
			int y = calculatePositionInGrid(coord.getY(), originY);
			int z = calculatePositionInGrid(coord.getZ(), originZ);
		}
		candidates.swap(neighbours);
	}

	delete [] candidate;

	return cavity;
}

/*vector<int*> Grid::getCandidatesToSeedPoint(int x, int y , int z, float cutoff){
	Coordinates c;
	vector<int*> candidates;

	for(int i=ca_limits[0]; i<ca_limits[1];i++){
		for(int j=ca_limits[2]; j<ca_limits[3];j++){
			for(int k=ca_limits[4]; k<ca_limits[5];k++){
				if(grid[i][j][k]==0){
					int *aux = new int[3];
					aux[0]=i; aux[1]=j; aux[2]=k;
					candidates.push_back(aux);
				}
			}
		}
	}

	if(candidates.size()==0){
		for(int i=0; i<width-1;i++){
			for(int j=0; j<height-1;j++){
				for(int k=0; k<depth-1;k++){
					if(grid[i][j][k]==0){
						int *aux = new int[3];
						aux[0]=i; aux[1]=j; aux[2]=k;
						candidates.push_back(aux);
					}
				}
			}
		}
	}

	return candidates;
}*/

/*bool Grid::isInsideCavity(int x, int y, int z){
	bool result = false;

	int max_aux = max((ca_limits[1]-ca_limits[0]), (ca_limits[3]-ca_limits[2]));
	int max_final = max(max_aux, (ca_limits[5]-ca_limits[4]));

	int pos = (max_final)/2;

	int protein_neighbours = surroundedBy(x,y,z, 1, pos);

	if(protein_neighbours>=5) {
		result = true;
	}else{
		int cavity_neighbours = surroundedBy(x,y,z, 2, pos);
		if(cavity_neighbours>=3){
			result = true;
		}
	}

	return result;
}*/

//delete?
/*int Grid::surroundedBy(int x, int y, int z, int element, int pos){

	return surroundedByX(x,y,z,element,pos) + surroundedByY(x,y,z,element,pos) + surroundedByZ(x,y,z,element,pos);
}*/

bool Grid::isInsideCavity(int x, int y, int z, int pos, int index, bool tunnel){
	bool result = false;

	int x_pos = surroundedByXPositive(x,y,z,1,pos);
	int x_neg = surroundedByXNegative(x,y,z,1,pos);
	int y_pos = surroundedByYPositive(x,y,z,1,pos);
	int y_neg = surroundedByYNegative(x,y,z,1,pos);
	int z_pos = surroundedByZPositive(x,y,z,1,pos);
	int z_neg = surroundedByZNegative(x,y,z,1,pos);

	int dimensions = 0;

	if(!tunnel){
		dimensions = this->getDimensions(x_pos, x_neg, y_pos, y_neg, z_pos, z_neg);
	}else{
		dimensions = this->getDimensionsTunnel(x_pos, x_neg, y_pos, y_neg, z_pos, z_neg);
	}

	if(dimensions>=this->dimensionsSearch){
		result = true;
	}else{
		result = doRefill(x,y,z,pos, index);
	}

	return result;
}

int Grid::getDimensions(int x_pos, int x_neg, int y_pos, int y_neg, int z_pos, int z_neg){
	int dimensions = 0;

	if (x_pos>0) dimensions ++;
	if (x_neg>0) dimensions ++;
	if (y_pos>0) dimensions ++;
	if (y_neg>0) dimensions ++;
	if (z_pos>0) dimensions ++;
	if (z_neg>0) dimensions ++;

	return dimensions;
}

int Grid::getDimensionsTunnel(int x_pos, int x_neg, int y_pos, int y_neg, int z_pos, int z_neg){
	int dimensions = 0;

	bool odd = false;

	if (dimensionsSearch%2==1){
		odd=true;
	}

	if(x_pos>0 && x_neg>0){
		dimensions = dimensions + 2;
	}else{
		if(odd==true && (x_pos>0 || x_neg>0)){
			dimensions = dimensions + 1;
			odd = false;
		}
	}

	if(y_pos>0 && y_neg>0){
		dimensions = dimensions + 2;
	}else{
		if(odd==true && (y_pos>0 || y_neg>0)){
			dimensions = dimensions + 1;
			odd = false;
		}
	}

	if(z_pos>0 && z_neg>0){
		dimensions = dimensions + 2;
	}else{
		if(odd==true && (z_pos>0 || z_neg>0)){
			dimensions = dimensions + 1;
			odd = false;
		}
	}

	return dimensions;
}

bool Grid::doRefill(int x, int y, int z, int pos, int index){

	bool result = false;

	int x_pos = surroundedByXPositive(x,y,z,index,pos);
	int x_neg = surroundedByXNegative(x,y,z,index,pos);
	int y_pos = surroundedByYPositive(x,y,z,index,pos);
	int y_neg = surroundedByYNegative(x,y,z,index,pos);
	int z_pos = surroundedByZPositive(x,y,z,index,pos);
	int z_neg = surroundedByZNegative(x,y,z,index,pos);

	int count = 0;

	if (x_pos>0 && x_pos<=3) count++;
	if (x_neg>0 && x_neg<=3) count++;
	if (y_pos>0 && y_pos<=3) count++;
	if (y_neg>0 && y_neg<=3) count++;
	if (z_pos>0 && z_pos<=3) count++;
	if (z_neg>0 && z_neg<=3) count++;

	if(count>=4) result = true;

	return result;
}

/*bool Grid::isInsideCavityNew(int x, int y, int z, int pos, int index){

	bool result = false;

	int boxed_result = 0;
	bool odd = false;

	bool sx,sy,sz = false;

	if (dimensionsSearch%2==1){
		odd=true;
	}

	int x_pos = surroundedByXPositive(x,y,z,1,pos);
	int x_neg = surroundedByXNegative(x,y,z,1,pos);
	int y_pos = surroundedByYPositive(x,y,z,1,pos);
	int y_neg = surroundedByYNegative(x,y,z,1,pos);
	int z_pos = surroundedByZPositive(x,y,z,1,pos);
	int z_neg = surroundedByZNegative(x,y,z,1,pos);

	if(x_pos>0 && x_neg>0){
		boxed_result = boxed_result + 2;
		sx=true;
	}else{
		if(odd==true && (x_pos>0 || x_neg>0)){
			boxed_result = boxed_result + 1;
			odd = false;
		}
	}

	if(y_pos>0 && y_neg>0){
		boxed_result = boxed_result + 2;
		sy=true;
	}else{
		if(odd==true && (y_pos>0 || y_neg>0)){
			boxed_result = boxed_result + 1;
			odd = false;
		}
	}

	if(z_pos>0 && z_neg>0){
		boxed_result = boxed_result + 2;
		sz=true;
	}else{
		if(odd==true && (z_pos>0 || z_neg>0)){
			boxed_result = boxed_result + 1;
			odd = false;
		}
	}

	if(boxed_result>=dimensionsSearch){
		result = true;
	}else{
		int boxed_result_2 = 0;

		int x_pos_2 = surroundedByXPositive(x,y,z,index,pos);
		int x_neg_2 = surroundedByXNegative(x,y,z,index,pos);
		int y_pos_2 = surroundedByYPositive(x,y,z,index,pos);
		int y_neg_2 = surroundedByYNegative(x,y,z,index,pos);
		int z_pos_2 = surroundedByZPositive(x,y,z,index,pos);
		int z_neg_2 = surroundedByZNegative(x,y,z,index,pos);

		if (dimensionsSearch%2==1){
			odd=true;
		}

		if(sx==false){
			if((x_pos_2>0 && x_neg_2>0) || (x_pos>0 && x_neg_2>0) || (x_pos_2>0 && x_neg>0)){
				boxed_result_2 = boxed_result_2 + 2;
			}else{
				if(odd==true && (x_pos_2>0 || x_neg_2>0)){
					boxed_result_2 = boxed_result_2 + 1;
					odd = false;
				}
			}
		}

		if(sy==false){
			if((y_pos_2>0 && y_neg_2>0) || (y_pos>0 && y_neg_2>0) || (y_pos_2>0 && y_neg>0)){
				boxed_result_2 = boxed_result_2 + 2;
			}else{
				if(odd==true && (y_pos_2>0 || y_neg_2>0)){
					boxed_result_2 = boxed_result_2 + 1;
					odd = false;
				}
			}
		}

		if(sz==false){
			if((z_pos_2>0 && z_neg_2>0) || (z_pos>0 && z_neg_2>0) || (z_pos_2>0 && z_neg>0)){
				boxed_result_2 = boxed_result_2 + 2;
			}else{
				if(odd==true && (z_pos_2>0 || z_neg_2>0)){
					boxed_result_2 = boxed_result_2 + 1;
					odd = false;
				}
			}
		}


		if(boxed_result_2+boxed_result>=dimensionsSearch){
			result = true;
		}else{
			//refill
			int count = 0;
			if (x_pos_2>0 && x_pos_2<=3) count++;
			if (x_neg_2>0 && x_neg_2<=3) count++;
			if (y_pos_2>0 && y_pos_2<=3) count++;
			if (y_neg_2>0 && y_neg_2<=3) count++;
			if (z_pos_2>0 && z_pos_2<=3) count++;
			if (z_neg_2>0 && z_neg_2<=3) count++;
			if(count>=4) result = true;
		}
	}

	return result;
}*/

int Grid::surroundedByX(int x, int y, int z, int element, int pos){

	int start_x = x - pos;
	int end_x = x + pos;

	int count = 0;

	if(start_x<0) start_x=0;
	if(end_x>=width) end_x = width-1;

	for(int i=x; i>start_x;i--){
		if(grid[i][y][z]==element){
			count ++;
			break;
		}
	}

	for(int i=x; i<end_x;i++){
		if(grid[i][y][z]==element){
			count ++;
			break;
		}
	}

	return count;
}

int Grid::surroundedByY(int x, int y, int z, int element, int pos){

	int start_y = y - pos;
	int end_y = y + pos;

	int count = 0;

	if(start_y<0) start_y=0;
	if(end_y>=height) end_y = height-1;

	for(int j=y; j>start_y; j--){
		if(grid[x][j][z]==element){
			count ++;
			break;
		}
	}
	for(int j=y; j<end_y; j++){
		if(grid[x][j][z]==element){
			count ++;
			break;
		}
	}
	return count;
}

int Grid::surroundedByZ(int x, int y, int z, int element, int pos){

	int start_z = z - pos;
	int end_z = z + pos;

	if(start_z<0) start_z=0;
	if(end_z>=depth) end_z = depth-1;

	int count = 0;

	for(int k=z; k>start_z; k--){
		if(grid[x][y][k]==element){
			count ++;
			break;
		}
	}
	for(int k=z; k<end_z; k++){
		if(grid[x][y][k]==element){
			count ++;
			break;
		}
	}

	return count;
}

int Grid::surroundedByXPositive(int x, int y, int z, int element, int pos){

	int end = x + pos;

	int count = 0;

	if(end>=width) end = width-1;

	for(int i=x+1; i<=end;i++){
		//count ++;
		if(element>0 && grid[i][y][z]>0){
			//protein
			count = 1;
			break;
		}else{
			if(grid[i][y][z]==element){
				//cavity
				count = 1;
				break;
			}
		}
	}

	return count;
}

int Grid::surroundedByXNegative(int x, int y, int z, int element, int pos){

	int start = x - pos;

	int count = 0;

	if(start<0) start=0;

	for(int i=x-1; i>=start;i--){
		//count ++;
		if(element>0 && grid[i][y][z]>0){
			count = 1;
			break;
		}else{
			if(grid[i][y][z]==element){
				count = 1;
				break;
			}
		}
	}

	return count;
}

int Grid::surroundedByYPositive(int x, int y, int z, int element, int pos){

	int end = y + pos;

	int count = 0;

	if(end>=height) end = height-1;

	for(int i=y+1; i<=end;i++){
		//count ++;
		if(element>0 && grid[x][i][z]>0){
			count = 1;
			break;
		}else{
			if(grid[x][i][z]==element){
				count = 1;
				break;
			}
		}
	}

	return count;
}

int Grid::surroundedByYNegative(int x, int y, int z, int element, int pos){

	int start = y - pos;

	int count = 0;

	if(start<0) start=0;

	for(int i=y-1; i>=start;i--){
		//count ++;
		if(element>0 && grid[x][i][z]>0){
			count = 1;
			break;
		}else{
			if(grid[x][i][z]==element){
				count = 1;
				break;
			}
		}
	}

	return count;
}

int Grid::surroundedByZPositive(int x, int y, int z, int element, int pos){

	int end = z + pos;

	int count = 0;

	if(end>=depth) end = depth-1;

	for(int i=z+1; i<=end;i++){
		//count ++;
		if(element>0 && grid[x][y][i]>0){
			count = 1;
			break;
		}else{
			if(grid[x][y][i]==element){
				count = 1;
				break;
			}
		}
	}

	return count;
}

int Grid::surroundedByZNegative(int x, int y, int z, int element, int pos){

	int start = z - pos;

	int count = 0;

	if(start<0) start=0;

	for(int i=z-1; i>=start;i--){
		//count ++;
		if(element>0 && grid[x][y][i]>0){
			count = 1;
			break;
		}else{
			if(grid[x][y][i]==element){
				count = 1;
				break;
			}
		}
	}

	return count;
}

void Grid::getNeighbourCandidate(int x, int y, int z, vector<int*>& neighbours, bool rectify, Grid statistics, int index){

	if(grid[x][y][z]==0){

		if(rectify){
			if((x>0) && (x<statistics.getWidth()) && (y>0) && (y<statistics.getHeight()) && (z>0) && (z<statistics.getDepth())){
				float percent = (float)statistics.getGrid()[x][y][z]/50.0;
				if(percent>0.1){
					if(isInsideCavity(x,y,z,cutoff,index,this->tunnel)){
						int* n1 = new int[3];
						n1[0] = x; n1[1] = y; n1[2] = z;
						neighbours.push_back(n1);
						//statistics.getGrid()[x][y][z] = statistics.getGrid()[x][y][z]+1;
						grid[x][y][z]=index;
					}
				}
			}
		}else{
			if(isInsideCavity(x,y,z,cutoff,index, this->tunnel)){
				int* n1 = new int[3];
				n1[0] = x; n1[1] = y; n1[2] = z;
				neighbours.push_back(n1);
				grid[x][y][z]=index;
			}
		}
	}
}

vector<int*> Grid::getNeighbours(int* seed, bool rectify, Grid statistics, int index){
	vector<int*> neighbours;

	int start_x = seed[0]-1;
	if(start_x<ca_limits[0])  start_x = ca_limits[0];

	int end_x = seed[0]+1;
	if(end_x>ca_limits[1])  end_x = ca_limits[1];

	int start_y = seed[1]-1;
	if(start_y<ca_limits[2])  start_y = ca_limits[2];

	int end_y = seed[1]+1;
	if(end_y>ca_limits[3])  end_y = ca_limits[3];

	int start_z = seed[2]-1;
	if(start_z<ca_limits[4])  start_z = ca_limits[4];

	int end_z = seed[2]+1;
	if(end_z>ca_limits[5])  end_z = ca_limits[5];

	for(int x = start_x; x<=end_x; x++){
		for(int y = start_y; y<=end_y; y++){
			for(int z = start_z; z<=end_z; z++){
				getNeighbourCandidate(x,y,z,neighbours, rectify, statistics, index);
			}
		}
	}
	return neighbours;
}

void Grid::setAxisDependentInfo(int axis, int &ilimit, int &jlimit, int &klimit, float &origin, int &a, int &b){
	if(axis==0){
		//X axis
		ilimit = height;
		jlimit = depth;
		klimit = width;
		origin = originX;
		a = 1;
		b = 2;
	}else{
		if(axis==1){
			//Y axis
			ilimit = width;
			jlimit = depth;
			klimit = height;
			origin = originY;
			a = 0;
			b = 2;
		}else{
			//Z axis
			ilimit = width;
			jlimit = height;
			klimit = depth;
			origin = originZ;
			a = 0;
			b = 1;
		}
	}
}

map<int,float> Grid::calculateBottleneckArea(int index,int axis, vector<Coordinates>& tunnel, float &bottleneck, float axis_value ,bool calculate_sector, bool radius){

	int i,j,k,a,b;
	int ilimit, jlimit, klimit;
	float origin;
	unsigned int size, max = 0;

	map<int, float> sector;

	vector<Coordinates> line;
	//vector<vector<Coordinates> > lines;
	vector<Coordinates> max_line;
	Coordinates c;
	int plane = 0;
	float round_value = 0;

	setAxisDependentInfo(axis, ilimit, jlimit, klimit, origin, a, b);

	for(i=0; i<ilimit; i++){
		for (j=0; j<jlimit; j++){
			line.clear();
			for (k=0;k<klimit; k++){
				if(grid[i][j][k]>0){
					break;
				}else{
					if(grid[i][j][k]==index){
						int *coor = new int[3];
						coor[0]=i; coor[1]=j; coor[2]=k;
						c = calculateCoordinateInGrid(coor);
						line.push_back(c);
						delete [] coor;
					}
				}
			}

			size = line.size();

			if(size>max){
				max_line.clear();
				max_line = line;
				max = size;
			}

			//If the cavity line gets to the opposite plane it is a tunnel
			if((k==klimit)  && size>klimit/2){
				//if(line.size()>max) max = (int)line.size();
				tunnel.insert(tunnel.end(), line.begin(), line.end());
				plane ++;
			}
		}
	}

	bottleneck=plane*(spacing*spacing);

	if(!tunnel.size()>0){
		tunnel = max_line;
	}

	if(tunnel.size()>0){

		sector = getSectorTunnel(index,axis,a,b,tunnel,bottleneck, radius, axis_value);

		/*if(axis_value!=NULL){

			int inspect_grid = this->calculatePositionInGrid(axis_value, origin);

			if(axis_value<0){
				round_value = (int)(inspect_grid-0.5);
			}else{
				round_value = (int)(inspect_grid+0.5);
			}

			int cont = 0;
			while(sector[round_value+cont]==0.0){
				cont++;
			}
			bottleneck = sector[round_value+cont];

			//for(unsigned t=0; t<tunnel.size(); t++){
			//	if(tunnel[t].getGridCoordinate()[axis]==inspect_grid) {
			//		bottleneck = getAreaAxis(tunnel[t].getGridCoordinate(), index, axis, a, b);
			//		break;
			//	}
			//}
		}*/
	}

	return sector;
}

map<int, float> Grid::getSectorTunnel(int index, int axis, int a, int b, vector<Coordinates>& tunnel, float &bottleneck, bool radius, float inspect_axis){

	vector<Coordinates> tunnel_extension;
	map<pair<int,int>,int > looked_up;
	map<int, float> sector;
	vector<int*> candidates;
	Coordinates coordinates;
	float refined_area;
	float botleneck_candidate = 0.0;

	int grid_value, round_value;
	int *member;
	float axis_value;
	float bottleneck_round=0.0;

	for(unsigned l=0; l<tunnel.size(); l++){

		axis_value = tunnel[l].getCoordinate(axis);
		//I round it to int, it makes it easier to compare keys
		if(axis_value<0){
			round_value = (int)(axis_value-0.5);
		}else{
			round_value = (int)(axis_value+0.5);
		}


		if(sector.find(round_value)==sector.end()){

			refined_area = 0.0;

			while(!candidates.empty()){ delete candidates.back(); candidates.pop_back();}
			if(looked_up.size()>0) looked_up.clear();
			if(tunnel_extension.size()>0) tunnel_extension.clear();

			candidates.push_back(tunnel[l].getGridCoordinate());

			while(candidates.size()>0){

				member = candidates.back();
				candidates.pop_back();

				looked_up[make_pair(member[a],member[b])]=0;

				for(int i=-1; i<=1; i++){
					for(int j=-1; j<=1; j++){
						if(looked_up.find(make_pair(member[a]+i, member[b]+j))==looked_up.end()){
							if(axis==0){
								grid_value = grid[member[0]][member[1]+i][member[2]+j];
							}else{
								if(axis==1){
									grid_value = grid[member[0]+i][member[1]][member[2]+j];
								}else{
									grid_value = grid[member[0]+i][member[1]+j][member[2]];
								}
							}

							if(grid_value==index){
								int *c = new int[3];
								c[a]=member[a]+i; c[b]=member[b]+j; c[axis]=member[axis];
								candidates.push_back(c);
								looked_up[make_pair(c[a],c[b])]=0;

							}
						}
					}
				}
				coordinates = calculateCoordinateInGrid(member);
				tunnel_extension.push_back(coordinates);
				refined_area = refined_area + refineArea(coordinates.getX(), coordinates.getY(), coordinates.getZ(), axis);
			}

			if(!radius){
				sector[round_value] = refined_area;
			}else{
				sector[round_value] = sqrt(refined_area/Pi);
			}

			if(inspect_axis!=NULL){
				if(axis_value<0){
					bottleneck_round = (int)(inspect_axis-0.5);
					if((botleneck_candidate==0.0 && (round_value<=(bottleneck_round+1) && round_value>=(bottleneck_round-1))) || ((botleneck_candidate>sector[round_value]) && (round_value<=(bottleneck_round+1) && round_value>=(bottleneck_round-1)))){
						botleneck_candidate = sector[round_value];
					}
				}else{
					bottleneck_round = (int)(inspect_axis+0.5);
					if((botleneck_candidate==0.0 && (round_value>=(bottleneck_round+1) && round_value<=(bottleneck_round-1))) || ((botleneck_candidate>sector[round_value]) && (round_value>=(bottleneck_round+1) && round_value<=(bottleneck_round-1)))){
						botleneck_candidate = sector[round_value];
					}
				}
			}else{
				if(botleneck_candidate==0.0|| botleneck_candidate>sector[round_value]){
					botleneck_candidate = sector[round_value];
				}
			}
			tunnel.insert(tunnel.end(), tunnel_extension.begin(), tunnel_extension.end());
		}
	}

	bottleneck = botleneck_candidate;

	return sector;
}

/*map<int, float> Grid::getSectorTunnel(int index, int axis, int a, int b, vector<Coordinates>& tunnel, int plane, bool radius){

	vector<Coordinates> tunnel_extension;
	map<pair<int,int>,int > looked_up;
	map<int, float> sector;
	vector<int*> candidates;
	Coordinates coordinates;
	float radius_value, refined_area;

	int sector_acc, grid_value, round_value;
	int *member;
	float axis_value;
	//set<int> diam1, diam2;
	//int protrusion, erase,;

	for(unsigned l=0; l<tunnel.size(); l++){

		axis_value = tunnel[l].getCoordinate(axis);
		//I round it to int, it makes it easier to compare keys
		if(axis_value<0){
			round_value = (int)(axis_value-0.5);
		}else{
			round_value = (int)(axis_value+0.5);
		}

		//diam1.clear();
		//diam2.clear();

		if(sector.find(round_value)==sector.end()){

			refined_area = 0.0;

			while(!candidates.empty()){ delete candidates.back(); candidates.pop_back();}
			if(looked_up.size()>0) looked_up.clear();
			if(tunnel_extension.size()>0) tunnel_extension.clear();

			sector_acc = 0;

			candidates.push_back(tunnel[l].getGridCoordinate());

			while(candidates.size()>0){

				member = candidates.back();
				candidates.pop_back();

				looked_up[make_pair(member[a],member[b])]=0;
				//diam1.insert(member[a]);
				//diam2.insert(member[b]);
				//protrusion =0;
				//erase = candidates.size();


				//I get a list of candidates that I don't consider a protrusion (to extend the bottleneck tunnel that already exists)
				for(int i=-1; i<=1; i++){
					for(int j=-1; j<=1; j++){
						if(looked_up.find(make_pair(member[a]+i, member[b]+j))==looked_up.end()){
							if(axis==0){
								grid_value = grid[member[0]][member[1]+i][member[2]+j];
							}else{
								if(axis==1){
									grid_value = grid[member[0]+i][member[1]][member[2]+j];
								}else{
									grid_value = grid[member[0]+i][member[1]+j][member[2]];
								}
							}

							if(grid_value==index){
								int *c = new int[3];
								c[a]=member[a]+i; c[b]=member[b]+j; c[axis]=member[axis];
								candidates.push_back(c);
								looked_up[make_pair(c[a],c[b])]=0;
								//diam1.insert(member[a]+i);
								//diam2.insert(member[b]+j);
								//protrusion++;
							}
						}//else{
							//protrusion++;
						//}
					}
				}

				sector_acc ++;
				coordinates = calculateCoordinateInGrid(member);
				tunnel_extension.push_back(coordinates);
				refined_area = refined_area + refineArea(coordinates.getX(), coordinates.getY(), coordinates.getZ(), axis);

				//Probaly to delete
				//if(protrusion>8){ //square-like section only (for radius calculation)
				//	sector_acc ++;
				//	coordinates = calculateCoordinateInGrid(member);
				//	tunnel_extension.push_back(coordinates);
				//}else{
				//	candidates.erase(candidates.begin()+erase,candidates.end());
				//}
			}

			if(!radius){
				//if(sector_acc<plane) sector_acc =plane;
				//sector_acc*(spacing*spacing);
				sector[round_value] = refined_area ;
			}else{
				//radius_value = min(diam1.size(),diam2.size())/2.0;
				//sector[round_value] = radius_value;
				sector[round_value] = sqrt(refined_area/3.1416);
			}

			tunnel.insert(tunnel.end(), tunnel_extension.begin(), tunnel_extension.end());
		}
	}

	return sector;
}*/

/*float Grid::getAreaAxis(vector<Coordinates> cavity, int axis, float value){

	float area = 0;

	if(!cavity.empty() && cavity.size()>0){

		int ilimit, jlimit, klimit, a, b;
		map<pair<int,int>,int > looked_up;
		vector<int*> candidates;
		int* member;
		int grid_value = 0;

		int index = grid[cavity[0].getGridCoordinate()[0]][cavity[0].getGridCoordinate()[1]][cavity[0].getGridCoordinate()[2]];

		setAxisDependentInfo(axis, ilimit, jlimit, klimit, a, b);

		int inspect_grid = this->calculatePositionInGrid(value, originZ);
		for(unsigned c=0; c<cavity.size();c++){
			if(cavity[c].getGridCoordinate()[axis]==inspect_grid) {
				candidates.push_back(cavity[c].getGridCoordinate());
				break;
			}
		}

		while(candidates.size()>0){

			member = candidates.back();
			candidates.pop_back();
			area ++;

			looked_up[make_pair(member[a],member[b])]=0;

			for(int i=-1; i<=1; i++){
				for(int j=-1; j<=1; j++){
					if(looked_up.find(make_pair(member[a]+i, member[b]+j))==looked_up.end()){
						if(axis==0){
							grid_value = grid[member[0]][member[1]+i][member[2]+j];
						}else{
							if(axis==1){
								grid_value = grid[member[0]+i][member[1]][member[2]+j];
							}else{
								grid_value = grid[member[0]+i][member[1]+j][member[2]];
							}
						}

						if(grid_value==index){
							int *c = new int[3];
							c[a]=member[a]+i; c[b]=member[b]+j; c[axis]=member[axis];
							candidates.push_back(c);
							looked_up[make_pair(c[a],c[b])]=0;
						}
					}
				}
			}
		}
	}

	return area*(spacing*spacing);

}*/


/*float Grid::getAreaAxis(int* seed, int index, int axis, int a, int b){

	float bottleneck = 0;
	int grid_value = 0;

	map<pair<int,int>,int > looked_up;
	vector<int*> candidates;
	int* member;

	candidates.push_back(seed);

	while(candidates.size()>0){

		member = candidates.back();
		candidates.pop_back();
		bottleneck ++;

		looked_up[make_pair(member[a],member[b])]=0;

		for(int i=-1; i<=1; i++){
			for(int j=-1; j<=1; j++){
				if(looked_up.find(make_pair(member[a]+i, member[b]+j))==looked_up.end()){
					if(axis==0){
						grid_value = grid[member[0]][member[1]+i][member[2]+j];
					}else{
						if(axis==1){
							grid_value = grid[member[0]+i][member[1]][member[2]+j];
						}else{
							grid_value = grid[member[0]+i][member[1]+j][member[2]];
						}
					}

					if(grid_value==index){
						int *c = new int[3];
						c[a]=member[a]+i; c[b]=member[b]+j; c[axis]=member[axis];
						candidates.push_back(c);
						looked_up[make_pair(c[a],c[b])]=0;
					}
				}
			}
		}
	}

	return bottleneck*(spacing*spacing);
}*/

float Grid::refineArea(float xf, float yf, float zf, int axis){

	int x = calculatePositionInGrid(xf, originX);
	int y = calculatePositionInGrid(yf, originY);
	int z = calculatePositionInGrid(zf, originZ);

	float x_length = spacing;
	float y_length = spacing;
	float z_length = spacing;

	float vdw, container, result = 0.0;

	if(grid[x][y][z]<0){
		if(axis!=0){
			if(x-1>0 && grid[x-1][y][z]>0){
				vdw = atoms[grid[x-1][y][z]-1].getVdwRadius();
				container = xf - (atoms[grid[x-1][y][z]-1].getCoordinates().getX()+vdw);
				x_length = x_length + container;
			}
			if(x+1<width && grid[x+1][y][z]>0){
				vdw = atoms[grid[x+1][y][z]-1].getVdwRadius();
				container = ((atoms[grid[x+1][y][z]-1].getCoordinates().getX()-vdw)-spacing) - xf;
				x_length = x_length + container;
			}
		}
		if(axis!=1){
			if(y-1>0 && grid[x][y-1][z]>0){
				vdw = atoms[grid[x][y-1][z]-1].getVdwRadius();
				container = yf - (atoms[grid[x][y-1][z]-1].getCoordinates().getY()+vdw);
				y_length = y_length + container;
			}
			if(y+1<height && grid[x][y+1][z]>0){
				vdw = atoms[grid[x][y+1][z]-1].getVdwRadius();
				container = ((atoms[grid[x][y+1][z]-1].getCoordinates().getY()-vdw)-spacing) - yf;
				y_length = y_length + container;
			}
		}

		/*if(x_length==spacing && y_length==spacing){
			if(x-1>0 && y-1>0 && grid[x-1][y-1][z]>0){
				vdw = atoms[grid[x-1][y-1][z]-1].getVdwRadius();
				container = xf - (atoms[grid[x-1][y-1][z]-1].getCoordinates().getX()+vdw);
				x_length = x_length + container;
				container = yf - (atoms[grid[x-1][y-1][z]-1].getCoordinates().getY()+vdw);
				y_length = y_length + container;
			}
			if(x-1>0 && y+1<height && grid[x-1][y+1][z]>0){
				vdw = atoms[grid[x-1][y+1][z]-1].getVdwRadius();
				container = xf - (atoms[grid[x-1][y+1][z]-1].getCoordinates().getX()+vdw);
				x_length = x_length + container;
				container = ((atoms[grid[x-1][y+1][z]-1].getCoordinates().getY()-vdw)-spacing) - yf;
				y_length = y_length + container;
			}
			if(x+1<width && y-1>0 && grid[x+1][y-1][z]>0){
				vdw = atoms[grid[x+1][y-1][z]-1].getVdwRadius();
				container = ((atoms[grid[x+1][y-1][z]-1].getCoordinates().getX()-vdw)-spacing) - xf;
				x_length = x_length + container;
				container = yf - (atoms[grid[x+1][y-1][z]-1].getCoordinates().getY()+vdw);
				y_length = y_length + container;
			}
			if(x+1<width && y+1<height && grid[x+1][y+1][z]>0){
				vdw = atoms[grid[x+1][y+1][z]-1].getVdwRadius();
				container = ((atoms[grid[x+1][y+1][z]-1].getCoordinates().getX()-vdw)-spacing) - xf;
				x_length = x_length + container;
				container = ((atoms[grid[x+1][y+1][z]-1].getCoordinates().getY()-vdw)-spacing) - yf;
				y_length = y_length + container;
			}
		}*/

		if(axis!=2){
			if(z-1>0 && grid[x][y][z-1]>0){
				vdw = atoms[grid[x][y][z-1]-1].getVdwRadius();
				container = zf - (atoms[grid[x][y][z-1]-1].getCoordinates().getZ()+vdw);
				z_length = z_length + container;
			}
			if(z+1<depth && grid[x][y][z+1]>0){
				vdw = atoms[grid[x][y][z+1]-1].getVdwRadius();
				container = ((atoms[grid[x][y][z+1]-1].getCoordinates().getZ()-vdw)-spacing) - zf;
				z_length = z_length + container;
			}
		}
	}

	switch(axis){
		case 0:
			result = y_length * z_length;
			break;
		case 1:
			result = x_length * z_length;
			break;
		case 2:
			result = x_length * y_length;
			break;
	}

	return result;
}

/*float Grid::refineVolume(float xf, float yf, float zf){

	int x = calculatePositionInGrid(xf, originX);
	int y = calculatePositionInGrid(yf, originY);
	int z = calculatePositionInGrid(zf, originZ);
	float vdw;
	float result, container;
	float x_length = spacing;
	float y_length = spacing;
	float z_length = spacing;


	if(grid[x][y][z]<0){
		if(x-1>0 && grid[x-1][y][z]>0){
			vdw = atoms[grid[x-1][y][z]-1].getVdwRadius();
			container = xf - (atoms[grid[x-1][y][z]-1].getCoordinates().getX()+vdw);
			x_length = x_length + container;
		}
		if(x+1<width && grid[x+1][y][z]>0){
			vdw = atoms[grid[x+1][y][z]-1].getVdwRadius();
			container = ((atoms[grid[x+1][y][z]-1].getCoordinates().getX()-vdw)-spacing) - xf;
			x_length = x_length + container;
		}
		if(y-1>0 && grid[x][y-1][z]>0){
			vdw = atoms[grid[x][y-1][z]-1].getVdwRadius();
			container = yf - (atoms[grid[x][y-1][z]-1].getCoordinates().getY()+vdw);
			y_length = y_length + container;
		}
		if(y+1<height && grid[x][y+1][z]>0){
			vdw = atoms[grid[x][y+1][z]-1].getVdwRadius();
			container = ((atoms[grid[x][y+1][z]-1].getCoordinates().getY()-vdw)-spacing) - yf;
			y_length = y_length + container;
		}
		if(z-1>0 && grid[x][y][z-1]>0){
			vdw = atoms[grid[x][y][z-1]-1].getVdwRadius();
			container = zf - (atoms[grid[x][y][z-1]-1].getCoordinates().getZ()+vdw);
			z_length = z_length + container;
		}
		if(z+1<depth && grid[x][y][z+1]>0){
			vdw = atoms[grid[x][y][z+1]-1].getVdwRadius();
			container = ((atoms[grid[x][y][z+1]-1].getCoordinates().getZ()-vdw)-spacing) - zf;
			z_length = z_length + container;
		}
	}

	result = x_length * y_length * z_length;

	return result;
}*/

void Grid::deleteGrid(){
	for (int x = 0; x < width; ++x) {
		for (int y = 0; y < height; ++y){
			delete[] grid[x][y];
		}

		delete[] grid[x];
	 }

	delete[] grid;
}

int Grid::getWidth(){
	return width;
}

int Grid::getHeight(){
	return height;
}

int Grid::getDepth(){
	return depth;
}

float Grid::getOriginX(){
	return originX;
}

float Grid::getOriginY(){
	return originY;
}

float Grid::getOriginZ(){
	return originZ;
}

float Grid::getSpacing(){
	return spacing;
}

int*** Grid::getGrid(){
	return grid;
}

void Grid::setTStartStatisitics(float tstat){
	tStartStatistics = tstat;
}

float Grid::getTStartStatisitics(){
	return tStartStatistics;
}

void Grid::setDimensionsSearch(int dim){
	dimensionsSearch = dim;
}

int Grid::getDimensionsSearch(){
	return dimensionsSearch;
}

vector<Atom> Grid::getAtoms(){
    return atoms;
}

void Grid::setAtoms(vector<Atom> atoms){
    this->atoms = atoms;
}

void Grid::setTunnel(bool tunnel){
	this->tunnel = tunnel;
}

bool Grid::getTunnel(){
	return this->tunnel;
}
Grid::~Grid() {
}
