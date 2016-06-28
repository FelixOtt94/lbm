#include <iostream>
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <vector>
#include <math.h>

#include "imageClass/GrayScaleImage.h"
#include "imageClass/lodepng.h"

#define TIMESTEP 0.001
#define SPACEING 0.001

using namespace std;

static bool scenario1 = false;
static double viscosity_l = 0.0;
static double timeToSimulate = 0.0;
static double acceleration = 0.0;
static int resolution = 0;
static double timestep = 0.0;
static double spaceing = 0.0;
static double omega = 0.0;
static double force_l = 0.0;
static int numCellsX = 0;
static int numCellsY = 0;

// template class for a 2D grid
template<typename T> class grid_lattice {

private:
    size_t lengthInX;
    size_t lengthInY;
    size_t lengthInF;
    std::vector<T> data;

public:
    grid_lattice() {
        lengthInX = 0;
        lengthInY = 0;
        lengthInF = 0;
    }
    grid_lattice(size_t xDim, size_t yDim, size_t f=9) {
        data = std::vector<T>((xDim+1) * (yDim+1) * f, (T) (0.0));
        lengthInX = xDim+1;
        lengthInY = yDim+1;
        lengthInF = f;
    }		// Standart-constructor
    grid_lattice(size_t xDim, size_t yDim, T value, size_t f=9) {
        data = std::vector<T>((xDim+1) * (yDim+1) * f, value);
        lengthInX = xDim+1;
        lengthInY = yDim+1;
        lengthInF = f;
    }	// initalisation-constructor
    virtual ~grid_lattice() {
    }		// Destructor

    size_t lengthX() {
        return lengthInX;
    }	// returns number of elements in x direction
    size_t lengthY() {
        return lengthInX;
    }	// returns number of elements in y direction
    size_t lengthF() {
        return lengthInF;
    }	// returns number of elements in f direction

    T& operator()(size_t i, size_t j, size_t f=0) {
        assert(i < lengthInX);
        assert(j < lengthInY);
        assert(f < lengthInF);
        return data[j * ((int)lengthInX*(int)lengthInF) + i*(int)lengthInF +(int)f];
    }

    void getCopy( grid_lattice &source ){
        for( int i=0; i<(int)lengthInX; i++){
            for( int j=0; j<(int)lengthInY; j++){
                for( int f=0; f<(int)lengthInF; f++){
                    data[j * ((int)lengthInX*(int)lengthInF) + i*(int)lengthInF +f] = source( i,j,f);
                }
            }
        }
    }

};

static grid_lattice<double> grid;
static grid_lattice<double> gridCopy;
static grid_lattice<int> isBoundary;

void streamStep(){
    gridCopy.getCopy( grid );

    int lengthX = (int)grid.lengthX();
    int lengthY = (int)grid.lengthY();

    // inner grid points without cylinder
    for( int i=1; i<lengthX-1; i++ ){
        for( int j=1; j<lengthY-1; j++ ){
            if( !isBoundary(i,j) && i>0 && i<numCellsX-1 ){
                grid(i+1, j,   1) = gridCopy( i,j,1);
                grid(i+1, j+1, 2) = gridCopy( i,j,2);
                grid(i,   j+1, 3) = gridCopy( i,j,3);
                grid(i-1, j+1, 4) = gridCopy( i,j,4);
                grid(i-1, j,   5) = gridCopy( i,j,5);
                grid(i-1, j-1, 6) = gridCopy( i,j,6);
                grid(i,   j-1, 7) = gridCopy( i,j,7);
                grid(i+1, j-1, 8) = gridCopy( i,j,8);
            }
        }
    }
    // periodic RB
    for (int i=1; i<lengthY-1; i++){
        // linke RB
        grid(lengthY-2, i, 4) = gridCopy( 0,i,4);
        grid(lengthY-2, i, 5) = gridCopy( 0,i,5);
        grid(lengthY-2, i, 6) = gridCopy( 0,i,6);
        // rechte RB
        grid(1, i, 1) = gridCopy( lengthX-1,i,1);
        grid(1, i, 2) = gridCopy( lengthX-1,i,2);
        grid(1, i, 8) = gridCopy( lengthX-1,i,8);
    }
    for (int i=0; i<lengthY; i++){
        // obere RB
        grid(i, lengthY-2, 6) = gridCopy( i,lengthY-1,2);
        grid(i, lengthY-2, 7) = gridCopy( i,lengthY-1,3);
        grid(i, lengthY-2, 8) = gridCopy( i,lengthY-1,4);
        // untere RB
        grid(i, 1, 4) = gridCopy( i,0,8);
        grid(i, 1, 3) = gridCopy( i,0,7);
        grid(i, 1, 2) = gridCopy( i,0,6);
    }

    // Cylinder cells
    for( int i=1; i<lengthX-1; i++ ){
        for( int j=1; j<lengthY-1; j++ ){
            if( isBoundary(i,j) ){
                grid(i-1, j,   5) = gridCopy( i,j,1);
                grid(i, j+1,   3) = gridCopy( i,j,7);
                grid(i, j-1,   7) = gridCopy( i,j,3);
                grid(i+1, j,   1) = gridCopy( i,j,5);
                grid(i-1, j-1, 6) = gridCopy( i,j,2);
                grid(i+1, j-1, 8) = gridCopy( i,j,4);
                grid(i+1, j+1, 2) = gridCopy( i,j,6);
                grid(i-1, j+1, 4) = gridCopy( i,j,8);
            }
        }
    }
}

void initBoundBoolean(){

    int lengthX = (int)grid.lengthX();
    int lengthY = (int)grid.lengthY();

    for(int i=0; i<lengthX; ++i){
        isBoundary(i,0) = 1;
        isBoundary(i,lengthY-1) = 1;
    }
    for(int i=0; i<lengthY; ++i){
        isBoundary(0,i) = 1;
        isBoundary(lengthX-1,i) = 1;
    }
    for(int i=0; i<lengthY; ++i){
        for(int j=0; j<lengthX; ++j){
            if((sqrt((i*spaceing-0.008)*(i*spaceing-0.008) + (j*spaceing-0.02)*(j*spaceing-0.02))) <= 0.0025 ){
                isBoundary(j,i) = 1;
            }
            else{
                isBoundary(j,i) = 0;
            }
        }
    }
}

void collideStep(){





}



void testStream(){



}


int main( int args, char** argv ){

    if( args != 2 ){
        cout << "USAGE: lbm scenario" << endl;
        exit( EXIT_SUCCESS );
    }

    if( strcmp( argv[1], "scenario1") == 0 ){
        timestep = TIMESTEP;
        spaceing = SPACEING;
        scenario1 = true;
        viscosity_l = 1e-6 * (timestep/(spaceing*spaceing));
        timeToSimulate = 3.0;
        acceleration = 0.01;
        resolution = 30;
        omega = 1.0 / (3.0*viscosity_l + 0.5 );
        //force_l =
        numCellsX = 0.06 / spaceing;
        numCellsY = 0.02 / spaceing;
        grid = grid_lattice<double>( numCellsX, numCellsY );
        gridCopy = grid_lattice<double>( numCellsX, numCellsY );
        isBoundary = grid_lattice<int>( numCellsX, numCellsY, 0, 1 );
    }else if( strcmp( argv[1], "scenario2") == 0){
        timestep = TIMESTEP;
        spaceing = SPACEING;
        scenario1 = false;
        viscosity_l = 1e-6 * (timestep/(spaceing*spaceing));
        timeToSimulate = 5.0;
        acceleration = 0.016;
        resolution = 60;
        omega = 1.0 / (3.0*viscosity_l + 0.5 );
        //force_l =
        numCellsX = 0.06 / spaceing;
        numCellsY = 0.02 / spaceing;
        grid = grid_lattice<double>( numCellsX, numCellsY );
        gridCopy = grid_lattice<double>( numCellsX, numCellsY );
        isBoundary = grid_lattice<int>( numCellsX, numCellsY, 0, 1 );
    }else{
        cout << "scenario has to be 1 or 2" << endl;
        exit( EXIT_SUCCESS );
    }



    cout << "timestep " << timestep << endl;
    cout << "spaceing " << spaceing << endl;
    cout << "force " << force_l << endl;
    cout << "omega " << omega << endl;



    exit( EXIT_SUCCESS );
}
