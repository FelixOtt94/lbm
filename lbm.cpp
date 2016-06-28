#include <iostream>
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <vector>

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

static grid_lattice<double> grid;
static grid_lattice<double> gridCopy;
static grid_lattice<bool> isBoundary;


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
    grid_lattice(size_t xDim, size_t yDim,  size_t f=9, T value) {
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
        return data[j * (lengthInX*lengthInF) + i*lengthInF +f];
    }

    void getCopy( grid_lattice &source ){
        for( int i=0; i<lengthInX; i++){
            for( int j=0; j<lengthInY; j++){
                for( int f=0; f<lengthInF; f++){
                    data[j * (lengthInX*lengthInF) + i*lengthInF +f] = source( i,j,f);
                }
            }
        }
    }

};


void streamStep(){
    gridCopy.getCopy( grid );

    for( int i=0; i<grid.lengthX(); i++ ){
        for( int j=0; j<grid.lengthY(); j++ ){
            if( isBoundary(i,j) ){
                // HANDLE BOUNDARY CONDITIONS
            }else{
                grid(i+1, j,   1) = girdCopy( i,j,1);
                grid(i+1, j+1, 2) = girdCopy( i,j,2);
                grid(i,   j+1, 3) = girdCopy( i,j,3);
                grid(i-1, j+1, 4) = girdCopy( i,j,4);
                grid(i-1, j,   5) = girdCopy( i,j,5);
                grid(i-1, j-1, 6) = girdCopy( i,j,6);
                grid(i,   j-1, 7) = girdCopy( i,j,7);
                grid(i+1, j-1, 8) = girdCopy( i,j,8);
            }
        }
    }

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
        viscosity_l = 1e-6 * (timestep/(spacing*spacing));
        timeToSimulate = 3.0;
        acceleration = 0.01;
        resolution = 30;
        omega = 1.0 / (3.0*viscosity_l + 0.5 );
        //force_l =
        numCellsX = 0.06 / spaceing;
        numCellsY = 0.02 / spaceing;
        grid = grid_lattice::grid_lattice( numCellsX, numCellsY );
        gridCopy = grid_lattice::grid_lattice( numCellsX, numCellsY );
        isBoundary = grid_lattice::grid_lattice( numCellsX, numCellsY, 1 );
    }else if( strcmp( argv[1], "scenario2") == 0){
        timestep = TIMESTEP;
        spaceing = SPACEING;
        scenario1 = false;
        viscosity_l = 1e-6 * (timestep/(spacing*spacing));
        timeToSimulate = 5.0;
        acceleration = 0.016;
        resolution = 60;
        omega = 1.0 / (3.0*viscosity_l + 0.5 );
        //force_l =
        numCellsX = 0.06 / spaceing;
        numCellsY = 0.02 / spaceing;
        grid = grid_lattice::grid_lattice( numCellsX, numCellsY );
        gridCopy = grid_lattice::grid_lattice( numCellsX, numCellsY );
        isBoundary = grid_lattice::grid_lattice( numCellsX, numCellsY, 1 );
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
