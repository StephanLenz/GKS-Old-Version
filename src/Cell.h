
#ifndef RECTCELL2D_H
#define RECTCELL2D_H

//#include "Interface.h"
class Interface;
#include "Types.h"
#include "BoundaryCondition.h"
#include <string>

using namespace std;

class Cell
{
private:
	// Cell center
	double centerX;
	double centerY;

	// Cell size
	double dx;
	double dy;

    FluidParameter fluidParam;

	// links to interfaces
	//    -----------
	//    |    3    |
	//    | 0     2 |
	//    |    1    |
	//    -----------
	Interface** InterfaceList;

	// Primary Variables
	double prim[4];
    double primOld[4];

	// Conseved Variables
	double cons[3];
    double consOld[3];

    // Boundary Cell
    BoundaryCondition* BoundaryContitionPointer;

public:
	Cell();
	Cell(double centerX, double centerY, double dx, double dy, BoundaryCondition* BC);

	~Cell();

	void addInterface(Interface* newInterface, int direction);

	void update(double dt, double G0, double beta, double Tave);

    void storeOldValues();

    void applyBoundaryCondition();

	void setValues(double rho, double u, double v, double T);

    void computePrim();

    void computeCons();

    double getLocalTimestep(double nu);

	float2 getCenter();

    PrimaryVariable getPrim();
    PrimaryVariable getPrimOld();

    ConservedVariable getCons();
    ConservedVariable getConsOld();

    float2 getDx();

    Cell* getNeighborCell(int i);

    Cell* findNeighborInDomain();

    bool isGhostCell();

	string toString();

    string valuesToString();

	string writeNodes();
};

#endif