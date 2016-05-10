
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

	void updateMassMomentum(double dt, double G0, double beta, double Tave);
    void updateTemperature();

    void storeOldValues();

    void applyBoundaryCondition();

	void setValues(double rho, double u, double v, double T);

    void computePrim();

    void computeCons();

    double getLocalTimestep(Parameters param);

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

	string writeNodes();
};

#endif