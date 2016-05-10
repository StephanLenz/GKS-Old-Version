
#ifndef GKSMESH_H
#define GKSMESH_H

#include "Cell.h"
#include "Interface.h"
#include "BoundaryCondition.h"
#include "Types.h"
#include <vector>
#include <string>

using namespace std;

using namespace std;

class GKSMesh
{
private:
	vector<Cell*>		CellList;
	vector<Interface*>	InterfaceList;
    vector<BoundaryCondition*> BoundaryConditionList;

    Parameters param;

	double lengthX;		
	double lengthY;

    double dt;
    unsigned int iter;

public:
	GKSMesh();

    GKSMesh(Parameters param);

	~GKSMesh();

	void generateRectMesh(double lengthX, double lengthY, int nx, int ny);

	void initMeshConstant(double rho, double u, double v, double T);

	void initMeshLinearTemperature(double rho, double u, double v, double * T);

    void addBoundaryCondition(  int rhoType, int UType, int VType, int TType,
                                double rho, double U, double V, double T);

    void applyBoundaryCondition();

    void computeGlobalTimestep();

    void timeStep();

    void iterate();

	string toString();

	void writeVTKFile(string filename, bool data = true, bool BC = false);

private:

    void writeGeometry(ofstream& file);

    void writeData(ofstream& file);
};

#endif