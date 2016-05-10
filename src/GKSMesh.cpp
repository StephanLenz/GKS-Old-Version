
#include "GKSMesh.h"
#include "Cell.h"
#include "Types.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>    //min()

using namespace std;

GKSMesh::GKSMesh()
{
}

GKSMesh::GKSMesh(Parameters param)
{
    this->param = param;
    this->iter = 0;
}


GKSMesh::~GKSMesh()
{
}

void GKSMesh::generateRectMesh(double lengthX, double lengthY, int nx, int ny)
{
	double dx = lengthX / (double)nx;
	double dy = lengthY / (double)ny;

	this->lengthX = lengthX;
	this->lengthY = lengthY;

	Cell*		tmpCell;
	Interface*  tmpInterface;
    float2      normal;

	//=========================================================================
	//=========================================================================
	//		Cell generation
	//			including ghost cells
	//=========================================================================
	//=========================================================================
    BoundaryCondition* currentBC = NULL;
	for (int i = -1; i < ny + 1; i++)       // Y-Direction
	{

		for (int j = -1; j < nx + 1; j++)   // X-Direction
		{
            if (j == -1)         currentBC = BoundaryConditionList[0];
            else if (j == nx)    currentBC = BoundaryConditionList[2];
            else if (i == -1)    currentBC = BoundaryConditionList[1];
            else if (i == ny)    currentBC = BoundaryConditionList[3];
            else                 currentBC = NULL;

			//                      cell centerX         cell centerY
			tmpCell = new Cell(((double)j + 0.5)*dx, ((double)i + 0.5)*dy, dx, dy, currentBC);
			// add interface to list
			this->CellList.push_back(tmpCell);
		}
	}

	//=========================================================================
	//=========================================================================
	//						F interface generation
	//=========================================================================
	//=========================================================================
    normal.x = 1;
    normal.y = 0;
	for (int i = 0; i <= ny + 1; i++)       // Y-Direction
	{
		for (int j = 0; j < nx + 1; j++)    // X-Direction
		{
			// create a new interface with the adjacent cells
			tmpInterface = new Interface(CellList[i*(nx + 2) + j], CellList[i*(nx + 2) + (j + 1)], 0, normal);
			// add itnerface to list
			this->InterfaceList.push_back(tmpInterface);
		}
	}

	//=========================================================================
	//=========================================================================
	//						G interface generation
	//=========================================================================
	//=========================================================================
    normal.x = 0;
    normal.y = 1;
	for (int i = 0; i < ny + 1; i++)        // Y-Direction
	{
		for (int j = 0; j <= nx + 1; j++)   // X-Direction
		{
			// create a new interface with the adjacent cells
			tmpInterface = new Interface(CellList[i*(nx + 2) + j], CellList[(i + 1)*(nx + 2) + j], 1, normal);
			// add itnerface to list
			this->InterfaceList.push_back(tmpInterface);
		}
	}


	return;
}

void GKSMesh::initMeshConstant(double rho, double u, double v, double T)
{
	for (vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
	{
		(*i)->setValues(rho, u, v, T);
	}
}

void GKSMesh::initMeshLinearTemperature(double rho, double u, double v, double* T)
{
	// Temprature definition
	//    ------------
	//    |   T[1]   |
	//    |          |
	//    |   T[0]   |
	//    ------------
	double interpolatedT;
	float2 center;
	for (vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
	{
		center = (*i)->getCenter();
		//interpolatedT = ( center.y*T[0] + (lengthY - center.y)*T[1] ) / lengthY;
        interpolatedT = T[0] + center.y*(T[1] - T[0]) / this->lengthY;
		(*i)->setValues(rho, u, v, interpolatedT);
	}
}

void GKSMesh::addBoundaryCondition( int rhoType, int UType, int VType, int TType, 
                                    double rho, double U, double V, double T)
{
    BoundaryCondition* tmp = new BoundaryCondition( rhoType, UType, VType, TType,
                                                    rho, U, V, T);
    BoundaryConditionList.push_back(tmp);
}

void GKSMesh::applyBoundaryCondition()
{
    for (vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i)
    {
        if ( ( (*i)->isGhostCell() ) )
            (*i)->applyBoundaryCondition();
    }
}

void GKSMesh::computeGlobalTimestep()
{
    this->dt = 1.0e99;
    for (vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i)
    {
        if (!((*i)->isGhostCell()))
        {
            this->dt = min( (*i)->getLocalTimestep(this->param), this->dt );
        }
    }
    this->dt *= this->param.CFL;
}

void GKSMesh::timeStep()
{
    this->iter++;

    if (this->param.verbose) cout << "Iterration: " << this->iter << endl;

    if(this->param.verbose) cout << "  Compute Timestep ..." << endl;
    this->computeGlobalTimestep();
    if (this->param.verbose) cout << "    dt = " << this->dt << endl;

    // ========================================================================

    if (this->param.verbose) cout << "  Apply Boundary Conditions ..." << endl;
    this->applyBoundaryCondition();

    // ========================================================================

    if (this->param.verbose) cout << "  Compute Mass- and Momentum Fluxes ..." << endl;
    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
        if( !(*i)->isGhostInterface() )
            (*i)->computeMassMomentumFlux(this->dt, this->param.tauMassMomentum);
    }

    if (this->param.verbose) cout << "  Update Mass and Momentum in Cells ..." << endl;
    for (vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i)
    {
        if( !(*i)->isGhostCell() )
            (*i)->updateMassMomentum(this->dt, this->param.G0, this->param.beta, this->param.TAve);
    }

    // ========================================================================

    if (this->param.verbose) cout << "  Apply Boundary Conditions ..." << endl;
    this->applyBoundaryCondition();

    // ========================================================================

    if (this->param.verbose) cout << "  Compute Heat Flux ..." << endl;
    for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
    {
        if ( ! (*i)->isGhostInterface() )
            (*i)->computeHeatFlux(this->dt, this->param.tauHeat);
    }

    if (this->param.verbose) cout << "  Update Temperature in Cells" << endl;
    for (vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i)
    {
        if ( ! (*i)->isGhostCell() )
            (*i)->updateTemperature();
    }

    // ========================================================================

    if (this->param.verbose) cout << "  Apply Boundary Conditions ..." << endl;
    this->applyBoundaryCondition();
    
    
}

void GKSMesh::iterate()
{
    ostringstream filename;
    filename << "out/result_0.vtk";
    writeVTKFile(filename.str(), true, false);

    while (this->iter < this->param.numberOfIterations)
    {
        this->timeStep();

        if (this->iter%this->param.outputInterval == 0)
        {
            ostringstream filename;
            filename << "out/result_" << this->iter << ".vtk";
            writeVTKFile(filename.str(), true, false);
        }
        //cout << this->toString();
        //cout << this->cellValuesToString();
    }
}

string GKSMesh::toString()
{
	ostringstream tmp;
	tmp << "The Mesh has following interfaces:\n";

	for (vector<Interface*>::iterator i = InterfaceList.begin(); i != InterfaceList.end(); ++i)
	{
        if( ! (*i)->isGhostInterface() )
		tmp << (*i)->toString();
	}

	return tmp.str();
}

string GKSMesh::cellValuesToString()
{
    ostringstream tmp;
    for (vector<Cell*>::iterator i = this->CellList.begin(); i != this->CellList.end(); ++i)
    {
        if( ! (*i)->isGhostCell()  )
            tmp << (*i)->valuesToString() << "\n";
    }
    return tmp.str();
}

void GKSMesh::writeVTKFile(string filename, bool data, bool BC)
{
    cout << "Wrinting file " << filename << " ... ";
	// open file stream
	ofstream file;
	file.open(filename.c_str());

	if (!file.is_open()) {
		cout << " File cound not be opened.\n\nERROR TERMINATION!\n\n\n";
		return;
	}

    this->writeGeometry(file);

    if (data) this->writeData(file);

    file.close();

    cout << "done!" << endl;
}

void GKSMesh::writeGeometry(ofstream& file)
{

    // write VTK Header
    file << "# vtk DataFile Version 1.0\n";
    file << "by Stephan Lenz\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    // write nodes
    //( one dummy node with the ID 0 must be written )
    file << "POINTS " << 4 * this->CellList.size() + 1 << " float\n";
    file << "0.0 0.0 0.0 \n";

    for (vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i)
    {
        file << (*i)->writeNodes();
    }

    // write elements
    file << "CELLS " << this->CellList.size() << " " << 5 * this->CellList.size() << endl;
    for (int i = 0; i < this->CellList.size(); ++i)
    {
        file << 4 << " " << i * 4 + 1
            << " " << i * 4 + 2
            << " " << i * 4 + 3
            << " " << i * 4 + 4 << endl;
    }

    // write element tyes( 9 = quad element )
    file << "CELL_TYPES " << this->CellList.size() << endl;
    for (int i = 0; i < this->CellList.size(); ++i)
    {
        file << "9" << endl;
    }
}

void GKSMesh::writeData(ofstream& file)
{
    // write cell data ( ID and stress )
    file << "CELL_DATA " << this->CellList.size() << endl;
    file << "FIELD Lable 5\n";

    file << "rho 1 " << this->CellList.size() << " float\n";
    for (vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i)
    {
        file << (*i)->getPrim().rho << endl;
    }
    file << "U 1 " << this->CellList.size() << " float\n";
    for (vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i)
    {
        file << (*i)->getPrim().U << endl;
    }
    file << "V 1 " << this->CellList.size() << " float\n";
    for (vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i)
    {
        file << (*i)->getPrim().V << endl;
    }
    file << "T 1 " << this->CellList.size() << " float\n";
    for (vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i)
    {
        file << (*i)->getPrim().T << endl;
    }
    file << "GhostCell 1 " << this->CellList.size() << " int\n";
    for (vector<Cell*>::iterator i = CellList.begin(); i != CellList.end(); ++i)
    {
        if ((*i)->isGhostCell())
            file << 1 << endl;
        else
            file << 0 << endl;
    }
}

