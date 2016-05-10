
#include "Cell.h"
#include "Interface.h"
#include <sstream>
#include <iostream>
#include <algorithm>    // min()

using namespace std;

Cell::Cell()
{
	this->InterfaceList = new Interface*[4];
    memset(InterfaceList, NULL, 4 * sizeof(Interface*));
}

Cell::Cell(double centerX, double centerY, double dx, double dy, BoundaryCondition* BC)
{
	this->InterfaceList = new Interface*[4];
    memset(InterfaceList, NULL, 4 * sizeof(Interface*));

	this->centerX = centerX;
	this->centerY = centerY;

	this->dx = dx;
	this->dy = dy;

    this->BoundaryContitionPointer = BC;
}


Cell::~Cell()
{
	delete [] InterfaceList;
}

void Cell::updateMassMomentum(double dt, double G0, double beta, double Tave)
{
    this->storeOldValues();

    // negative interfaces = in flux
    // positive interfaces = out flux
    this->cons[0] += ( this->InterfaceList[0]->getMassMomentumFlux().rho
                     + this->InterfaceList[1]->getMassMomentumFlux().rho
                     - this->InterfaceList[2]->getMassMomentumFlux().rho
                     - this->InterfaceList[3]->getMassMomentumFlux().rho
                     ) / (this->dx*this->dy);

    this->cons[1] += ( this->InterfaceList[0]->getMassMomentumFlux().rhoU
                     + this->InterfaceList[1]->getMassMomentumFlux().rhoU
                     - this->InterfaceList[2]->getMassMomentumFlux().rhoU
                     - this->InterfaceList[3]->getMassMomentumFlux().rhoU
                     ) / (this->dx*this->dy);

    this->cons[2] += ( this->InterfaceList[0]->getMassMomentumFlux().rhoV
                     + this->InterfaceList[1]->getMassMomentumFlux().rhoV
                     - this->InterfaceList[2]->getMassMomentumFlux().rhoV
                     - this->InterfaceList[3]->getMassMomentumFlux().rhoV
                     ) / (this->dx*this->dy)
                   + this->getPrimOld().rho*beta*G0*(this->getPrimOld().T - Tave)*dt;

    this->computePrim();
}

void Cell::updateTemperature()
{
    this->prim[4] = ( this->InterfaceList[0]->getHeatFlux()
                    + this->InterfaceList[1]->getHeatFlux()
                    - this->InterfaceList[2]->getHeatFlux()
                    - this->InterfaceList[3]->getHeatFlux()
                    ) / (this->dx*this->dy);
}

void Cell::storeOldValues()
{
    for (int i = 0; i < 4; i++)
        primOld[i] = prim[i];

    for (int i = 0; i < 3; i++)
        consOld[i] = cons[i];
}

void Cell::applyBoundaryCondition()
{

    Cell* neighborCell = this->findNeighborInDomain();

    // if no neighbor was found
    // TODO: Corner cell, no BC implemented right now.....
    if (neighborCell == NULL)
    {
        // search any inerface and take the neigbor 
        for (int i = 0; i < 4; i++)
        {
            if (this->InterfaceList[i] != NULL)
            {
                neighborCell = InterfaceList[i]->getNeigborCell(this)->findNeighborInDomain();
                break;
            }
        }
    }

    // loop over all primitive variables
    for (int i = 0; i < 4; i++)
    {
        int type =     this->BoundaryContitionPointer->getType(i);
        double value = this->BoundaryContitionPointer->getValue(i);

        if (type == 0)
        {
            this->prim[i] = 2*value - neighborCell->prim[i];
        }
        else if (type == 1)
        {
            this->prim[i] = neighborCell->prim[i];
        }
    }
    
    this->computeCons();
}

void Cell::addInterface(Interface* newInterface, int direction)
{
	this->InterfaceList[direction] = newInterface;
}

void Cell::setValues(double rho, double u, double v, double T)
{
	this->prim[0] = rho;
	this->prim[1] = u;
	this->prim[2] = v;
	this->prim[3] = T;

    this->computeCons();
}

void Cell::computePrim()
{
    this->prim[0] = this->cons[0];
    this->prim[1] = this->cons[1] / this->cons[0];
    this->prim[2] = this->cons[2] / this->cons[0];
}

void Cell::computeCons()
{
    this->cons[0] = this->prim[0];
    this->cons[1] = this->prim[0] * this->prim[1]; 
    this->cons[2] = this->prim[0] * this->prim[2];
}

double Cell::getLocalTimestep(Parameters param)
{
    double velocitySquare = this->getPrim().U*this->getPrim().U
                          + this->getPrim().V*this->getPrim().V;
    double localTimestep =  min(dx, dy) 
                         / ( velocitySquare
                           + 1.0/sqrt(3.0) 
                           + 2.0*param.nu/min(dx, dy) );
    return localTimestep;
}

float2 Cell::getCenter()
{
	float2 center;
	center.x = this->centerX;
	center.y = this->centerY;
	return center;
}

PrimaryVariable Cell::getPrim()
{
    PrimaryVariable tmp;
    tmp.rho = this->prim[0];
    tmp.U   = this->prim[1];
    tmp.V   = this->prim[2];
    tmp.T   = this->prim[3];
    return tmp;
}

PrimaryVariable Cell::getPrimOld()
{
    PrimaryVariable tmp;
    tmp.rho = this->primOld[0];
    tmp.U   = this->primOld[1];
    tmp.V   = this->primOld[2];
    tmp.T   = this->primOld[3];
    return tmp;
}

ConservedVariable Cell::getCons()
{
    ConservedVariable tmp;
    tmp.rho  = this->cons[0];
    tmp.rhoU = this->cons[1];
    tmp.rhoV = this->cons[2];
    return tmp;
}

ConservedVariable Cell::getConsOld()
{
    ConservedVariable tmp;
    tmp.rho  = this->consOld[0];
    tmp.rhoU = this->consOld[1];
    tmp.rhoV = this->consOld[2];
    return tmp;
}

float2 Cell::getDx()
{
    float2 tmp;
    tmp.x = this->dx;
    tmp.y = this->dy;
    return tmp;
}

Cell * Cell::getNeighborCell(int i)
{
    return this->InterfaceList[i]->getNeigborCell(this);
}

Cell * Cell::findNeighborInDomain()
{
    Cell* neighborCell = NULL;
    // find Interface to domain
    for (int i = 0; i < 4; i++)
    {
        if ((this->InterfaceList[i] != NULL) && !this->InterfaceList[i]->getNeigborCell(this)->isGhostCell())
            neighborCell = this->InterfaceList[i]->getNeigborCell(this);
    }
    return neighborCell;
}

bool Cell::isGhostCell()
{
    return !(BoundaryContitionPointer== NULL);
}

string Cell::toString()
{
	ostringstream tmp;
	tmp << "Center = ( ";
	tmp << this->centerX;
	tmp << " , ";
	tmp << this->centerY;
	tmp << " )";
	return tmp.str();
}

string Cell::writeNodes()
{
	ostringstream tmp;
	tmp << this->centerX - 0.5*dx << " ";
	tmp << this->centerY - 0.5*dy << " 0.0\n";

	tmp << this->centerX + 0.5*dx << " ";
	tmp << this->centerY - 0.5*dy << " 0.0\n";

	tmp << this->centerX + 0.5*dx << " ";
	tmp << this->centerY + 0.5*dy << " 0.0\n";

	tmp << this->centerX - 0.5*dx << " ";
	tmp << this->centerY + 0.5*dy << " 0.0\n";
	return tmp.str();
}
