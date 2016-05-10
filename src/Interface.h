
#ifndef INTERFACE_H
#define INTERFACE_H

#include "Cell.h"
//class Cell;
#include <string>

using namespace std;

class Interface
{
private:
	Cell* negCell;
	Cell* posCell;

    float2 normal;
    int axis;

    double MassMomentumFlux[3];
    double HeatFlux;
public:
	Interface();
	Interface(Cell* negCell, Cell* posCell, int axis, float2 normal);
	~Interface();

	void computeMassMomentumFlux(double dt, double tau);
    void computeHeatFlux(double dt, double tau);

    Cell* getNeigborCell(Cell* askingCell);
    ConservedVariable getMassMomentumFlux();
    double getHeatFlux();

    bool isGhostInterface();

	string toString();

private:

    void interpolatePrim(double* prim);
    void differentiateCons(double* normalGradCons, double* tangentialGradCons, double* prim);
    void differentiateTemperature(double& normalGradTemperatur, double& tangentialGradTemperatur);

    void rotate(double* vector);

    void computeMicroSlope(double* prim, double* macroSlope, double* microSlope);
    void computeMomentU(double* prim, double* MomentU, double* MomentV, int numberMoments);
    void computeMomentUV(double* MomentU, double* MomentV, int alpha, int beta, double* MomentUV);
    void computeMomentAU(double* a, double* MomentU, double* MomentV, int alpha, int beta, double* MomentAU);
};

#endif