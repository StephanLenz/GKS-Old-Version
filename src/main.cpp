
#include "GKSMesh.h"
#include "BoundaryCondition.h"
#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{

	double H = 1.0;
	double W = 2.0;


    // ========================================================================
    Parameters param;
	
	param.numberOfIterations = 1000000;
    param.outputInterval = 10000;
	param.CFL = 0.5;

	param.Pr = 1.0;
	param.Ra = 50000.0;

	param.G0 = 1.0;
	param.beta = 0.1;

	param.Tbot = 1.0;
	param.Ttop = 0.0;
	param.TAve = 0.5*(param.Tbot + param.Ttop);
	param.TInf = param.Ttop;

	param.rhoInf = 1.0;
	param.pInf = param.rhoInf / 3.0;

	param.nu = sqrt(param.beta*(param.Tbot - param.Ttop)*H*H*H*param.Pr / param.Ra);
	param.k = param.nu / param.Pr;
	param.uInf = param.k / H;

	param.Re = param.uInf*H / param.nu;
	param.Ma = sqrt(3.0)*param.uInf;

    param.tauMassMomentum = 3.0*param.nu;
    param.tauHeat         = 3.0*param.k;

    param.verbose = false;

    // ========================================================================

    GKSMesh* mesh = new GKSMesh(param);

    // Define Boundary Conditions
    //    -----------
    //    |    3    |
    //    | 0     2 |
    //    |    1    |
    //    -----------
    ///*
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, 0.0, 0.0, 0.0);
    mesh->addBoundaryCondition(1, 0, 0, 0, 0.0, 0.0, 0.0, param.Tbot);
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, 0.0, 0.0, 0.0);
    mesh->addBoundaryCondition(1, 0, 0, 0, 0.0, 0.0, 0.0, param.Ttop);
    //*/
    /*
    mesh->addBoundaryCondition(1, 0, 1, 1, 0.0, 1.0, 0.0, 10.0); // left
    mesh->addBoundaryCondition(1, 0, 1, 1, 0.0, 0.0, 0.0, 10.0); // bottom
    mesh->addBoundaryCondition(1, 0, 1, 1, 0.0, 1.0, 0.0, 10.0); // right
    mesh->addBoundaryCondition(1, 0, 1, 1, 0.0, 0.0, 0.0, 10.0); // top
    */

    // Generate Mesh
    mesh->generateRectMesh(W, H, 80, 40);
    //mesh->generateRectMesh(H, W, 40, 80);

    //cout << mesh->toString();

    // Initialize Values
	double T[2] = { param.Tbot, param.Ttop };
	mesh->initMeshLinearTemperature(param.rhoInf, 0.0, 0.0, T);
    //mesh->initMeshConstant(param.rhoInf, 1.0, 0.0, 10.0);
    
    mesh->applyBoundaryCondition();

    mesh->iterate();

	//char a; cin >> a;
}