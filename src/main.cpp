// ============================================================================
//
//                      Incompressible Thermal GKS
//
// ============================================================================
#include "GKSMesh.h"
#include "BoundaryCondition.h"
#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
    // ========================================================================
    //
    //                  Rayleigh-Bernard
    //
    // ========================================================================
    /*
    Parameters param;

    double H = 1.0;
    double W = 2.0;
    	
	param.numberOfIterations = 1000000;
    param.outputInterval = 10000;
    param.CFL = 0.5;
    
	double Pr = 1.0;
	double Ra = 50000.0;
    
    param.G0 = 1.0;
    param.beta = 0.1;
    
	double Tbot = 1.0;
	double Ttop = 0.0;
	param.TAve = 0.5*(Tbot + Ttop);
	double TInf = Ttop;
    
	double rhoInf = 1.0;
	double pInf = rhoInf / 3.0;
    
	param.nu = sqrt(param.beta*(Tbot - Ttop)*H*H*H*Pr / Ra);
	param.k = param.nu / Pr;
	double uInf = param.k / H;
    
	double Re = uInf*H / param.nu;
	double Ma = sqrt(3.0)*uInf;
    
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
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, 0.0, 0.0, 0.0);
    mesh->addBoundaryCondition(1, 0, 0, 0, 0.0, 0.0, 0.0, Tbot);
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, 0.0, 0.0, 0.0);
    mesh->addBoundaryCondition(1, 0, 0, 0, 0.0, 0.0, 0.0, Ttop);

    // Generate Mesh
    mesh->generateRectMesh(W, H, 80, 40);

    // Initialize Values
	double T[2] = { Tbot, Ttop };
	mesh->initMeshLinearTemperature(rhoInf, 0.0, 0.0, T);
    */


    // ========================================================================
    //
    //                  Driven-Cavity
    //
    // ========================================================================

    Parameters param;

    double H = 1.0;
    double W = 1.0;

    param.numberOfIterations = 100;
    param.outputInterval = 1;
    param.CFL = 0.5;
    
    param.G0 = 1.0;
    param.beta = 0.1;

    param.TAve = 1.0;

    param.nu = 0.1;
    param.k = 0.1;
    
    double uTop = 0.01;
    double Re = uTop*H / param.nu;

    param.tauMassMomentum = 3.0*param.nu;
    param.tauHeat = 3.0*param.k;

    param.verbose = false;

    // ========================================================================

    GKSMesh* mesh = new GKSMesh(param);

    // Define Boundary Conditions
    //    -----------
    //    |    3    |
    //    | 0     2 |
    //    |    1    |
    //    -----------
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, 0.0,  0.0, param.TAve);
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, 0.0,  0.0, param.TAve);
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, 0.0,  0.0, param.TAve);
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, uTop, 0.0, param.TAve);

    // Generate Mesh
    mesh->generateRectMesh(W, H, 10, 10);

    // Initialize Values
    mesh->initMeshConstant(1.0, 0.0, 0.0, param.TAve);

    // ========================================================================
    // ========================================================================
    // ========================================================================
    // ========================================================================
    // ========================================================================

    mesh->applyBoundaryCondition();

    mesh->iterate();
    //char a; cin >> a;
}