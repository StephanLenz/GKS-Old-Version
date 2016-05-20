// ============================================================================
//
//                      Compressible Thermal GKS
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
    //                  Driven-Cavity
    //
    // ========================================================================

    Parameters param;

    double H = 1.0;
    double W = 1.0;

    param.numberOfIterations = 100;
    param.outputInterval = 1;
    param.CFL = 0.5;

    param.verbose = false;

    // ========================================================================

    FluidParameter fluidParam;

    fluidParam.K  = 0;
    fluidParam.nu = 1e-3;
    fluidParam.R = 200.0;
    
    double uTop = 0.01;
    double TAve = 10.0;

    // ========================================================================

    GKSMesh* mesh = new GKSMesh(param, fluidParam);

    // Define Boundary Conditions
    //    -----------
    //    |    3    |
    //    | 0     2 |
    //    |    1    |
    //    -----------
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, 0.0,  0.0, 0.0);
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, 0.0,  0.0, 0.0);
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, 0.0,  0.0, 0.0);
    mesh->addBoundaryCondition(1, 0, 0, 1, 0.0, uTop, 0.0, 0.0);

    // Generate Mesh
    mesh->generateRectMesh(W, H, 10, 10);

    // Initialize Values
    mesh->initMeshConstant(10.0, 0.0, 0.0, TAve);

    // ========================================================================
    // ========================================================================
    // ========================================================================
    // ========================================================================
    // ========================================================================

    mesh->applyBoundaryCondition();

    mesh->iterate();
    //char a; cin >> a;
}