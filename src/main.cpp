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