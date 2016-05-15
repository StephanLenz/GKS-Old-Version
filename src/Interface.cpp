

#include "Interface.h"
#include <sstream>

Interface::Interface()
{
}

Interface::Interface(Cell* negCell, Cell* posCell, int axis, float2 normal)
{
	this->negCell = negCell;
	this->posCell = posCell;
    
    // links to interfaces
    //    ---------------------
    //    |    3    |    3    |
    //    | 0     2 | 0     2 |
    //    |    1    |    1    |
    //    ---------------------
    //     neg Cell   pos Cell
    //
    //    -----------
    //    |    3    |
    //    | 0     2 |   pos Cell
    //    |    1    |
    //    -----------
    //    |    3    |
    //    | 0     2 |   neg Cell
    //    |    1    |
    //    -----------

	this->negCell->addInterface(this,axis+2);
	this->posCell->addInterface(this,axis+0);

    this->normal = normal;
    this->axis = axis;
}

Interface::~Interface()
{
}

void Interface::computeFlux(double dt, double tau)
{
    const int NUMBER_OF_MOMENTS = 5;

    double prim[4];
    double normalGradCons[4];
    double tangentialGradCons[4];
    double timeGrad[4];

    double a[4];
    double b[4];
    double A[4];

    double MomentU[NUMBER_OF_MOMENTS];
    double MomentV[NUMBER_OF_MOMENTS];
    double MomentXi[NUMBER_OF_MOMENTS];

    double MomentUV_1_0[3];

    double normalMomentAU_1_0[3];
    double normalMomentAU_2_0[3];
    double tangentialMomentAU_0_1[3];
    double tangentialMomentAU_1_1[3];
    double timeMomentAU_1_0[3];

    // compute the length of the interface
    double dy = this->posCell->getDx().x * normal.y
              + this->posCell->getDx().y * normal.x;

    // time integration Coefficients
    double timeCoefficients[6] = {  dt + tau * exp(-dt/tau) ,
                                  - dt * tau * ( 1.0 + exp(-dt/tau),
                                    0.5 * dt*dt - dt*tau - exp(-dt / tau),
                                  - tau * exp(-dt / tau),
                                  - tau * (1 +2.0*tau) * exp(-dt / tau),
                                    tau * tau * exp(-dt / tau) };

    this->interpolatePrim(prim);
    this->differentiateCons(normalGradCons, tangentialGradCons, prim);

    // in case of horizontal interface (G interface), swap velocity directions
    if (this->axis == 1)
    {
        this->rotate(prim);
        this->rotate(normalGradCons);
        this->rotate(tangentialGradCons);
    }

    // spacial micro slopes a = a1 + a2 u + a3 v
    //                      b = b1 + b2 u + b3 v

    this->computeMicroSlope(prim, normalGradCons,     a);
    this->computeMicroSlope(prim, tangentialGradCons, b);

    // temporal micro slopes A = A1 + A2 u + A3 v

    this->computeMoments(prim, MomentU, MomentV, MomentXi, NUMBER_OF_MOMENTS);

    this->computeTimeDerivative(prim, MomentU, MomentV, MomentXi, a, b, timeGrad);

    this->computeMicroSlope(prim, timeGrad, A);

    // compute mass and momentum fluxes

    this->assembleFlux(MomentU, MomentV, MomentXi, a, b, A, timeCoefficients);

    // in case of horizontal interface (G interface), swap velocity fluxes
    if (this->axis == 1)
        this->rotate(this->Flux);

}

Cell * Interface::getNeigborCell(Cell * askingCell)
{
    if (posCell == askingCell)
        return negCell;
    else
        return posCell;
}

ConservedVariable Interface::getFlux()
{
    ConservedVariable tmp;
    tmp.rho  = this->Flux[0];
    tmp.rhoU = this->Flux[1];
    tmp.rhoV = this->Flux[2];
    tmp.rhoE = this->Flux[3];
    return tmp;
}

bool Interface::isGhostInterface()
{
    return this->posCell->isGhostCell() && this->negCell->isGhostCell();
}

string Interface::toString()
{
	ostringstream tmp;
	tmp << "Interface between: \n";
	tmp << this->negCell->toString();
	tmp << "\n";
	tmp << this->posCell->toString();
	tmp << "\n";
    tmp << this->Flux[0] << " " << this->Flux[1] << " " << this->Flux[2] << " " << this->Flux[3];
    tmp << "\n";
    tmp << "\n";
	return tmp.str();
}

void Interface::interpolatePrim(double * prim)
{

    prim[0] = 0.5*( this->negCell->getPrim().rho
                  + this->posCell->getPrim().rho );
    prim[1] = 0.5*( this->negCell->getPrim().U
                  + this->posCell->getPrim().U   );
    prim[2] = 0.5*( this->negCell->getPrim().V
                  + this->posCell->getPrim().V   );
    prim[3] = 0.5*( this->negCell->getPrim().L
                  + this->posCell->getPrim().L   );
}

void Interface::differentiateCons(double* normalGradCons, double* tangentialGradCons, double* prim)
{
    // ========================================================================
    // normal direction
    // ========================================================================

    // compute the distance between 
    double dn = ( ( this->posCell->getDx().x + this->negCell->getDx().x ) * normal.x
                + ( this->posCell->getDx().y + this->negCell->getDx().y ) * normal.y ) * 0.5;

    normalGradCons[0] = ( this->posCell->getCons().rho
                        - this->negCell->getCons().rho  ) / dn;

    normalGradCons[1] = ( this->posCell->getCons().rhoU
                        - this->negCell->getCons().rhoU ) / dn;

    normalGradCons[2] = ( this->posCell->getCons().rhoV
                        - this->negCell->getCons().rhoV ) / dn;

    normalGradCons[3] = ( this->posCell->getCons().rhoE
                        - this->negCell->getCons().rhoE ) / dn;

    // ========================================================================
    // tangential direction
    // ========================================================================

    // The tangential derivative is computed by finite difference between the
    // values at the edge of the interface (A, B in fig).
    // These are computed by interpolation.
    //
    //  A = 0.5 (pos + pos pos + neg + neg pos)
    //
    //  ---------------------------------
    //  |               |               |
    //  |    neg pos    |    pos pos    |
    //  |               |               |
    //  --------------- A --------------
    //  |               |               |
    //  |               |               |
    //  |      neg      |      pos      |
    //  |               |               |
    //  |               |               |
    //  --------------- B ---------------
    //  |               |               |
    //  |    neg neg    |    pos neg    |
    //  |               |               |
    //  ---------------------------------

    // get the indieces of the perpendicular interfaces for tangential derivative
    int posIdx;
    int negIdx;
    if (this->axis == 0)
    {
        posIdx = 3;
        negIdx = 1;
    }
    else
    {
        posIdx = 2;
        negIdx = 0;
    }

    // compute the tangential distance (length of the interface)
    double dt = this->posCell->getDx().x * normal.y
              + this->posCell->getDx().y * normal.x;

    tangentialGradCons[0] = ( ( this->posCell->getNeighborCell(posIdx)->getCons().rho
                              + this->negCell->getNeighborCell(posIdx)->getCons().rho
                              + this->posCell->getCons().rho 
                              + this->negCell->getCons().rho 
                              ) * 0.25
                            - ( this->posCell->getNeighborCell(negIdx)->getCons().rho
                              + this->negCell->getNeighborCell(negIdx)->getCons().rho 
                              + this->posCell->getCons().rho 
                              + this->negCell->getCons().rho
                              ) * 0.25
                            ) / dt;

    tangentialGradCons[1] = ( ( this->posCell->getNeighborCell(posIdx)->getCons().rhoU
                              + this->negCell->getNeighborCell(posIdx)->getCons().rhoU
                              + this->posCell->getCons().rhoU 
                              + this->negCell->getCons().rhoU 
                              ) * 0.25
                            - ( this->posCell->getNeighborCell(negIdx)->getCons().rhoU
                              + this->negCell->getNeighborCell(negIdx)->getCons().rhoU 
                              + this->posCell->getCons().rhoU 
                              + this->negCell->getCons().rhoU 
                              ) * 0.25
                            ) / dt;

    tangentialGradCons[2] = ( ( this->posCell->getNeighborCell(posIdx)->getCons().rhoV
                              + this->negCell->getNeighborCell(posIdx)->getCons().rhoV
                              + this->posCell->getCons().rhoV 
                              + this->negCell->getCons().rhoV 
                              ) * 0.25
                            - ( this->posCell->getNeighborCell(negIdx)->getCons().rhoV
                              + this->negCell->getNeighborCell(negIdx)->getCons().rhoV 
                              + this->posCell->getCons().rhoV 
                              + this->negCell->getCons().rhoV 
                              ) * 0.25
                            ) / dt;

    tangentialGradCons[3] = ( ( this->posCell->getNeighborCell(posIdx)->getCons().rhoE
                              + this->negCell->getNeighborCell(posIdx)->getCons().rhoE
                              + this->posCell->getCons().rhoE 
                              + this->negCell->getCons().rhoE 
                              ) * 0.25
                            - ( this->posCell->getNeighborCell(negIdx)->getCons().rhoE
                              + this->negCell->getNeighborCell(negIdx)->getCons().rhoE 
                              + this->posCell->getCons().rhoE 
                              + this->negCell->getCons().rhoE 
                              ) * 0.25
                            ) / dt;

}

void Interface::computeTimeDerivative(double * prim, double * MomentU, double * MomentV, double * MomentXi,
                                      double* a, double* b, double * timeGrad)
{

    timeGrad[0] = a[0] * MomentU[1] * MomentV[0]
                + a[1] * MomentU[2] * MomentV[0]
                + a[2] * MomentU[1] * MomentV[1]
                + a[3] * ( MomentU[3] * MomentV[0] + MomentU[1] * MomentV[2] + MomentU[1] * MomentV[0] * MomentXi[2] )
                + b[0] * MomentU[0] * MomentV[1]
                + b[1] * MomentU[1] * MomentV[1]
                + b[2] * MomentU[0] * MomentV[2]
                + b[3] * ( MomentU[2] * MomentV[1] + MomentU[0] * MomentV[3] + MomentU[0] * MomentV[1] * MomentXi[2] ) ;

    timeGrad[1] = a[0] * MomentU[2] * MomentV[0]
                + a[1] * MomentU[3] * MomentV[0]
                + a[2] * MomentU[2] * MomentV[1]
                + a[3] * ( MomentU[4] * MomentV[0] + MomentU[2] * MomentV[2] + MomentU[2] * MomentV[0] * MomentXi[2] )
                + b[0] * MomentU[1] * MomentV[1]
                + b[1] * MomentU[2] * MomentV[1]
                + b[2] * MomentU[1] * MomentV[2]
                + b[3] * ( MomentU[3] * MomentV[1] + MomentU[1] * MomentV[3] + MomentU[1] * MomentV[1] * MomentXi[2] );

    timeGrad[2] = a[0] * MomentU[1] * MomentV[1]
                + a[1] * MomentU[2] * MomentV[1]
                + a[2] * MomentU[1] * MomentV[2]
                + a[3] * ( MomentU[3] * MomentV[1] + MomentU[1] * MomentV[3] + MomentU[1] * MomentV[1] * MomentXi[2] )
                + b[0] * MomentU[0] * MomentV[2]
                + b[1] * MomentU[1] * MomentV[2]
                + b[2] * MomentU[0] * MomentV[3]
                + b[3] * ( MomentU[2] * MomentV[2] + MomentU[0] * MomentV[4] + MomentU[0] * MomentV[2] * MomentXi[2] );

    timeGrad[3] = a[0] * 0.50 * ( MomentU[3] * MomentV[0] + MomentU[1] * MomentV[2] + MomentU[1] * MomentV[0] * MomentXi[2] )
                + a[1] * 0.50 * ( MomentU[4] * MomentV[0] + MomentU[2] * MomentV[2] + MomentU[2] * MomentV[0] * MomentXi[2] )
                + a[2] * 0.50 * ( MomentU[3] * MomentV[1] + MomentU[1] * MomentV[3] + MomentU[1] * MomentV[1] * MomentXi[2] )
                + a[4] * 0.25 * ( MomentU[5] + MomentU[1]* ( MomentV[4] + MomentXi[4] )
                                + 2.0 * MomentU[3] * MomentV[2]
                                + 2.0 * MomentU[3] * MomentXi[2]
                                + 2.0 * MomentU[1] * MomentV[2] * MomentXi[2] )
                + b[0] * 0.50 * ( MomentU[2] * MomentV[1] + MomentU[0] * MomentV[3] + MomentU[0] * MomentV[1] * MomentXi[2] )
                + b[1] * 0.50 * ( MomentU[3] * MomentV[1] + MomentU[1] * MomentV[3] + MomentU[1] * MomentV[1] * MomentXi[2] )
                + b[2] * 0.50 * ( MomentU[2] * MomentV[2] + MomentU[0] * MomentV[4] + MomentU[0] * MomentV[2] * MomentXi[2] )
                + b[4] * 0.25 * ( MomentV[5] + MomentV[1] * ( MomentU[4] + MomentXi[4] )
                                + 2.0 * MomentU[2] * MomentV[3]
                                + 2.0 * MomentU[2] * MomentV[1] * MomentXi[2]
                                + 2.0 * MomentV[3] * MomentXi[2] );

    timeGrad[0] /= -prim[0];
    timeGrad[1] /= -prim[1];
    timeGrad[2] /= -prim[2];
    timeGrad[3] /= -prim[3];
}

void Interface::assembleFlux(double * MomentU, double * MomentV, double * MomentXi, double * a, double * b, double * A, double * timeCoefficients)
{
    this->Flux[0] = 
}

void Interface::rotate(double * vector)
{
    double tmp = vector[1];
    vector[1] = vector[2];
    vector[2] = tmp;
}

void Interface::computeMicroSlope(double * prim, double * macroSlope, double * microSlope)
{
    double A, B, C, D, U_2_V_2;

    U_2_V_2 = prim[1] * prim[1] + prim[2] * prim[2];

    A = 2.0*macroSlope[3] - ( U_2_V_2  + (this->fluidParam.K + 2.0) / (2.0*prim[3]) * macroSlope[0] );

    // the product rule of derivations is used here!
    B = macroSlope[1] - prim[1] * macroSlope[0];
    C = macroSlope[2] - prim[2] * macroSlope[0];

    // compute micro slopes of primitive variables from macro slopes of conservatice variables
    microSlope[3] = 2.0 * (4.0 * prim[3]*prim[3])/(this->fluidParam.K + 2.0)
                        * ( A - 2.0*prim[1]*B - 2.0*prim[2]*C );

    microSlope[2] = 2.0 * prim[3] * C - prim[2] * microSlope[3];

    microSlope[1] = 2.0 * prim[3] * B - prim[1] * microSlope[3];

    microSlope[0] = macroSlope[0] - prim[1]*microSlope[1] - prim[2]*microSlope[2] 
                                  - 0.5 * ( U_2_V_2 + (this->fluidParam.K + 2.0) / (2.0*prim[3]) )* microSlope[3];
}

void Interface::computeMoments(double * prim, double * MomentU, double* MomentXi, double * MomentV, int numberMoments)
{
    //==================== U Moments ==========================================
    MomentU[0] = 1.0;
    MomentU[1] = prim[1];
    for (int i = 2; i < numberMoments; i++)
        MomentU[i] = prim[1] * MomentU[i - 1] + 0.5*(i - 1)*MomentU[i - 2] * 2.0 / 3.0;

    //==================== V Moments ==========================================
    MomentV[0] = 1.0;
    MomentV[1] = prim[2];
    for (int i = 2; i < numberMoments; i++)
        MomentV[i] = prim[2] * MomentV[i - 1] + 0.5*(i - 1)*MomentV[i - 2] * 2.0 / 3.0;

    //==================== Xi Moments =========================================
    MomentXi[0] = 1.0;
    MomentXi[1] = 0.0;
    MomentXi[2] = this->fluidParam.K / (2.0 * prim[3]);
    MomentXi[3] = 0.0;
    MomentXi[4] = ( 2.0*this->fluidParam.K + 1.0*this->fluidParam.K*this->fluidParam.K ) / (4.0 * prim[3] * prim[3]);
    MomentXi[5] = 0.0;
}

void Interface::computeMomentUV(double * MomentU, double * MomentV, int alpha, int beta, double * MomentUV)
{
    // compute <u^alpha v^beta>         for rho
    //         <u^alpha+1 v^beta>       for rhoU
    //         <u^alpha v^beta + 1>     for rhoV
    MomentUV[0] = MomentU[alpha]   * MomentV[beta];
    MomentUV[1] = MomentU[alpha+1] * MomentV[beta];
    MomentUV[2] = MomentU[alpha]   * MomentV[beta+1];
    MomentUV[3] = MomentU[alpha+2] * MomentV[beta]
                + MomentU[alpha]   * MomentV[beta+2]
}

void Interface::computeMomentAU(double * a, double * MomentU, double * MomentV, int alpha, int beta, double* MomentAU)
{
    double MomentUV_0[3];
    double MomentUV_1[3];
    double MomentUV_2[3];

    this->computeMomentUV(MomentU, MomentV, alpha  , beta  , MomentUV_0);
    this->computeMomentUV(MomentU, MomentV, alpha+1, beta  , MomentUV_1);
    this->computeMomentUV(MomentU, MomentV, alpha  , beta+1, MomentUV_2);

    for (int i = 0; i < 3; i++)
        MomentAU[i] = a[0] * MomentUV_0[i] 
                    + a[1] * MomentUV_1[i] 
                    + a[2] * MomentUV_2[i]
                    + a[3] * ;
}

