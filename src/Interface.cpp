

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

void Interface::computeMassMomentumFlux(double dt, double tau)
{
    const int NUMBER_OF_MOMENTS = 5;

    double prim[4];
    double normalGradCons[3];
    double tangentialGradCons[3];
    double timeGrad[3];

    double normalA[3];
    double tangentialB[3];
    double timeA[3];

    double MomentU[NUMBER_OF_MOMENTS];
    double MomentV[NUMBER_OF_MOMENTS];

    double MomentUV_1_0[3];

    double normalMomentAU_1_0[3];
    double normalMomentAU_2_0[3];
    double tangentialMomentAU_0_1[3];
    double tangentialMomentAU_1_1[3];
    double timeMomentAU_1_0[3];

    // compute the length of the interface
    double dy = this->posCell->getDx().x * normal.y
              + this->posCell->getDx().y * normal.x;

    // time integration Coefficients in Eq. 10
    double timeCoefficients[3] = { dt, -tau*dt, 0.5*dt*dt - tau*dt };

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

    this->computeMicroSlope(prim, normalGradCons,     normalA    );
    this->computeMicroSlope(prim, tangentialGradCons, tangentialB);

    // temporal micro slopes A = A1 + A2 u + A3 v

    this->computeMomentU(prim, MomentU, MomentV, NUMBER_OF_MOMENTS);
    this->computeMomentAU(normalA,     MomentU, MomentV, 1, 0, normalMomentAU_1_0    );
    this->computeMomentAU(tangentialB, MomentU, MomentV, 0, 1, tangentialMomentAU_0_1);

    for (int i = 0; i < 3; i++)
        timeGrad[i] = - (normalMomentAU_1_0[i] + tangentialMomentAU_0_1[i]);

    this->computeMicroSlope(prim, timeGrad, timeA);

    // compute mass and momentum fluxes

    this->computeMomentUV(MomentU, MomentV, 1, 0, MomentUV_1_0);

    this->computeMomentAU(normalA,     MomentU, MomentV, 2, 0, normalMomentAU_2_0    );
    this->computeMomentAU(tangentialB, MomentU, MomentV, 1, 1, tangentialMomentAU_1_1);

    this->computeMomentAU(timeA, MomentU, MomentV, 1, 0, timeMomentAU_1_0);

    for(int i = 0; i < 3; i++)
        this->MassMomentumFlux[i] = ( timeCoefficients[0] * MomentUV_1_0[i]
                                    + timeCoefficients[1] * ( normalMomentAU_2_0[i] + tangentialMomentAU_1_1[i] )
                                    + timeCoefficients[2] * timeMomentAU_1_0[i]
                                    ) * prim[0] * dy;

        // in case of horizontal interface (G interface), swap velocity fluxes
    if (this->axis == 1)
        this->rotate(this->MassMomentumFlux);

}

void Interface::computeHeatFlux(double dt, double tau)
{
    double prim[4];
    double normalGradTemperature;
    double tangentialGradTemperature;

    double normalA;
    double tangentialB;
    double timeA;

    // compute the length of the interface
    double dy = this->posCell->getDx().x * normal.y
              + this->posCell->getDx().y * normal.x;

    // Coefficients in Eq. 10
    double timeCoefficients[3] = { dt, -tau*dt, 0.5*dt*dt - tau*dt };

    this->interpolatePrim(prim);
    this->differentiateTemperature(normalGradTemperature, tangentialGradTemperature);

    // in case of horizontal interface, swap velocity directions
    if (this->axis == 1)
    {
        this->rotate(prim);
    }

    normalA = normalGradTemperature;
    tangentialB = tangentialGradTemperature;
    timeA = -(normalA*prim[1] + tangentialB*prim[2]);

    this->HeatFlux = ( timeCoefficients[0] * prim[1] * prim[3]
                     + timeCoefficients[1] * ( normalA * ( prim[1] * prim[1] + 1.0/3.0 )
                                             + tangentialB * prim[1] * prim[2] )
                     + timeCoefficients[2] * timeA * prim[1]
                     ) * dy;

}

Cell * Interface::getNeigborCell(Cell * askingCell)
{
    if (posCell == askingCell)
        return negCell;
    else
        return posCell;
}

ConservedVariable Interface::getMassMomentumFlux()
{
    ConservedVariable tmp;
    tmp.rho  = this->MassMomentumFlux[0];
    tmp.rhoU = this->MassMomentumFlux[1];
    tmp.rhoV = this->MassMomentumFlux[2];
    return tmp;
}

double Interface::getHeatFlux()
{
    return this->HeatFlux;
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
    tmp << this->MassMomentumFlux[0] << " " << this->MassMomentumFlux[1] << " " << this->MassMomentumFlux[2];
    tmp << "\n";
    tmp << this->HeatFlux;
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
    prim[3] = 0.5*( this->negCell->getPrim().T
                  + this->posCell->getPrim().T   );
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
                        - this->negCell->getCons().rho  ) / dn / prim[0];

    normalGradCons[1] = ( this->posCell->getCons().rhoU
                        - this->negCell->getCons().rhoU ) / dn / prim[0];

    normalGradCons[2] = ( this->posCell->getCons().rhoV
                        - this->negCell->getCons().rhoV ) / dn / prim[0];

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
                            ) / dt / prim[0];

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
                            ) / dt / prim[0];

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
                            ) / dt / prim[0];

}

void Interface::differentiateTemperature(double& normalGradTemperatur, double& tangentialGradTemperatur )
{
    // ========================================================================
    // normal direction
    // ========================================================================

    // compute the distance between 
    double dn = ( ( this->posCell->getDx().x + this->negCell->getDx().x ) * normal.x
                + ( this->posCell->getDx().y + this->negCell->getDx().y ) * normal.y ) * 0.5;

    normalGradTemperatur = ( this->posCell->getPrim().T
                           - this->negCell->getPrim().T ) / dn;

    // ========================================================================
    // tangential direction
    // ========================================================================

    // get the indieces of the perpendicular interfaces
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

    tangentialGradTemperatur = ( ( this->posCell->getNeighborCell(posIdx)->getPrim().T
                                 + this->negCell->getNeighborCell(posIdx)->getPrim().T
                                 + this->posCell->getPrim().T 
                                 + this->negCell->getPrim().T 
                                 ) * 0.25
                               - ( this->posCell->getNeighborCell(negIdx)->getPrim().T
                                 + this->negCell->getNeighborCell(negIdx)->getPrim().T 
                                 + this->posCell->getPrim().T 
                                 + this->negCell->getPrim().T
                                 ) * 0.25
                               ) / dn;

}

void Interface::rotate(double * vector)
{
    double tmp = vector[1];
    vector[1] = vector[2];
    vector[2] = tmp;
}

void Interface::computeMicroSlope(double * prim, double * macroSlope, double * microSlope)
{
    // compute micro slopes of primitive variables from macro slopes of conservatice variables
    // the product rule of derivations is used here!
    microSlope[2] = 3.0 * ( macroSlope[2] - prim[2] * macroSlope[0] );
    microSlope[1] = 3.0 * ( macroSlope[1] - prim[1] * macroSlope[0] );
    microSlope[0] = macroSlope[0] - prim[1] * microSlope[1] - prim[2] * microSlope[2];
}

void Interface::computeMomentU(double * prim, double * MomentU, double * MomentV, int numberMoments)
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
}

void Interface::computeMomentUV(double * MomentU, double * MomentV, int alpha, int beta, double * MomentUV)
{
    // compute <u^alpha v^beta>         for rho
    //         <u^alpha+1 v^beta>       for rhoU
    //         <u^alpha v^beta + 1>     for rhoV
    MomentUV[0] = MomentU[alpha] * MomentV[beta];
    MomentUV[1] = MomentU[alpha+1] * MomentV[beta];
    MomentUV[2] = MomentU[alpha] * MomentV[beta+1];
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
        MomentAU[i] = a[0] * MomentUV_0[i] + a[1] * MomentUV_1[i] + a[2] * MomentUV_2[i];
}

