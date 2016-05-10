
#ifndef TYPES_H
#define TYPES_H

struct PrimaryVariable
{
	double rho;
	double U;
	double V;
	double T;
};

struct ConservedVariable
{
	double rho;
	double rhoU;
	double rhoV;
};

struct float2
{
	double x;
	double y;
};

struct Parameters
{
	unsigned int numberOfIterations;
    unsigned int outputInterval;

	double nu;
	double k;
    double tauMassMomentum;
    double tauHeat;

	double G0;
	double beta;
	double TAve;

	double CFL;

    bool verbose;
};

#endif