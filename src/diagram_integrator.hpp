#ifndef _DIAGRAM_INTEGRATOR_H
#define _DIAGRAM_INTEGRATOR_H
#include "interpolation.hpp"
#include "vector.hpp"
#include "proton.hpp"

#include <gsl/gsl_monte.h>
#include <gsl/gsl_rng.h>

enum Diagram
{
    DIAG_2A,
};

enum IntegrationMethod
{
    MISER,
    VEGAS
};

class DiagramIntegrator
{
public:
    DiagramIntegrator();
    double IntegrateDiagram(Diagram diag, Vec q1, Vec q2 );
    
    double GetAlpha() { return alpha;}
    double GetMf() { return mf; }
    Proton& GetProton() { return proton; }
    void SetMCIntPoints(unsigned int n) { MCINTPOINTS = n; }
    void UseInterpolator(bool s) { use_interpolator = s; }
    Interpolator* InitializeInterpolator();
    bool UseInterpolator() { return use_interpolator; }
private:
    double alpha;   // Cutoff
    double mf;      // Quark mass
    IntegrationMethod intmethod;
    Proton proton;
    gsl_rng *rng;
    unsigned int MCINTPOINTS;
    bool use_interpolator;
};


#endif
