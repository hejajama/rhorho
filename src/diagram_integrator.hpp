#ifndef _DIAGRAM_INTEGRATOR_H
#define _DIAGRAM_INTEGRATOR_H
#include "interpolation.hpp"
#include "functions.hpp"
#include "vector.hpp"
#include "proton.hpp"

#include <gsl/gsl_monte.h>
#include <gsl/gsl_rng.h>

enum Diagram
{
    DIAG_2A,
    DIAG_2B,
    DIAG_3A,
    DIAG_3B,
    DIAG_5A,
    DIAG_5C
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
    ~DiagramIntegrator();
    double IntegrateDiagram(Diagram diag, Vec q1, Vec q2 );
    
    double GetMf() { return mf; }
    Proton& GetProton() { return proton; }
    void SetMCIntPoints(unsigned int n) { MCINTPOINTS = n; }
    void UseInterpolator(bool s) { use_interpolator = s; }
    Interpolator* InitializeInterpolator();
    bool UseInterpolator() { return use_interpolator; }
    double GetX() { return x; }
    F_worker* GetF_worker() { return F; }
    
private:
    double mf;      // Quark mass
    IntegrationMethod intmethod;
    Proton proton;
    gsl_rng *rng;
    F_worker *F;
    unsigned int MCINTPOINTS;
    bool use_interpolator;
    double x; // cutoff
};


#endif
