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
    DIAG_3A,
    DIAG_3B,
    DIAG_5A,
    DIAG_5C,
    
    // Finite
    DIAG_2B,
    DIAG_3C,
    DIAG_3D,
    DIAG_6E_1,
    DIAG_6E_2,
    DIAG_6F_1,
    DIAG_6F_2,
    DIAG_6G_1,
    DIAG_6G_2,
    DIAG_7H,
    DIAG_7J,
    DIAG_7K,
    DIAG_7L,
    DIAG_8H_1,
    DIAG_8H_2,
    DIAG_8J_1,
    DIAG_8J_2,
    DIAG_8L_1,
    DIAG_8L_2,
    
    DIAG_LO
    
};
const int NUM_OF_DIAGRAMS = 25;
const std::string DIAGRAM_STRINGS[NUM_OF_DIAGRAMS] = {"2a", "3a", "3b", "5a", "5c",
        "2b", "3c", "3d", "6e_1", "6e_2", "6f_1",
        "6f_2", "6g_1", "6g_2", "7h", "7j", "7k", "7l", "8h_1",
    "8h_2", "8j_1", "8j_2", "8l_1", "8l_2", "LO"};
const Diagram DIAGRAMS[NUM_OF_DIAGRAMS] = {
    DIAG_2A,
    DIAG_3A,
    DIAG_3B,
    DIAG_5A,
    DIAG_5C,
    
    // Finite
    DIAG_2B,
    DIAG_3C,
    DIAG_3D,
    DIAG_6E_1,
    DIAG_6E_2,
    DIAG_6F_1,
    DIAG_6F_2,
    DIAG_6G_1,
    DIAG_6G_2,
    DIAG_7H,
    DIAG_7J,
    DIAG_7K,
    DIAG_7L,
    DIAG_8H_1,
    DIAG_8H_2,
    DIAG_8J_1,
    DIAG_8J_2,
    DIAG_8L_1,
    DIAG_8L_2,
    
    DIAG_LO
};

const int NC=3;
const double CF = (NC*NC-1.)/(2.*NC);

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
    
    Diagram DiagramType(std::string str); // Map str to Diagram
    
    bool Add_Q1Q2_exchange(Diagram diag); // Should I add q1<->q2
    
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


struct inthelper_diagint
{
    DiagramIntegrator *integrator;
    Vec q1;
    Vec q2;
    Interpolator *F_B_interpolator;
    Diagram diag;
};

#endif
