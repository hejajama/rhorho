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
    DIAG_3A_2,
    DIAG_3B,
    DIAG_3B_2,
    DIAG_5A,
    DIAG_5C,
    DIAG_5C_1,
    
    // Finite
    DIAG_2B,
    DIAG_3C,
    DIAG_3C_2,
    DIAG_3D,
    DIAG_3D_2,
    //DIAG_6E, // cancel
    //DIAG_6E_1,
    DIAG_6E_2,
    DIAG_6F_1,
    //DIAG_6F_2, // cancels with 6f
    DIAG_6G_1,
    //DIAG_6G_2, // cancels with 6g
    DIAG_7H,
    DIAG_7I,
    DIAG_7J,
    DIAG_7K,
    DIAG_7L,
    DIAG_7M,
    DIAG_8H_1,
    DIAG_8H_2,
    DIAG_8I_1,
    DIAG_8I_2,
    DIAG_8J_1,
    DIAG_8J_2,
    DIAG_8K_1,
    DIAG_8K_2,
    DIAG_8L_1,
    DIAG_8L_2,
    DIAG_8M_1,
    DIAG_8M_2,
    
    DIAG_LO
    
};
const int NUM_OF_DIAGRAMS = 35;
const std::string DIAGRAM_STRINGS[NUM_OF_DIAGRAMS] = {"2a", "3a", "3a_2", "3b", "3b_2", "5a", "5c", "5c_1",
        "2b", "3c", "3c_2", "3d", "3d_2", "6e_2", "6f_1",
         "6g_1", "7h", "7i", "7j", "7k", "7l", "7m", "8h_1",
    "8h_2", "8i_1", "8i_2", "8j_1", "8j_2", "8k_1", "8k_2", "8l_1", "8l_2", "8m_1", "8m_2", "LO"};
const Diagram DIAGRAMS[NUM_OF_DIAGRAMS] = {
    DIAG_2A,
    DIAG_3A,
    DIAG_3A_2,
    DIAG_3B,
    DIAG_3B_2,
    DIAG_5A,
    DIAG_5C,
    DIAG_5C_1,
    
    // Finite
    DIAG_2B,
    DIAG_3C,
    DIAG_3C_2,
    DIAG_3D,
    DIAG_3D_2,
    //DIAG_6E, // cancel
    //DIAG_6E_1,
    DIAG_6E_2,
    DIAG_6F_1,
    //DIAG_6F_2, // cancels with 6f
    DIAG_6G_1,
    //DIAG_6G_2, // cancels with 6g
    DIAG_7H,
    DIAG_7I,
    DIAG_7J,
    DIAG_7K,
    DIAG_7L,
    DIAG_7M,
    DIAG_8H_1,
    DIAG_8H_2,
    DIAG_8I_1,
    DIAG_8I_2,
    DIAG_8J_1,
    DIAG_8J_2,
    DIAG_8K_1,
    DIAG_8K_2,
    DIAG_8L_1,
    DIAG_8L_2,
    DIAG_8M_1,
    DIAG_8M_2,
    
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
