#ifndef _DIAGRAM_INTEGRATOR_H
#define _DIAGRAM_INTEGRATOR_H
#include "interpolation.hpp"
#include "functions.hpp"
#include "vector.hpp"
#include "proton.hpp"

#include <gsl/gsl_monte.h>
#include <gsl/gsl_rng.h>
#include <string>



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
    DIAG_6E, 
    DIAG_6E_1,
    DIAG_6E_2,
    DIAG_6F,   
    DIAG_6F_1,
    DIAG_6F_2, 
    //DIAG_6G, // With fixed symmetry factors 6g+6g'+6g''=0
    //DIAG_6G_1,
    //DIAG_6G_2, 
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
    
    DIAG_LO,
    
    DIAG_FINITE_SUM,
    DIAG_UV_SUM
    
};
const int FIRST_UV_FINITE = 8; 
const int LAST_UV_FINITE = 36; 
const int FIRST_UV_DIV = 0;
const int LAST_UV_DIV = 7;
const int NUM_OF_DIAGRAMS = 37+3;
const std::string DIAGRAM_STRINGS[NUM_OF_DIAGRAMS] = {"2a", "3a", "3a_2", "3b", "3b_2", "5a", "5c", "5c_1",
        "2b", "3c", "3c_2", "3d", "3d_2",
        "6e","6e_1", // Added when fixing symmetry factors
         "6e_2", 
         "6f", // Added when fixing symmetry factors
         "6f_1",
         "6f_2", // Added when fixing symmetry factors
         "7h", "7i", "7j", "7k", "7l", "7m", "8h_1",
    "8h_2", "8i_1", "8i_2", "8j_1", "8j_2", "8k_1", "8k_2", "8l_1", "8l_2", "8m_1", "8m_2", "LO", "finite_sum","uv_sum"};
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
    DIAG_6E, 
    DIAG_6E_1,
    DIAG_6E_2,
    DIAG_6F,   
    DIAG_6F_1,
    DIAG_6F_2, 
//    DIAG_6G,
//    DIAG_6G_1,
//    DIAG_6G_2, 
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
    
    DIAG_LO,
    
    DIAG_FINITE_SUM,
    DIAG_UV_SUM
};

const int NC=3;
const double CF = (NC*NC-1.)/(2.*NC);

// Color factors
const double COLOR_ADJ = 3.;
const double COLOR_FUND = 0.5;

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
    void SetX(double x_) { x=x_; }
    double GetX() { return x; }
    F_worker* GetF_worker() { return F; }
    
    double DipoleAmplitudeBruteForce(Diagram diag, Vec r, Vec b);
    double MixedSpaceBruteForce(Diagram diag, Vec q12, Vec b);
    
    void SetPerturbativeMass(double mf_) { mf=mf_; }
    
    Diagram DiagramType(std::string str); // Map str to Diagram
    
    bool Add_Q1Q2_exchange(Diagram diag); // Should I add q1<->q2
              
    std::string InfoStr();
    
    void SetSmallX(bool sx) { small_x = x; }
    bool SmallXLimit() { return small_x; }
    
    bool CollinearCutoffUVFinite() { return collinear_cutoff_uv_finite; }
    void SetCollinearCutoffUVFinite(bool s){ collinear_cutoff_uv_finite=s; }
    
private:
    double mf;      // Quark mass
    IntegrationMethod intmethod;
    Proton proton;
    gsl_rng *rng;
    F_worker *F;
    unsigned int MCINTPOINTS;
    bool use_interpolator;
    double x; // cutoff
    bool small_x;
    bool collinear_cutoff_uv_finite; // Use collinear cutoff (m) in UV finite diagrams also
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
