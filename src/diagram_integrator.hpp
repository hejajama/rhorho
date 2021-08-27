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
    
    DIAG_LO,
    
    // ODDERON_DIAG_13  vanishes
    ODDERON_DIAG_14,
    ODDERON_DIAG_15,
    ODDERON_DIAG_16,
    ODDERON_DIAG_17,
    ODDERON_DIAG_20,
    ODDERON_DIAG_21,
    ODDERON_DIAG_22,
    ODDERON_DIAG_36,
    ODDERON_DIAG_37,
    
    
    
    ODDERON_DIAG_18,
    ODDERON_DIAG_29,
    ODDERON_DIAG_32,
    ODDERON_DIAG_38,
    ODDERON_DIAG_39,
    ODDERON_DIAG_40,
    
    ODDERON_DIAG_19,
    ODDERON_DIAG_30,
    ODDERON_DIAG_33,
    ODDERON_DIAG_41,
    ODDERON_DIAG_42,
    ODDERON_DIAG_48,
    
    // Rest that do not contribute to ward in q3->0
    ODDERON_DIAG_31,
    ODDERON_DIAG_34,
    ODDERON_DIAG_43,
    ODDERON_DIAG_49,
    
    // UV finite
    ODDERON_DIAG_69,
    ODDERON_DIAG_72,
    ODDERON_DIAG_75,
    
    ODDERON_DIAG_70,
    ODDERON_DIAG_78,
    ODDERON_DIAG_81,
    ODDERON_DIAG_90,
    
    ODDERON_DIAG_71,
    ODDERON_DIAG_79,
    ODDERON_DIAG_82,
    ODDERON_DIAG_91,
    
    ODDERON_DIAG_74,
    ODDERON_DIAG_85,
    ODDERON_DIAG_88,
    ODDERON_DIAG_97,
    
    ODDERON_DIAG_89,
    ODDERON_DIAG_117,
    ODDERON_DIAG_118,
    ODDERON_DIAG_119,
    
    ODDERON_DIAG_76,
    ODDERON_DIAG_93,
    ODDERON_DIAG_99,
    
    
};
const int NUM_OF_DIAGRAMS = 35+9+6+6+4+3+4+4+4+4+3;
const std::string DIAGRAM_STRINGS[] = {"2a", "3a", "3a_2", "3b", "3b_2", "5a", "5c", "5c_1",
        "2b", "3c", "3c_2", "3d", "3d_2", "6e_2", "6f_1",
         "6g_1", "7h", "7i", "7j", "7k", "7l", "7m", "8h_1",
    "8h_2", "8i_1", "8i_2", "8j_1", "8j_2", "8k_1", "8k_2", "8l_1", "8l_2", "8m_1", "8m_2", "LO",
    
    "odderon_14","odderon_15","odderon_16","odderon_17","odderon_20", "odderon_21","odderon_22","odderon_36","odderon_37",
    
    "odderon_18", "odderon_29", "odderon_32", "odderon_38", "odderon_39", "odderon_40",
    
    "odderon_19","odderon_30", "odderon_33", "odderon_41", "odderon_42", "odderon_48",
    
    "odderon_31","odderon_34","odderon_43","odderon_49",
    
    "odderon_69", "odderon_72", "odderon_75",
    "odderon_70", "odderon_78", "odderon_81", "odderon_90",
    "odderon_71", "odderon_79", "odderon_82", "odderon_91",
    "odderon_74", "odderon_85", "odderon_88", "odderon_97",
    "odderon_89", "odderon_117", "odderon_118", "odderon_119",
    "odderon_76", "odderon_93", "odderon_99",
    
};
const Diagram DIAGRAMS[] = {
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
    
    DIAG_LO,
    
    
    // 1st UV diag combination that vanishes in the q3->0 limit (but not q1->0!!!)
    ODDERON_DIAG_14,
    ODDERON_DIAG_15,
    ODDERON_DIAG_16,
    ODDERON_DIAG_17,
    ODDERON_DIAG_20,
    ODDERON_DIAG_21,
    ODDERON_DIAG_22,
    ODDERON_DIAG_36,
    ODDERON_DIAG_37,
    
    
    // 2nd uv divergent
    ODDERON_DIAG_18,
    ODDERON_DIAG_29,
    ODDERON_DIAG_32,
    ODDERON_DIAG_38,
    ODDERON_DIAG_39,
    ODDERON_DIAG_40,
    
    // 3rd uv div
    ODDERON_DIAG_19,
    ODDERON_DIAG_30,
    ODDERON_DIAG_33,
    ODDERON_DIAG_41,
    ODDERON_DIAG_42,
    ODDERON_DIAG_48,
    
    //
    ODDERON_DIAG_31,
    ODDERON_DIAG_34,
    ODDERON_DIAG_43,
    ODDERON_DIAG_49,
    
    //
    ODDERON_DIAG_69,
    ODDERON_DIAG_72,
    ODDERON_DIAG_75,
    
    ODDERON_DIAG_70,
    ODDERON_DIAG_78,
    ODDERON_DIAG_81,
    ODDERON_DIAG_90,
    
    ODDERON_DIAG_71,
    ODDERON_DIAG_79,
    ODDERON_DIAG_82,
    ODDERON_DIAG_91,
    
    ODDERON_DIAG_74,
    ODDERON_DIAG_85,
    ODDERON_DIAG_88,
    ODDERON_DIAG_97,
    
    ODDERON_DIAG_89,
    ODDERON_DIAG_117,
    ODDERON_DIAG_118,
    ODDERON_DIAG_119,
    
    ODDERON_DIAG_76,
    ODDERON_DIAG_93,
    ODDERON_DIAG_99,
    
    
};

const double NC=3;
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
    double IntegrateDiagram(Diagram diag, Vec q1, Vec q2, Vec q3=Vec(0,0) );
    
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
    Vec q3;
    Interpolator *F_B_interpolator;
    Diagram diag;
};

#endif
