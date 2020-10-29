
#include "diagram_integrator.hpp"
#include "functions.hpp"
#include "interpolation.hpp"
#include "proton.hpp"
#include <vector>
#include <cmath>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte_vegas.h>

using namespace std;

struct inthelper_diagint
{
    DiagramIntegrator *integrator;
    Vec q1;
    Vec q2;
    Interpolator *F_B_interpolator;
};

/*
 * Vec components are [k1x,k1y,k2x,k2y,x1,x2]
 */
double inthelperf_mc_diag2a(double *vec, size_t dim, void* p)
{
    if (dim != 6) exit(1);
    inthelper_diagint *par = (inthelper_diagint*)p;
    Vec k1(vec[0],vec[1]);
    Vec k2(vec[2], vec[3]);
    Vec q1 = par->q1;
    Vec q2 = par->q2;
    
    double x1=vec[4];
    double x2=vec[5];
    //double x3 = 1.-x1-x2;
    
    //Vec k3 = (k1+k2)*(-1);
    
    double wf1 = par->integrator->GetProton().WaveFunction( k1, k2, x1,  x2);
    
    // Momentum vectors to the wf after probing the proton
    // k1' = k1 - (1-x1)*(q1+q2)
    Vec k12 = q1+q2; k12*=(1.-x1); k12 = k1 - k12;
    // k2' = k2 + x2(q2+q2)
    Vec k22 = q1+q2; k22*=x2; k22 = k2 + k22;
    // k3' = k3 + x3(q1+q2); fixed by the fact that no pT in the initial state (in the conjugate amplitude)
    
    //Vec k32 = q1+q2; k32*=x3; k32 = k3 + k32;
    
    double wf2 = par->integrator->GetProton().WaveFunction( k12, k22, x1, x2);
    
    double alpha = par->integrator->GetAlpha(); double mf =  par->integrator->GetMf();
    
    if (alpha < 1e-8 or mf < 1e-8)
    {
        cerr << "Invalid alpha=" << alpha << ", mf=" << mf << " GeV!" << endl;
        exit(1);
    }
    
    Vec l = q1+q2; Vec l1(0,0,0);
    
    double fintb = 0;
    if (par->integrator->UseInterpolator() == true)
        fintb = par->F_B_interpolator->Evaluate(l.Len());
    else
        fintb =F_int_B0(l, l1, alpha, mf*mf);
    double result = wf1*wf2*fintb;
    
    if (isinf(result) or isnan(result))
    {
        cerr << "Result "<< result << " k1=" << k1 <<", k2=" << k2 << " wf1 " << wf1 << " wf2 " << wf2 << endl;
    }
   
    return 2.0*std::pow(M_PI,3.)*result; // A21 gives 2pi^3
    
}
double DiagramIntegrator::IntegrateDiagram(Diagram diag, Vec q1, Vec q2 )
{
    inthelper_diagint helper;
    helper.q1=q1; helper.q2=q2; helper.integrator=this;
    gsl_monte_function F;
    F.dim=6;
    F.params = &helper;
    const double KLIM=9;
    double xlim=0.001;
    
    double lower[6]={-KLIM,-KLIM,-KLIM,-KLIM,xlim,xlim};
    double upper[6]={KLIM,KLIM,KLIM,KLIM,1.-xlim,1.-xlim};
    
    switch (diag) {
        case DIAG_2A:
            F.f=inthelperf_mc_diag2a;
            break;
            
        default:
            std::cerr << "Unknown diagram!" << std::endl;
            return 0;
    }
    
    Interpolator *F_b_interp;
    if (use_interpolator)
    {
        F_b_interp = InitializeInterpolator();
        helper.F_B_interpolator = F_b_interp;
    }
    
    double result,error;
    if (intmethod == MISER)
    {
        gsl_monte_miser_state *s = gsl_monte_miser_alloc(F.dim);
        gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS, rng, s, &result, &error);
        cout << "# Miser result " << result << " err " << error << " relerr " << std::abs(error/result) << endl;
        gsl_monte_miser_free(s);
    }
    else if (intmethod == VEGAS)
    {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(F.dim);
        gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/2, rng, s, &result, &error);
        cout << "# vegas warmup " << result << " +/- " << error << endl;
        int iter=0;
        do
        {
            gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS, rng, s, &result, &error);
            cout << "# Vegas interation " << result << " +/- " << error << " chisqr " << gsl_monte_vegas_chisq(s) << endl;
            iter++;
        } while ((fabs( gsl_monte_vegas_chisq(s) - 1.0) > 0.4 or iter < 2) and iter < 5);
        gsl_monte_vegas_free(s);
    }
    else
    {
        result=0;error=0;
    }
    
    
    if (use_interpolator)
        delete F_b_interp;
    
    return result;
    
}

Interpolator* DiagramIntegrator::InitializeInterpolator()
{
    double minq=0.01;
    double maxq=10;
    int npoints=100;
    
    std::vector<double> qvals; std::vector<double> F;
    for (double q=minq; q<=maxq; q+=(maxq-minq)/npoints)
    {
        qvals.push_back(q);
        Vec l(0,q);
        F.push_back(F_int_B0(l, Vec(0,0), alpha, mf*mf));
    }
    Interpolator *interp = new Interpolator(qvals,F);
    interp->SetFreeze(true); interp->SetUnderflow(0); interp->SetOverflow(0);
    return interp;
    
}


DiagramIntegrator::DiagramIntegrator()
{
    mf=0.1;
    alpha=0.01;
    intmethod = VEGAS;
    
    proton.SetBeta(0.55);
    proton.SetM(0.26);
    cout <<"# Initializing proton..." << endl;
    proton.ComputeWFNormalizationCoefficient();
    cout << "#... done" << endl;
    gsl_rng_env_setup ();

    const gsl_rng_type *T = gsl_rng_default;
    rng = gsl_rng_alloc (T);
}


