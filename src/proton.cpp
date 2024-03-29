
#include "proton.hpp"
#include "vector.hpp"
#include <cmath>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>
using namespace std;
inline double SQR(double x){ return x*x; }
/*
 * Evaluate proton wave function
 *
 * Assuems that k1+k2+k3=0, x1+x2+x3=0
 */
double Proton::WaveFunction(Vec k1, Vec k2, double x1, double x2)
{
    if (wf_normalization < 0)
    {
        std::cerr << "Wave function normalization is not defined, exiting..." << std::endl;
        exit(1);
    }
    Vec k3 = (k1+k2)*(-1);
    double x3 = 1.-x1-x2;
    
    if (x3 < 0 or x3 > 1) return 0;
    
    // Invariant mass of the qqq state
    double M2 = (k1.LenSqr()+mq*mq)/x1 + (k2.LenSqr()+mq*mq)/x2 + (k3.LenSqr()+mq*mq)/x3;
    
    if (wave_function ==HarmoinicOscillator)
        return wf_normalization*std::sqrt(x1*x2*x3)*std::exp(-M2/(2.0*beta*beta));
    if (wave_function == Power)
        return wf_normalization*sqrt(x1*x2*x3)*std::pow(1.+M2/(beta*beta), -wf_power);
    
    std::cerr << "Unknown proton wave function!" << std::endl;
    return 0;
}

Proton::Proton()
{
    mq=0.14;
    beta=0.55;
    wf_power=3.5;
    wave_function=HarmoinicOscillator;
    wf_normalization=-1;
}


double Close(double a, double b, double eps=1e-6)
{
    if (std::abs(a-b)<eps) return true;
    else return false;
}
void Proton::ScaleWaveFunctionParameters(double alphas, double x)
{
    if (!Close(beta,0.55))
    {
        cerr << "Proton::ScaleWaveFunctionParameters only works if default value for beta (0.55) is used, but here we have " << beta << endl;
        exit(1);
    }
    if (wave_function != HarmoinicOscillator)
    {
        cerr << "Proton::ScaleWaveFunctionParameters only works  for harmonic oscillator wave funtion" << endl;
        exit(1);
    }

    // Table by Drumitru, Stebel, email 31.1.2024
    if ( Close(alphas, 0.2) and Close(x, 0.3) )
        beta *= 1.014;
    else if ( Close(alphas, 0.25) and Close(x, 0.3) )
        beta *= 1.018;
    else if ( Close(alphas, 0.3) and Close(x, 0.3) )
        beta *= 1.022;
    else if ( Close(alphas, 0.2) and Close(x, 0.1) )
        beta *= 1.057;
    else if ( Close(alphas, 0.25) and Close(x, 0.1))
        beta *= 1.076;
    else if ( Close(alphas, 0.3) and Close(x, 0.1))
        beta *= 1.096;
    else if ( Close(alphas, 0.2) and Close(x, 0.05))
        beta *= 1.1;
    else if (Close(alphas, 0.25) and Close(x,0.05))
        beta *= 1.14;
    else if (Close(alphas, 0.05) and Close(x,0.3))
        beta *= 1.19;
    else
    {
        cerr << "Proton::ScaleWaveFunctionParameters does not support alphas="<< alphas <<", x=" << x << endl;
        exit(1);
    }
}


struct inthelper_norm
{
    Proton* proton;
};



/*
 * Proton wave function normalization
 * https://arxiv.org/pdf/2010.11245.pdf eq 13
 */

double inthelperf_mc_proton_norm(double *vec, size_t dim, void* p)
{
    double k1 = vec[0];
    double k2=vec[1];
    double x1=vec[2];
    double x2=vec[3];
    if (x1+x2 > 1) return 0;
    
    Vec k1v(k1,0);
    Vec k2v(k2*std::cos(vec[4]), k2*std::sin(vec[4]));
    
    double res= 2.0*M_PI*k1*k2*SQR(((inthelper_norm*)p)->proton->WaveFunction(k1v,k2v, x1, x2));
    
    if (isnan(res) or isinf(res))
    {
        cout << res << " at k1=" << k1 << " k2 " << k2 << " x1 " << x1 << " x2 " << x2 << endl;
    }
    
    return 0.5*res/(x1*x2*(1.0-x1-x2)*8*std::pow(2.0*M_PI,6.));
}
double Proton::ComputeWFNormalizationCoefficient()
{
    wf_normalization = 1; // Use this to evaluate wf^2 here
    
    gsl_rng_env_setup ();
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng* rng = gsl_rng_alloc (T);
    
    inthelper_norm helper;
    helper.proton=this;
    gsl_monte_function F;
    F.dim=5;
    F.params = &helper;
    F.f = inthelperf_mc_proton_norm;
    //cout << "Computing normalization, beta=" << beta << ", mq=" << mq << endl;
    
    const double KLIM=20*beta;
    double xlim=0.00001;
    
    double lower[5]={0,0,xlim,xlim,0};
    double upper[5]={KLIM,KLIM,1.-xlim,1.-xlim,2.0*M_PI};
    const int MCINTPOINTS = 1e7;
    
    double result,error;
    gsl_monte_miser_state *s = gsl_monte_miser_alloc(F.dim);
    gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS, rng, s, &result, &error);
    //std::cout << "# Miser result " << result << " err " << error << " relerr " << std::abs(error/result) << std::endl;
    gsl_monte_miser_free(s);
    
    wf_normalization = 1./std::sqrt(result);
    
    return wf_normalization;
    
}

std::string WaveFunctionString(ProtonWaveFunction wf)
{
    if (wf == HarmoinicOscillator)
        return "Harmonic oscillator";
    else if (wf == Power)
        return "Power law";
    else
        return "UNKNOWN!";
    
}
