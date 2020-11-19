#include "functions.hpp"
#include <gsl/gsl_integration.h>
#include <cmath>

inline double SQR(double x) { return x*x; }

/* B0 from PV reduction
 *
 * Delta = m^2(1-z_1)^2 [GeV^2]
 * hsqr [GeV^2]
 *
 */

double B0(double hsqr, double Delta)
{
    /*if (Delta < 1e-20)
    {
        //std::cerr << "B0 only works at finite Delta! got Delta=" << Delta << std::endl;
        return 0;
    }*/
    // Analtyical results
    double res =  -(2*std::log((2.*Delta)/(2.* Delta + hsqr + std::sqrt(hsqr* (4.*Delta + hsqr))))/std::sqrt(hsqr*(4.*Delta + hsqr)));

    if (isnan(res) or isinf(res))
    {
        std::cout << res << " , h^2=" << hsqr << ", Delta=" << Delta << std::endl;
    }
    return res;

}

F_worker::F_worker(const int int_divisions_, const double intaccuracy)
{
    int_accuracy=intaccuracy;
    int_divisions=int_divisions_;
    ws = gsl_integration_workspace_alloc(int_divisions);
}

F_worker::~F_worker()
{
    gsl_integration_workspace_free(ws);
    
}

/*
 * Integral of B0
 * Risto, Adrian (A22)
 * Note that special case l_1=0 simplifies!
 */

struct inthelper_F{ Vec l; Vec l1; double m2; };
double inthelperf_F(double z1, void* p)
{
    inthelper_F* par = (inthelper_F*)p;
    double hsqr = par->l1.LenSqr() + par->l.LenSqr() * SQR(1.-z1)
        - 2.*(1.-z1)* (par->l1*par->l);
    
    if (hsqr < 1e-7)
    {
        // We always have h^2 B0, and in the limit h^2 -> 0 we easily see
        // that h^2 B_0 -> 0, so we can just cut out this part of the phase space
        return 0;
    }
    
    double Delta = SQR(1.-z1)*par->m2;
    return 1.0/z1 * (1.0 + SQR(1.-z1)) * hsqr/2. * B0(hsqr, Delta);
}
double F_worker::F_int_B0(Vec l, Vec l1, double alpha, double m2)
{
    
    // use analytical result if l1=0
    if (l1.LenSqr() < 1e-7)
    {
        double l4l2m2 =SQR(l.LenSqr()) + 4.*l.LenSqr()*m2; //l^4 + 4l^2m^2
        double res = -l.LenSqr() * std::log(2.0*m2 / (l.LenSqr() + 2.*m2+std::sqrt( l4l2m2 ) ) )
            * (3. - 4.*alpha + SQR(alpha) + 4.*std::log(alpha))
        / (16.*SQR(M_PI)*std::sqrt(l4l2m2) );
        
        return res;
    }
    
    
    gsl_function fun;
    fun.function=inthelperf_F;
    inthelper_F helper; helper.l=l; helper.l1=l1; helper.m2=m2;
    fun.params=&helper;
    double result, abserr;
    
    int status = gsl_integration_qag(&fun, alpha, 1., 0, int_accuracy,
        int_divisions, GSL_INTEG_GAUSS41, ws, &result, &abserr);


    
    if (status)
    {
        std::cerr << "Warning: inaccurate F integral = " << result << " +/- " << abserr << std::endl;
    }
    
    if (isnan(result) or isinf(result))
    {
        std::cerr << "Result " << result << " l=" << l <<", l1=" << l1 << std::endl;
    }
    
    return -1./(8.*SQR(M_PI))*result;
     
}


