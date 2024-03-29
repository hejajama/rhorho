
#include "diagram_integrator.hpp"
#include "functions.hpp"
#include "interpolation.hpp"
#include "proton.hpp"
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte_vegas.h>

using namespace std;

double inthelperf_mc_finitesum(double *vec, size_t dim, void* p);
double inthelperf_mc_uvsum(double* vec, size_t dim, void* p);

/*
 * LO diagram
 * [k1x,k1y,k2x,k2y,x1,x2]
 */
double inthelperf_mc_lo(double *vec, size_t dim, void* p)
{
    if (dim != 6) exit(1);
    inthelper_diagint *par = (inthelper_diagint*)p;
    Vec k1(vec[0]*std::cos(vec[1]),vec[0]*std::sin(vec[1]) );
    Vec k2(vec[2]*std::cos(vec[3]), vec[2]*std::sin(vec[3]) );
    Vec q1 = par->q1;
    Vec q2 = par->q2;
    
    // Assume P=0
    Vec p1 = k1;
    Vec p2 = k2;
    
    double x1=vec[4];
    double x2=vec[5];
    
    double x3 = 1.-x1-x2;
    double x = par->integrator->GetX();
    if (x3 >= 1 or x3 < x) return 0;
    if (x1+x2 >=1) return 0;
    
    double wf1 = par->integrator->GetProton().WaveFunction( p1, p2, x1,  x2);
    
    
    // Risto (77)
    double wf2 =par->integrator->GetProton().WaveFunction(k1 - (q1+q2)*(1.-x1), k2+(q1+q2)*x2, x1, x2)
    - 0.5*par->integrator->GetProton().WaveFunction(k1+(q1+q2)*x1-q1, k2+(q1+q2)*x2-q2, x1, x2)
    - 0.5*par->integrator->GetProton().WaveFunction(k1+(q1+q2)*x1-q2, k2+(q1+q2)*x2-q1, x1, x2);
    
    if (isnan(wf1) or isnan(wf2))
        cerr <<"Note: WF NaN at x1 " << x1 << " x2 "<<  x2 << " k1 " << k1 << " k2 " << k2 << " q1 " << q1 << " q2 " << q2 << endl;
    
    double res = wf1*wf2 * 1./6.*3*COLOR_FUND; // 3 is the symmetry factor and 1/6 the normalization in Risto (77)
    
    // 16pi^3 because NLO diagrams do not include 1/(16pi^3) prefactor
    // Python analysis notebook divides by 1/16pi^3
    return 16.*std::pow(M_PI,3.) * res * vec[0] * vec[2] / (8.0*x1*x2*(1.-x1-x2)*std::pow(2.0*M_PI,6.0));
    

}

/*
 * UV finite diagrams
 * Vector components are
 * [k1, theta_1, k2, theta_2, x1, x2, xg, kg, theta_g]
 */
double inthelperf_mc_diag2b(double *vec, size_t dim, void* p)
{
    if (dim != 9) exit(1);
    inthelper_diagint *par = (inthelper_diagint*)p;
    
    // If the integrand is diag sum, use specific integrand funtions
    if (par->diag == DIAG_FINITE_SUM)
    {
         return inthelperf_mc_finitesum(vec, dim, p);
    }
    
    std::cerr << "NOTE: No color factors included when computing individual diagrams!" << endl;
    
    Vec k1(vec[0]*std::cos(vec[1]), vec[0]*std::sin(vec[1]));
    Vec k2(vec[2]*std::cos(vec[3]), vec[2]*std::sin(vec[3]));
    Vec kg(vec[7]*std::cos(vec[8]), vec[7]*std::sin(vec[8]));
    Vec q1 = par->q1;
    Vec q2 = par->q2;
    Vec q = q1+q2;
    Vec K = q1*(-1) - q2;
    
    double x1=vec[4];
    double x2=vec[5];
    
    double x3 = 1.-x1-x2;
    double x = par->integrator->GetX();
    if (x3 >= 1 or x3 < x) return 0;
    if (x1+x2 >=1) return 0;
    
    double xg = vec[6];
    
    // This is bit of a hack, but for now the small xg limit
    // is computed by performing exactly the same integrals, so we don't
    // do one simple integral analycially in order to keep the code structure
    // the same.
    double inv_xg = 1.0/xg;
    
    // Do not allow exactly the upper limit.
    if (xg > std::min(x1,1.-x2)-1e-4) return 0;
    
    if (par->integrator->SmallXLimit()== true)
    {
        xg = 0;
    }
    
    double z1,z2;
    z1 = xg/x1; z2 = xg /( x2+xg );
    
    
    /// We work in the frame where P=0
    Vec p1 = k1; Vec p2 = k2;
    
    Vec ktilde_1; Vec ktilde_2;
    Vec A,B;
    double f_xg=std::sqrt(x1*x2/((x1-xg)*(x2+xg))) * (1. - (z1+z2)/2. + z1*z2/6.);
    double norm=1; // normalization * symmetry factor
    
    switch (par->diag) {
        case DIAG_2B:
            ktilde_1 = k1 + q*x1  - kg + K*xg;
            ktilde_2 = k2 - q*(1.-x2) + kg - K*xg;
            A = p1*z1 - kg;
            B = p2*z2 - (kg-q1-q2)*(1.-z2);
            norm=-1./6. * 6;
            break;
        case DIAG_3C:
            ktilde_1 =k1 + q*x1 - q1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q2 + kg - K*xg;
            A = p1*z1 - kg;
            B = p2*z2 - (kg-q2)*(1.-z2);
            norm=1./12.*6;
            break;
        case DIAG_3C_2: // q1 <->q2 swapped
            ktilde_1 =k1 + q*x1 - q2 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q1 + kg - K*xg;
            A = p1*z1 - kg;
            B = p2*z2 - (kg-q1)*(1.-z2);
            norm=1./12.*6;
            break;
        case DIAG_3D:
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 - q*(1.-x2)+kg-K*xg;
            A = p1*z1 - kg;
            B = (p2-q1)*z2 - (kg-q2)*(1.-z2);
            norm=1./12.*6.;
            break;
        case DIAG_3D_2:
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 - q*(1.-x2)+kg-K*xg;
            A = p1*z1 - kg;
            B = (p2-q2)*z2 - (kg-q1)*(1.-z2);
            norm=1./12.*6.;
            break;
        case DIAG_6E: // (82): now included when symmetry factors fixed
            ktilde_1 = k1 - q*(1.-x1) - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = -CF/3.*6;
            break;
        case DIAG_6E_1: // (85), previously  6E + 6E_1 = 0, but now not true with fixed symmetry factors
            ktilde_1 = k1 - q*(1.-x1) - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = CF/3. * 3;
            break;
        case DIAG_6E_2: // Risto (86)
            ktilde_1 = k1 - q*(1.-x1) - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            A = (p1-q)*z1 -kg;
            B = p2*z2 - kg*(1.-z2);
            norm = CF/3. * 3;   // Symmetry factor fixed
            break;
        case DIAG_6F: // Cancelled with with 6f'' originally, not anymore with fixed symmetry factors
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 - q*(1.-x2) + kg - K*xg;
            A = p1*z1 - kg;
            B = (p2-q)*z2 - kg*(1.-z2);
            norm = -CF/3.*6;
            break;
        case DIAG_6F_1: // (87)
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 - q*(1.-x2) + kg - K*xg;
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = CF/3. * 3;
            break;
        case DIAG_6F_2: // Originally cancelled with 6f, not with fixed symmetry factors, (88)
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 - q*(1.-x2) + kg - K*xg;
            A = p1*z1-kg;
            B = (p2-q)*z2 - kg*(1.-z2);
            norm = CF/3. * 3;
            break;
        
/* 6g+6g'+6g''=0 with fixed symmetry factors
         case DIAG_6G: // Included with fixed symmetry factors, originally cancelled with 6G''
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg; 
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = -1./3. * CF * 6;
            break;
        case DIAG_6G_1: // (89)
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = CF*1./3. * 3; // Symmetry factor fixed
            break;
        case DIAG_6G_2:  // Originally cancelled with 6G, not with fixed symmetry factors
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = 1./3. * CF * 3; // Symmetry factor fixed
            break;
            */
        case DIAG_7H:
            ktilde_1 = k1 + q*x1 - q2 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q1 + kg - K*xg;
            A = p1*z1 - kg;
            B = (p2-q1)*z2 - kg*(1.-z2);
            norm = 1./3.*(0.5-CF) * 6;
            break;
        case DIAG_7I: // 7H, q1 <-> q2 swap
            ktilde_1 = k1 + q*x1 - q1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q2 + kg - K*xg;
            A = p1*z1 - kg;
            B = (p2-q2)*z2 - kg*(1.-z2);
            norm = 1./3.*(0.5-CF) * 6;
            break;
        case DIAG_7J:
            ktilde_1 = k1 + q*x1 - q2 - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = 1./3. * (CF-1./2. - 1./6.)*6;
            break;
        case DIAG_7K: // 7K, q1 <-> q2
            ktilde_1 = k1 + q*x1 - q1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = 1./3. * (CF-1./2. - 1./6.)*6;
            break;
        case DIAG_7L:
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q2 + kg - K*xg;
            A = p1*z1 - kg;
            B = (p2-q2)*z2 - kg*(1.-z2);
            norm = 1./3. * (CF-1./2.-1./6.)*6;
            break;
        case DIAG_7M: // 7L, q1 <-> q2
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q1 + kg - K*xg;
            A = p1*z1 - kg;
            B = (p2-q1)*z2 - kg*(1.-z2);
            norm = 1./3. * (CF-1./2.-1./6.)*6;
            break;
        case DIAG_8H_1:
            ktilde_1 = k1 + q*x1 - kg - q2 + K*xg;
            ktilde_2 = k2 + q*x2 + kg - q1 - K*xg;
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = -2./9.*3;
            break;
        case DIAG_8H_2:
            ktilde_1 = k1 + q*x1 - kg - q2 + K*xg;
            ktilde_2 = k2 + q*x2 + kg - q1 - K*xg;
            A = (p1-q2)*z1-kg;
            B = (p2-q1)*z2 - kg*(1.-z2);
            norm = -2./9. * 3;
            break;
        case DIAG_8I_1: // 8H_1, q1 <-> q2
            ktilde_1 = k1 + q*x1 - kg - q1 + K*xg;
            ktilde_2 = k2 + q*x2 + kg - q2 - K*xg;
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = -2./9. * 3;
            break;
        case DIAG_8I_2: // 8H_2, q1 <-> q2
            ktilde_1 = k1 + q*x1 - kg - q1 + K*xg;
            ktilde_2 = k2 + q*x2 + kg - q2 - K*xg;
            A = (p1-q1)*z1-kg;
            B = (p2-q2)*z2 - kg*(1.-z2);
            norm = -2./9. * 3;
            break;
        case DIAG_8J_1:
            ktilde_1 = k1 + q*x1 - kg - q2 + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = -1./3. * (CF-2./3.)*3;
            break;
        case DIAG_8J_2:
            ktilde_1 = k1 + q*x1 - kg - q2 + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            A = (p1-q2)*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = -1./3. * (CF-2./3.)*3;
            break;
        case DIAG_8K_1: // 8J_1, q1 <-> q2
            ktilde_1 = k1 + q*x1 - kg - q1 + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = -1./3. * (CF-2./3.)*3;
            break;
        case DIAG_8K_2: // 8J_2, q1 <-> q2
            ktilde_1 = k1 + q*x1 - kg - q1 + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            A = (p1-q1)*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = -1./3. * (CF-2./3.)*3;
            break;
        case DIAG_8L_1:
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - q2 - K*xg;
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = -1./3. * (CF-2./3)*3;
            break;
        case DIAG_8L_2:
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - q2 - K*xg;
            A = p1*z1 - kg;
            B = (p2-q2)*z2 - kg*(1.-z2);
            norm = -1./3. * (CF-2./3)*3;
            break;
        case DIAG_8M_1: // 8L_1, q1<->q2
            ktilde_1 =  k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - q1 - K*xg;
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = -1./3. * (CF-2./3)*3;
            break;
        case DIAG_8M_2: // 8L_2, q1<->q2
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - q1 - K*xg;
            A = p1*z1 - kg;
            B = (p2-q1)*z2 - kg*(1.-z2);
            norm = -1./3. * (CF-2./3)*3;
            break;
        default:
            cerr << "Unknown diagram in inthelperf_mc_diag2b: " << par->diag << endl;
            exit(1);
            break;
    }
    
    if (A.LenSqr() < 1e-15 or B.LenSqr() < 1e-15)
        return 0;
    
   
    double wf1 =par->integrator->GetProton().WaveFunction(k1, k2, x1, x2);
    double wf2 = par->integrator->GetProton().WaveFunction(ktilde_1,ktilde_2,x1-xg, x2+xg);
    
    double dotprod = 0;
    double mf = par->integrator->GetMf();
    if (par->integrator->CollinearCutoffUVFinite())
        dotprod = (A*B) / ( (A.LenSqr()+mf*mf)*(B.LenSqr()+mf*mf));
    else
        dotprod =(A*B)/(A.LenSqr()*B.LenSqr());
    
    double res = norm*wf1*wf2*f_xg*dotprod;
    
    res *= inv_xg; // same as res /= xg;
    
    // Jacobian
    res *= vec[0] * vec[2] * vec[7];
    res /= 8.0*x1*x2*(1.-x1-x2)*std::pow(2.0*M_PI,6.0);
    
    return res;
    
}



/*
 * Vec components are [k1, theta_1, k2; theta_2, x1, x2]
 */
double inthelperf_mc_diag2a(double *vec, size_t dim, void* p)
{
    if (dim != 6) exit(1);
    inthelper_diagint *par = (inthelper_diagint*)p;
    
    // If the integrand is diag sum, use specific integrand funtions
    if (par->diag == DIAG_UV_SUM)
    {
         return inthelperf_mc_uvsum(vec, dim, p);
    }
    
    std::cerr << "NOTE: No color factors included when computing individual diagrams!" << endl;
    Vec k1(vec[0]*std::cos(vec[1]), vec[0]*std::sin(vec[1]));
    Vec k2(vec[2]*std::cos(vec[3]), vec[2]*std::sin(vec[3]));
    Vec q1 = par->q1;
    Vec q2 = par->q2;
    
    double x1=vec[4];
    double x2=vec[5];
    
    double x3 = 1.-x1-x2;
    double x = par->integrator->GetX();
    if (x3 >= 1 or x3 < x) return 0;
    if (x1+x2 >=1) return 0;
    
    double wf1 = par->integrator->GetProton().WaveFunction( k1, k2, x1,  x2);
    
    Vec k12; Vec k22;
    Vec l; Vec l1;
    double norm=1; // Normalization factor * symmetry factor,not including g^4 / 16pi^3
    
    double alpha =  par->integrator->GetX() / x1;
    // Default, changed below if necessary
    switch(par->diag)
    {
        case DIAG_2A:
            l=q1+q2;
            l1=Vec(0,0);
            k12 =k1 - (q1+q2)*(1.-x1);
            k22=k2 + (q1+q2)*x2;
            norm = 2./3. * 3.;
            break;
        case DIAG_3A:
            l=q1+q2;
            l1=q1;
            k12=k1 - (q1+q2)*(1.-x1);
            k22=k2 + (q1+q2)*x2;
            norm = -1./3.*3;
            break;
        case DIAG_3A_2:
            l=q1+q2;
            l1=q2;
            k12=k1 - (q1+q2)*(1.-x1);
            k22=k2 + (q1+q2)*x2;
            norm = -1./3.*3;
            break;
        case DIAG_3B:
            l=q2;
            l1=Vec(0,0);
            k12 = k1 + (q1+q2)*x1-q2;
            k22 = k2 + (q1+q2)*x2-q1;
            norm = -1./6. * 6;
            break;
        case DIAG_3B_2:
            l=q1;
            l1=Vec(0,0);
            k12 = k1 + (q1+q2)*x1-q1;
            k22 = k2 + (q1+q2)*x2-q2;
            norm = -1./6. * 6;
            alpha = par->integrator->GetX() / x2;
            break;
        case DIAG_5A:
            l=q1+q2;
            l1=q1+q2;
            k12 = k1 - (q1+q2)*(1.-x1);
            k22 = k2 + (q1+q2)*x2;
            norm = 4.*CF/3.*3;
            break;
        case DIAG_5C:
            l=q2;
            l1=q2;
            k12 = k1+(q1+q2)*x1-q2;
            k22 = k2+(q1+q2)*x2-q1;
            norm = 2./(3.*6.)*6;
            break;
        case DIAG_5C_1:
            l=q1;
            l1=q1;
            k12 = k1+(q1+q2)*x1-q1;
            k22 = k2+(q1+q2)*x2-q2;
            norm = 2./(3.*6.)*6;
            alpha = par->integrator->GetX() / x2;
            break;
        default:
            std::cerr << "Unknown diagram " << par->diag << " encountered!" << std::endl;
            exit(1);
    }
    
    if (l1.LenSqr() < 1e-20 and l.LenSqr() < 1e-20)
        return 0;
     

    double wf2 = par->integrator->GetProton().WaveFunction( k12, k22, x1, x2);
    
    
    double mf =  par->integrator->GetMf();
    
    if (alpha < 1e-8 or mf < 1e-8)
    {
        cerr << "Invalid alpha=" << alpha << ", mf=" << mf << " GeV!" << endl;
        exit(1);
    }
    
    
    double fintb = 0;
    
    if (par->integrator->SmallXLimit())
    {
        // A22 in the small-x limit
        double hsqr = l1.LenSqr() + l.LenSqr() - 2.0*(l*l1);
        if (hsqr < 1e-7)
            return 0; // h^2 B0 -> 0
        
        fintb = -1.0/(4.0*M_PI*M_PI) * std::log(alpha) * std::log(mf*mf*alpha/hsqr);
    }
    else
    {
        if (par->integrator->UseInterpolator() == true)
            fintb = par->F_B_interpolator->Evaluate(l.Len());
        else
            fintb = par->integrator->GetF_worker()->F_int_B0(l, l1, alpha, mf*mf);
    }
    double result = norm*wf1*wf2*fintb;
        
    if (isinf(result) or isnan(result))
    {
        cerr << "Result "<< result << " k1=" << k1 <<", k2=" << k2 << " wf1 " << wf1 << " wf2 " << wf2 << endl;
    }
    
    // Jacobian
    result *= vec[0]*vec[2];
    result /= x1*x2*(1.-x1-x2)*8*std::pow(2.0*M_PI,6.);
   
    return 2.0*std::pow(M_PI,3.)*result; // A21 gives 2pi^3
    
}
double DiagramIntegrator::IntegrateDiagram(Diagram diag, Vec q1, Vec q2 )
{
    inthelper_diagint helper;
    helper.q1=q1; helper.q2=q2; helper.integrator=this;
    helper.diag = diag;
    gsl_monte_function F;
    
    F.params = &helper;
    const double KLIM=15;
    double xup=1.-0.0001;
    double xlow = x;
        
    double *lower;
    double* upper;
    
    switch (diag) {
        case DIAG_2A:
        case DIAG_3A:
        case DIAG_3A_2:
        case DIAG_3B:
        case DIAG_3B_2:
        case DIAG_5A:
        case DIAG_5C:
        case DIAG_5C_1:
        case DIAG_LO:
        case DIAG_UV_SUM:
            F.dim=6;
            lower = new double[F.dim];
            upper = new double [F.dim];
            lower[0]=lower[1]=lower[2]=lower[3]=0;
            lower[4]=lower[5]=xlow;
            upper[0]=upper[2]=KLIM;
            upper[1]=upper[3]=2.0*M_PI;
            upper[4]=upper[5]=xup;
            F.f=inthelperf_mc_diag2a;
            break;
        default: // Other diagrams are of this type
            F.dim=9;
            lower = new double[F.dim];
            upper = new double [F.dim];
            lower[0]=lower[1]=lower[2]=lower[3]=lower[7]=lower[8]=0;
            lower[4]=lower[5]=xlow; lower[6]=x;
            
            upper[0]=upper[2]=upper[7]=KLIM;
            upper[1]=upper[3]=upper[8]=2.0*M_PI;
            upper[4]=upper[5]=upper[6]=xup;
            F.f=inthelperf_mc_diag2b;
            break;
    }
    
    if (diag == DIAG_LO)
        F.f = inthelperf_mc_lo;
    
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
        //cout << "# vegas warmup " << result << " +/- " << error << endl;
        int iter=0;
        do
        {
            gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS, rng, s, &result, &error);
            //cout << "# Vegas interation " << result << " +/- " << error << " chisqr " << gsl_monte_vegas_chisq(s) << endl;
            iter++;
        } while ((fabs( gsl_monte_vegas_chisq(s) - 1.0) > 0.4 or iter < 2) and iter < 6);
        gsl_monte_vegas_free(s);
    }
    else
    {
        result=0;error=0;
    }
    
    
    
    delete lower;
    delete upper;
    
    if (use_interpolator)
        delete F_b_interp;
    
    return result;
    
}



///////////////
// Mixed space brute force
////
/////// Dipole ampiltude

struct mixed_space_helper
{
    DiagramIntegrator* integrator;
    Vec q12;
    Vec b;
    Diagram diag;
};

double inthelperf_mc_mixedspace(double *vec, size_t dim, void* p)
{
    /*if (dim != 4) exit(1);*/
    mixed_space_helper *par = (mixed_space_helper*)p;
    Vec q12 = par->q12;
    Vec b = par->b;
    Vec K;
    
    Vec qv1;
    Vec qv2;
    inthelper_diagint momspacehelper;
    momspacehelper.integrator=par->integrator;
    
    momspacehelper.diag = par->diag;
    
    
    double momspace=0;
    if (dim == 8) // LO or type a
    {
        K = Vec (vec[6]*std::cos(vec[7]),vec[6]*std::sin(vec[7]));
        qv1= q12*0.5 - K*0.5;
        qv2= q12*(-0.5) - K*0.5;
        momspacehelper.q1 = qv1;
        momspacehelper.q2 = qv2;
        
        double loparvec[6]={vec[0], vec[1], vec[2], vec[3], vec[4], vec[5]};
        
        if (par->diag == DIAG_LO)
            momspace = inthelperf_mc_lo(loparvec, 6, &momspacehelper);
        else
            momspace = inthelperf_mc_diag2a(loparvec, 6, &momspacehelper);
        
    }
    else
    {
        K =  Vec(vec[9]*std::cos(vec[10]),vec[9]*std::sin(vec[10]));
        qv1= q12*0.5 - K*0.5;
        qv2= q12*(-0.5) - K*0.5;
        momspacehelper.q1 = qv1;
        momspacehelper.q2 = qv2;
        
        double parvec[9] = {vec[0],vec[1],vec[2],vec[3],vec[4],vec[5],vec[6],vec[7],vec[8]};
        momspace = inthelperf_mc_diag2b(parvec, 9, &momspacehelper);
        
    }

    
    
    // Ward
    if (qv1.LenSqr() < 1e-15 or qv2.LenSqr() < 1e-15)
        return 0;
    
 
    // Ward limit
    //if ((q-K*0.5).LenSqr() < 1e-5 or (q+K*0.5).LenSqr() < 1e-5)
    //    return 0;
    
    
    double res = momspace / std::pow(2.0*M_PI,2.);
    
    
    res *= std::cos(b*K);
    
    
    // Jacobian
    res *= K.Len();
    
    if (isnan(res))
    {
        //return 0;
        cerr << "NaN with K " << K << " q12 " << q12 << endl;
        //cerr << "Diag is " << diag_momentumspace << endl;
        cerr << "Argumets" << endl;
        cerr << qv1 << endl;
        cerr << qv2 << endl;
        cerr << "This probably means that you need more MC integration points" << endl;
        
        
        cerr << endl;
    }
    
    return res;
}






double DiagramIntegrator::MixedSpaceBruteForce(Diagram diag, Vec q12, Vec b)
{
    /*if (diag != DIAG_LO)
    {
        cerr << "DipoleAmplitudeBruteForce only supports LO at the moment" << endl;
        exit(1);
    }*/
    
    // Integrate over k, ktheta
    /*
    double lower[4] = {0.01, 0, 0.01,0};
    double upper[4] = {10,2.0*M_PI,10,2.0*M_PI};
     */
    // Integrata over the same variables as in the LO diagram + k,ktheta
    double KLIM = 12;
    double xlow=x;
    double xup = 0.999;
    
    mixed_space_helper helper;
    helper.q12=q12; helper.b=b; helper.integrator=this;
    helper.diag = diag;
    gsl_monte_function F;
       
    double *lower;
    double *upper;
    
    switch (diag) {
        case DIAG_2A:
        case DIAG_3A:
        case DIAG_3A_2:
        case DIAG_3B:
        case DIAG_3B_2:
        case DIAG_5A:
        case DIAG_5C:
        case DIAG_5C_1:
        case DIAG_LO:
            F.dim=8;
            lower = new double[F.dim];
            upper = new double [F.dim];
            lower[0]=lower[1]=lower[2]=lower[3]=0;
            lower[4]=lower[5]=xlow;
            lower[6]=0.0; lower[7]=0; // minK  theta_k
            
            upper[0]=upper[2]=KLIM;
            upper[1]=upper[3]=2.0*M_PI;
            upper[4]=upper[5]=xup;
            upper[6]=KLIM; upper[7]=2.0*M_PI; // maxK theta_k
            break;
        default:
            F.dim=11;
            lower = new double[F.dim];
            upper = new double [F.dim];
            lower[0]=lower[1]=lower[2]=lower[3]=lower[7]=lower[8]=0;
            lower[4]=lower[5]=xlow; lower[6]=x;
            lower[9]=0.0; lower[10]=0;
            
            upper[0]=upper[2]=upper[7]=KLIM;
            upper[1]=upper[3]=upper[8]=2.0*M_PI;
            upper[4]=upper[5]=upper[6]=xup;
            upper[9]=KLIM; upper[10]=2.0*M_PI;
            
    };
    
    
    
    F.f = inthelperf_mc_mixedspace;
    F.params = &helper;
    

    
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
        //cout << "# vegas warmup " << result << " +/- " << error << endl;
        int iter=0;
        do
        {
            gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS, rng, s, &result, &error);
            //cout << "# Vegas interation " << result << " +/- " << error << " chisqr " << gsl_monte_vegas_chisq(s) << endl;
            iter++;
        } while ((fabs( gsl_monte_vegas_chisq(s) - 1.0) > 0.4 or iter < 2) and iter < 6);
        gsl_monte_vegas_free(s);
    }
    else
        return 0;
    
    delete[] upper;
    delete[] lower;
    
    return result;
}







////
/////// Dipole ampiltude

struct dipole_helper
{
    DiagramIntegrator* integrator;
    Vec r;
    Vec b;
    Diagram diag;
};

double inthelperf_mc_dipole(double *vec, size_t dim, void* p)
{
    /*if (dim != 4) exit(1);*/
    dipole_helper *par = (dipole_helper*)p;
    Vec r = par->r;
    Vec b = par->b;
    
    Vec K;
    Vec q;
    Vec qv1;
    Vec qv2;
    inthelper_diagint momspacehelper;
    momspacehelper.integrator=par->integrator;
    
    momspacehelper.diag = par->diag;
    
    
    double momspace=0;
    if (dim == 10) // LO or type a
    {
        K = Vec(vec[6]*std::cos(vec[7]),vec[6]*std::sin(vec[7]));
        q = Vec(vec[8]*std::cos(vec[9]),vec[8]*std::sin(vec[9]));
        
        qv1= q - K*0.5;
        qv2= q*(-1) - K*0.5;
        
        if (qv1.LenSqr() < 1e-15 or qv2.LenSqr() < 1e-15)
            return 0;
        
        momspacehelper.q1 = qv1;
        momspacehelper.q2 = qv2;
        
        double loparvec[6]={vec[0], vec[1], vec[2], vec[3], vec[4], vec[5]};
        
        if (par->diag == DIAG_LO)
            momspace = inthelperf_mc_lo(loparvec, 6, &momspacehelper);
        else
            momspace = inthelperf_mc_diag2a(loparvec, 6, &momspacehelper);
        
    }
    else
    {
        K =  Vec(vec[9]*std::cos(vec[10]),vec[9]*std::sin(vec[10]));
        q = Vec(vec[11]*std::cos(vec[12]),vec[11]*std::sin(vec[12]));
        qv1= q - K*0.5;
        qv2= q*(-1) - K*0.5;
        if (qv1.LenSqr() < 1e-15 or qv2.LenSqr() < 1e-15)
            return 0;
        momspacehelper.q1 = qv1;
        momspacehelper.q2 = qv2;
        
        double parvec[9] = {vec[0],vec[1],vec[2],vec[3],vec[4],vec[5],vec[6],vec[7],vec[8]};
        momspace = inthelperf_mc_diag2b(parvec, 9, &momspacehelper);
        
    }
    
    
    double Klen=K.Len();
    double qlen=q.Len();

    
    // Ward limit
    //if ((q-K*0.5).LenSqr() < 1e-5 or (q+K*0.5).LenSqr() < 1e-5)
    //    return 0;
    
    
    double res = momspace / std::pow(2.0*M_PI,4.);
    
    res /= (qv1.LenSqr() * qv2.LenSqr());
    
    res *= std::cos(b*K) * (std::cos(r*q) - std::cos((r*K)/2.));
    
    /*double diag_momentumspace = par->integrator->IntegrateDiagram(par->diag, q - K*0.5, q*(-1) - K*(0.5));
    res *= diag_momentumspace;*/
    
    // Jacobian
    res *= Klen*qlen;
    
    if (isnan(res))
    {
        //return 0;
        cerr << "NaN with K " << K << " q " << q << endl;
        //cerr << "Diag is " << diag_momentumspace << endl;
        cerr << "Argumets" << endl;
        cerr << qv1 << endl;
        cerr << qv2 << endl;
        cerr << "This probably means that you need more MC integration points" << endl;
        
        
        cerr << endl;
    }
    
    return res;
}

// Color factor -g^2/2 Cf not included
double DiagramIntegrator::DipoleAmplitudeBruteForce(Diagram diag, Vec r, Vec b)
{

    // Integrate over k, ktheta, q, qhteta
    /*
    double lower[4] = {0.01, 0, 0.01,0};
    double upper[4] = {10,2.0*M_PI,10,2.0*M_PI};
     */
    // Integrata over the same variables as in the lO diagram + k,ktheta,q,qtheta
    double KLIM = 15;
    double xlow=x;
    double xup = 0.999;
    
    double *lower;
    double *upper;
    
    gsl_monte_function Ff;
    
    switch (diag) {
        case DIAG_2A:
        case DIAG_3A:
        case DIAG_3A_2:
        case DIAG_3B:
        case DIAG_3B_2:
        case DIAG_5A:
        case DIAG_5C:
        case DIAG_5C_1:
        case DIAG_UV_SUM:
        case DIAG_LO:
            Ff.dim=10;
            lower = new double[Ff.dim];
            upper = new double [Ff.dim];
            lower[0]=lower[1]=lower[2]=lower[3]=0;
            lower[4]=lower[5]=xlow;
            lower[6]=0; lower[7]=0; // minK  theta_k
            lower[8]=0; lower[9]=0;
            upper[0]=upper[2]=KLIM;
            upper[1]=upper[3]=2.0*M_PI;
            upper[4]=upper[5]=xup;
            upper[6]=KLIM; upper[7]=2.0*M_PI; // maxK theta_k
            upper[8]=KLIM; upper[9]=2.0*M_PI;
            break;
        default:
            Ff.dim=13;
            lower = new double[Ff.dim];
            upper = new double [Ff.dim];
            lower[0]=lower[1]=lower[2]=lower[3]=lower[7]=lower[8]=0;
            lower[4]=lower[5]=xlow; lower[6]=x;
            lower[9]=0.; lower[10]=0;
            lower[11]=0; lower[12]=0;
            upper[0]=upper[2]=upper[7]=KLIM;
            upper[1]=upper[3]=upper[8]=2.0*M_PI;
            
            upper[4]=upper[5]=upper[6]=xup;
            upper[9]=KLIM; upper[10]=2.0*M_PI;
            upper[11]=KLIM; upper[12]=2.0*M_PI;
            
    };
    
   
    
    
    dipole_helper helper;
    helper.r=r; helper.b=b; helper.integrator=this;
    helper.diag = diag;
    
       
    Ff.params = &helper;
    Ff.f = inthelperf_mc_dipole;

    
    double result,error;
    if (intmethod == MISER)
    {
        gsl_monte_miser_state *s = gsl_monte_miser_alloc(Ff.dim);
        gsl_monte_miser_integrate(&Ff, lower, upper, Ff.dim, MCINTPOINTS, rng, s, &result, &error);
        cout << "# Miser result " << result << " err " << error << " relerr " << std::abs(error/result) << endl;
        gsl_monte_miser_free(s);
    }
    else if (intmethod == VEGAS)
    {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(Ff.dim);
        gsl_monte_vegas_integrate(&Ff, lower, upper, Ff.dim, MCINTPOINTS/2, rng, s, &result, &error);
        //cout << "# vegas warmup " << result << " +/- " << error << endl;
        int iter=0;
        do
        {
            gsl_monte_vegas_integrate(&Ff, lower, upper, Ff.dim, MCINTPOINTS, rng, s, &result, &error);
            //cout << "# Vegas interation " << result << " +/- " << error << " chisqr " << gsl_monte_vegas_chisq(s) << endl;
            iter++;
        } while ((fabs( gsl_monte_vegas_chisq(s) - 1.0) > 0.4 or iter < 2) and iter < 6);
        gsl_monte_vegas_free(s);
    }
    else
        return 0;
    
    delete[] upper;
    delete[] lower;
    
    return result;
}


Interpolator* DiagramIntegrator::InitializeInterpolator()
{
    
    std::cerr << "Interpolator not implemented! " << std::endl;
    
    /*
     
    std::vector<double> qvals; std::vector<double> F;
    for (double q=minq; q<=maxq; q+=(maxq-minq)/npoints)
    {
        qvals.push_back(q);
        Vec l(0,q);
        F.push_back(F_int_B0(l, Vec(0,0), alpha, mf*mf));
    }
    Interpolator *interp = new Interpolator(qvals,F);
    interp->SetFreeze(true); interp->SetUnderflow(0); interp->SetOverflow(0);
     */
    return 0;
    
    
}


DiagramIntegrator::DiagramIntegrator()
{
    mf=0.2;
    intmethod = VEGAS;
    
    proton.SetBeta(0.55);
    proton.SetM(0.26);
   
    gsl_rng_env_setup ();
    
    F = new F_worker(20, 0.0001); // divisions accuracy

    const gsl_rng_type *T = gsl_rng_default;
    rng = gsl_rng_alloc (T);
    x=0.01;
    
    small_x = false;
    collinear_cutoff_uv_finite = false;
}


DiagramIntegrator::~DiagramIntegrator()
{
    delete F;
}


Diagram DiagramIntegrator::DiagramType(std::string str)
{
    for (int i=0; i < NUM_OF_DIAGRAMS; i++)
    {
        if (DIAGRAM_STRINGS[i] == str)
        {
            return DIAGRAMS[i];
        }
    }
    std::cerr << "Unknown diagram " << str << std::endl;
    exit(1);
}

bool DiagramIntegrator::Add_Q1Q2_exchange(Diagram diag)
{
    return false;
}

std::string DiagramIntegrator::InfoStr()
{
    std::stringstream ss;
    ss << "# x=" << x << endl <<
    "# Perturbative m_f=" << mf << endl
    << "# Proton wave function: " << WaveFunctionString(proton.GetWaveFunction()) << endl
    << "# Proton wave function params: mq=" << proton.GetM() << "GeV, beta=" << proton.GetBeta() <<" GeV, p=" << proton.GetP() << endl;
    ss << "# small-x limit: "; if (small_x) ss << "true"; else ss << "false"; ss << endl;
    if (collinear_cutoff_uv_finite)
        ss << "# IR regulator included in UV finite diagrams" << endl;
    
    return ss.str();
}
