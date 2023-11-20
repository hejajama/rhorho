
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

double inthelperf_mc_finitesum(double *vec, size_t dim, void* p)
{
    if (dim != 9) exit(1);
    inthelper_diagint *par = (inthelper_diagint*)p;
    
   
    
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
        std::cerr << "NOTE: Using small-x limit, probably not ok???" << std::endl;
        xg = 0;
    }
    
    double z1,z2;
    z1 = xg/x1; z2 = xg /( x2+xg );
    
    
    /// We work in the frame where P=0
    Vec p1 = k1; Vec p2 = k2;
    
    double sum = 0;
    for (unsigned int di=FIRST_UV_FINITE; di <= LAST_UV_FINITE; di++)
    {
        Diagram diag = DIAGRAMS[di];
        
        Vec ktilde_1; Vec ktilde_2;
        Vec A,B;
        double f_xg=std::sqrt(x1*x2/((x1-xg)*(x2+xg))) * (1. - (z1+z2)/2. + z1*z2/6.);
        double norm=1; // normalization * symmetry factor

        // Note: this is a bit ugly, but the following is basically a copypaste from diagram_integrator.cpp
        // with the difference that color factors that would be otherwise taken into account by my 
        // notebook are also included here.
        
        switch (diag) {
            case DIAG_2B:
                ktilde_1 = k1 + q*x1  - kg + K*xg;
                ktilde_2 = k2 - q*(1.-x2) + kg - K*xg;
                A = p1*z1 - kg;
                B = p2*z2 - (kg-q1-q2)*(1.-z2);
                norm=-1./6. * 6 * COLOR_ADJ;
                break;
            case DIAG_3C:
                ktilde_1 =k1 + q*x1 - q1 - kg + K*xg;
                ktilde_2 = k2 + q*x2 - q2 + kg - K*xg;
                A = p1*z1 - kg;
                B = p2*z2 - (kg-q2)*(1.-z2);
                norm=1./12.*6*COLOR_ADJ;
                break;
            case DIAG_3C_2: // q1 <->q2 swapped
                ktilde_1 =k1 + q*x1 - q2 - kg + K*xg;
                ktilde_2 = k2 + q*x2 - q1 + kg - K*xg;
                A = p1*z1 - kg;
                B = p2*z2 - (kg-q1)*(1.-z2);
                norm=1./12.*6*COLOR_ADJ;
                break;
            case DIAG_3D:
                ktilde_1 = k1 + q*x1 - kg + K*xg;
                ktilde_2 = k2 - q*(1.-x2)+kg-K*xg;
                A = p1*z1 - kg;
                B = (p2-q1)*z2 - (kg-q2)*(1.-z2);
                norm=1./12.*6.*COLOR_ADJ;
                break;
            case DIAG_3D_2:
                ktilde_1 = k1 + q*x1 - kg + K*xg;
                ktilde_2 = k2 - q*(1.-x2)+kg-K*xg;
                A = p1*z1 - kg;
                B = (p2-q2)*z2 - (kg-q1)*(1.-z2);
                norm=1./12.*6.*COLOR_ADJ;
                break;
            case DIAG_6E: // (82): now included when symmetry factors fixed
                ktilde_1 = k1 - q*(1.-x1) - kg + K*xg;
                ktilde_2 = k2 + q*x2 + kg - K*xg;
                A = p1*z1 - kg;
                B = p2*z2 - kg*(1.-z2);
                norm = -CF/3.*6*COLOR_FUND;
                break;
            case DIAG_6E_1: // (85), previously  6E + 6E_1 = 0, but now not true with fixed symmetry factors
                ktilde_1 = k1 - q*(1.-x1) - kg + K*xg;
                ktilde_2 = k2 + q*x2 + kg - K*xg;
                A = p1*z1 - kg;
                B = p2*z2 - kg*(1.-z2);
                norm = CF/3. * 3*COLOR_FUND;
                break;
            case DIAG_6E_2: // Risto (86)
                ktilde_1 = k1 - q*(1.-x1) - kg + K*xg;
                ktilde_2 = k2 + q*x2 + kg - K*xg;
                A = (p1-q)*z1 -kg;
                B = p2*z2 - kg*(1.-z2);
                norm = CF/3. * 3*COLOR_FUND;   // Symmetry factor fixed
                break;

            case DIAG_6F: // Cancelled with with 6f'' originally, not anymore with fixed symmetry factors
                ktilde_1 = k1 + q*x1 - kg + K*xg;
                ktilde_2 = k2 - q*(1.-x2) + kg - K*xg;
                A = p1*z1 - kg;
                B = (p2-q)*z2 - kg*(1.-z2);
                norm = -CF/3.*6*COLOR_FUND;
                break;
            case DIAG_6F_1: // (87)
                ktilde_1 = k1 + q*x1 - kg + K*xg;
                ktilde_2 = k2 - q*(1.-x2) + kg - K*xg;
                A = p1*z1 - kg;
                B = p2*z2 - kg*(1.-z2);
                norm = CF/3. * 3*COLOR_FUND;
                break;
            case DIAG_6F_2: // Originally cancelled with 6f, not with fixed symmetry factors, (88)
                ktilde_1 = k1 + q*x1 - kg + K*xg;
                ktilde_2 = k2 - q*(1.-x2) + kg - K*xg;
                A = p1*z1-kg;
                B = (p2-q)*z2 - kg*(1.-z2);
                norm = CF/3. * 3*COLOR_FUND;
                break;
               // 6g+6g'+6g''=0 when symmetry factors are fixed now
            case DIAG_7H:
                ktilde_1 = k1 + q*x1 - q2 - kg + K*xg;
                ktilde_2 = k2 + q*x2 - q1 + kg - K*xg;
                A = p1*z1 - kg;
                B = (p2-q1)*z2 - kg*(1.-z2);
                norm = 1./3.*(0.5-CF) * 6*COLOR_FUND;
                break;
            case DIAG_7I: // 7H, q1 <-> q2 swap
                ktilde_1 = k1 + q*x1 - q1 - kg + K*xg;
                ktilde_2 = k2 + q*x2 - q2 + kg - K*xg;
                A = p1*z1 - kg;
                B = (p2-q2)*z2 - kg*(1.-z2);
                norm = 1./3.*(0.5-CF) * 6*COLOR_FUND;
                break;
            case DIAG_7J:
                ktilde_1 = k1 + q*x1 - q2 - kg + K*xg;
                ktilde_2 = k2 + q*x2 + kg - K*xg;
                A = p1*z1 - kg;
                B = p2*z2 - kg*(1.-z2);
                norm = 1./3. * (CF-1./2. - 1./6.)*6*COLOR_FUND;
                break;
            case DIAG_7K: // 7K, q1 <-> q2
                ktilde_1 = k1 + q*x1 - q1 - kg + K*xg;
                ktilde_2 = k2 + q*x2 + kg - K*xg;
                A = p1*z1 - kg;
                B = p2*z2 - kg*(1.-z2);
                norm = 1./3. * (CF-1./2. - 1./6.)*6*COLOR_FUND;
                break;
            case DIAG_7L:
                ktilde_1 = k1 + q*x1 - kg + K*xg;
                ktilde_2 = k2 + q*x2 - q2 + kg - K*xg;
                A = p1*z1 - kg;
                B = (p2-q2)*z2 - kg*(1.-z2);
                norm = 1./3. * (CF-1./2.-1./6.)*6*COLOR_FUND;
                break;
            case DIAG_7M: // 7L, q1 <-> q2
                ktilde_1 = k1 + q*x1 - kg + K*xg;
                ktilde_2 = k2 + q*x2 - q1 + kg - K*xg;
                A = p1*z1 - kg;
                B = (p2-q1)*z2 - kg*(1.-z2);
                norm = 1./3. * (CF-1./2.-1./6.)*6*COLOR_FUND;
                break;
            case DIAG_8H_1:
                ktilde_1 = k1 + q*x1 - kg - q2 + K*xg;
                ktilde_2 = k2 + q*x2 + kg - q1 - K*xg;
                A = p1*z1 - kg;
                B = p2*z2 - kg*(1.-z2);
                norm = -2./9.*3*COLOR_FUND;
                break;
            case DIAG_8H_2:
                ktilde_1 = k1 + q*x1 - kg - q2 + K*xg;
                ktilde_2 = k2 + q*x2 + kg - q1 - K*xg;
                A = (p1-q2)*z1-kg;
                B = (p2-q1)*z2 - kg*(1.-z2);
                norm = -2./9. * 3*COLOR_FUND;
                break;
            case DIAG_8I_1: // 8H_1, q1 <-> q2
                ktilde_1 = k1 + q*x1 - kg - q1 + K*xg;
                ktilde_2 = k2 + q*x2 + kg - q2 - K*xg;
                A = p1*z1 - kg;
                B = p2*z2 - kg*(1.-z2);
                norm = -2./9.*3*COLOR_FUND;
                break;
            case DIAG_8I_2: // 8H_2, q1 <-> q2
                ktilde_1 = k1 + q*x1 - kg - q1 + K*xg;
                ktilde_2 = k2 + q*x2 + kg - q2 - K*xg;
                A = (p1-q1)*z1-kg;
                B = (p2-q2)*z2 - kg*(1.-z2);
                norm = -2./9. * 3*COLOR_FUND;
                break;
            case DIAG_8J_1:
                ktilde_1 = k1 + q*x1 - kg - q2 + K*xg;
                ktilde_2 = k2 + q*x2 + kg - K*xg;
                A = p1*z1 - kg;
                B = p2*z2 - kg*(1.-z2);
                norm = -1./3. * (CF-2./3.)*3*COLOR_FUND;
                break;
            case DIAG_8J_2:
                ktilde_1 = k1 + q*x1 - kg - q2 + K*xg;
                ktilde_2 = k2 + q*x2 + kg - K*xg;
                A = (p1-q2)*z1 - kg;
                B = p2*z2 - kg*(1.-z2);
                norm = -1./3. * (CF-2./3.)*3*COLOR_FUND;
                break;
            case DIAG_8K_1: // 8J_1, q1 <-> q2
                ktilde_1 = k1 + q*x1 - kg - q1 + K*xg;
                ktilde_2 = k2 + q*x2 + kg - K*xg;
                A = p1*z1 - kg;
                B = p2*z2 - kg*(1.-z2);
                norm = -1./3. * (CF-2./3.)*3*COLOR_FUND;
                break;
            case DIAG_8K_2: // 8J_2, q1 <-> q2
                ktilde_1 = k1 + q*x1 - kg - q1 + K*xg;
                ktilde_2 = k2 + q*x2 + kg - K*xg;
                A = (p1-q1)*z1 - kg;
                B = p2*z2 - kg*(1.-z2);
                norm = -1./3. * (CF-2./3.)*3*COLOR_FUND;
                break;
            case DIAG_8L_1:
                ktilde_1 = k1 + q*x1 - kg + K*xg;
                ktilde_2 = k2 + q*x2 + kg - q2 - K*xg;
                A = p1*z1 - kg;
                B = p2*z2 - kg*(1.-z2);
                norm = -1./3. * (CF-2./3)*3*COLOR_FUND;
                break;
            case DIAG_8L_2:
                ktilde_1 = k1 + q*x1 - kg + K*xg;
                ktilde_2 = k2 + q*x2 + kg - q2 - K*xg;
                A = p1*z1 - kg;
                B = (p2-q2)*z2 - kg*(1.-z2);
                norm = -1./3. * (CF-2./3)*3*COLOR_FUND;
                break;
            case DIAG_8M_1: // 8L_1, q1<->q2
                ktilde_1 =  k1 + q*x1 - kg + K*xg;
                ktilde_2 = k2 + q*x2 + kg - q1 - K*xg;
                A = p1*z1 - kg;
                B = p2*z2 - kg*(1.-z2);
                norm = -1./3. * (CF-2./3)*3*COLOR_FUND;
                break;
            case DIAG_8M_2: // 8L_2, q1<->q2
                ktilde_1 = k1 + q*x1 - kg + K*xg;
                ktilde_2 = k2 + q*x2 + kg - q1 - K*xg;
                A = p1*z1 - kg;
                B = (p2-q1)*z2 - kg*(1.-z2);
                norm = -1./3. * (CF-2./3)*3*COLOR_FUND;
                break;
            default:
                cerr << "Unknown diagram in inthelperf_mc_diag2b: " << par->diag << endl;
                exit(1);
                break;
        }
        
        if (A.LenSqr() < 1e-15 or B.LenSqr() < 1e-15)
            return 0;
        
        
        
        double wf2 = par->integrator->GetProton().WaveFunction(ktilde_1,ktilde_2,x1-xg, x2+xg);
        
        double dotprod = 0;
        double mf = par->integrator->GetMf();
        if (par->integrator->CollinearCutoffUVFinite())
            dotprod = (A*B) / ( (A.LenSqr()+mf*mf)*(B.LenSqr()+mf*mf));
        else
            dotprod =(A*B)/(A.LenSqr()*B.LenSqr());
        
        //double res = norm*wf1*wf2*f_xg*dotprod;
        //res *= inv_xg; // same as res /= xg;
        sum +=norm*wf2*f_xg*dotprod;
    }
    
    double wf1 =par->integrator->GetProton().WaveFunction(k1, k2, x1, x2);
    
    double res = wf1*sum*inv_xg;
    // Jacobian
    res *= vec[0] * vec[2] * vec[7];
    res /= 8.0*x1*x2*(1.-x1-x2)*std::pow(2.0*M_PI,6.0);
    
    return res;
    
}

double inthelperf_mc_uvsum(double* vec, size_t dim, void* p)
{
    if (dim != 6) exit(1);
    inthelper_diagint *par = (inthelper_diagint*)p;
    
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
    
    double sum = 0;
    for (unsigned int di = FIRST_UV_DIV; di<= LAST_UV_DIV; di++)
    {
        Diagram diag = DIAGRAMS[di];
        
        Vec k12; Vec k22;
        Vec l; Vec l1;
        double norm=1; // Normalization factor * symmetry factor,not including g^4 / 16pi^3
        
        double alpha =  par->integrator->GetX() / x1;
        // Default, changed below if necessary
        switch(diag)
        {
            case DIAG_2A:
                l=q1+q2;
                l1=Vec(0,0);
                k12 =k1 - (q1+q2)*(1.-x1);
                k22=k2 + (q1+q2)*x2;
                norm = 2./3. * 3. * COLOR_ADJ;
                break;
            case DIAG_3A:
                l=q1+q2;
                l1=q1;
                k12=k1 - (q1+q2)*(1.-x1);
                k22=k2 + (q1+q2)*x2;
                norm = -1./3.*3 * COLOR_ADJ;
                break;
            case DIAG_3A_2:
                l=q1+q2;
                l1=q2;
                k12=k1 - (q1+q2)*(1.-x1);
                k22=k2 + (q1+q2)*x2;
                norm = -1./3.*3 * COLOR_ADJ;
                break;
            case DIAG_3B:
                l=q2;
                l1=Vec(0,0);
                k12 = k1 + (q1+q2)*x1-q2;
                k22 = k2 + (q1+q2)*x2-q1;
                norm = -1./6. * 6 * COLOR_ADJ;
                break;
            case DIAG_3B_2:
                l=q1;
                l1=Vec(0,0);
                k12 = k1 + (q1+q2)*x1-q1;
                k22 = k2 + (q1+q2)*x2-q2;
                norm = -1./6. * 6 * COLOR_ADJ;
                alpha = par->integrator->GetX() / x2;
                break;
            case DIAG_5A:
                l=q1+q2;
                l1=q1+q2;
                k12 = k1 - (q1+q2)*(1.-x1);
                k22 = k2 + (q1+q2)*x2;
                norm = 4.*CF/3.*3 * COLOR_FUND;
                break;
            case DIAG_5C:
                l=q2;
                l1=q2;
                k12 = k1+(q1+q2)*x1-q2;
                k22 = k2+(q1+q2)*x2-q1;
                norm = 2./(3.*6.)*6*COLOR_FUND;
                break;
            case DIAG_5C_1:
                l=q1;
                l1=q1;
                k12 = k1+(q1+q2)*x1-q1;
                k22 = k2+(q1+q2)*x2-q2;
                norm = 2./(3.*6.)*6*COLOR_FUND;
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
        sum += norm*wf2*fintb;
            
        /*if (isinf(sum) or isnan(sum))
        {
            cerr << "Result "<< result << " k1=" << k1 <<", k2=" << k2 << " wf1 " << wf1 << " wf2 " << wf2 << endl;
        }*/
        
    }
    
    sum *= wf1;
   
    // Jacobian
    sum *= vec[0]*vec[2];
    sum /= x1*x2*(1.-x1-x2)*8*std::pow(2.0*M_PI,6.);
   
    return 2.0*std::pow(M_PI,3.)*sum; // A21 gives 2pi^3
}
