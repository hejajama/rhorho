
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
    Diagram diag;
};

/*
 * LO diagram
 * [k1x,k1y,k2x,k2y,x1,x2,xg,kgx,kgy]
 */
double inthelperf_mc_lo(double *vec, size_t dim, void* p)
{
    if (dim != 6) exit(1);
    inthelper_diagint *par = (inthelper_diagint*)p;
    Vec k1(vec[0],vec[1]);
    Vec k2(vec[2], vec[3]);
    Vec q1 = par->q1;
    Vec q2 = par->q2;
    
    // Assume P=0
    Vec p1 = q1;
    Vec p2 = q2;
    
    double x1=vec[4];
    double x2=vec[5];
    if (x1+x2 >=1) return 0;
    
    //double x3 = 1.-x1-x2;
    
    //Vec k3 = (k1+k2)*(-1);
    
    double wf1 = par->integrator->GetProton().WaveFunction( k1, k2, x1,  x2);
    
    Vec K = (q1+q2)*(-1.);
    double wf2 = par->integrator->GetProton().WaveFunction(p1-q1-q2-K*x1, p2-K*x2, x1, x2)
        - par->integrator->GetProton().WaveFunction(p1-q1-K*x1, p2-q2-K*x2, x1, x2);
    
    double res = wf1*wf2;
    
    // 16pi^3 because NLO diagrams do not include 1/(16pi^3) prefactor
    // Python analysis notebook divides by 1/16pi^3
    return 16.*std::pow(M_PI,3.) * 1./2. * 1./4. * res / (8.0*x1*x2*(1.-x1-x2)*std::pow(2.0*M_PI,6.0));
    

}

/*
 * Diagram 2b
 * Vector components are
 * [k1x,k1y,k2x,k2y,x1,x2,xg,kgx,kgy]
 */
double inthelperf_mc_diag2b(double *vec, size_t dim, void* p)
{
    if (dim != 9) exit(1);
    inthelper_diagint *par = (inthelper_diagint*)p;
    
    Vec k1(vec[0],vec[1]);
    Vec k2(vec[2], vec[3]);
    Vec kg(vec[7],vec[8]);
    Vec q1 = par->q1;
    Vec q2 = par->q2;
    
    double x1=vec[4];
    double x2=vec[5];
    if (x1+x2 >= 1) return 0;
    
    double xg = vec[6];
    if (xg > std::min(x1,1.-x2)) return 0;
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
            ktilde_1 = k1 + (q1+q2)*x1-kg;
            ktilde_2 = k2 - (q1+q2)*(1.-x2) + kg;
            A = p1*z1 - kg;
            B = p2*z2 - (kg-q1-q2)*(1.-z2);
            norm=-1./6. * 6;
            
            break;
        case DIAG_3C:
            ktilde_1 =k1 + (q1+q2)*x1 - q1 - kg;
            ktilde_2 = k2 + (q1+q2)*x2 - q2 + kg;
            A = p1*z1 - kg;
            B = p2*z2 - (kg-q2)*(1.-z2);
            norm=1./12.*6;
            break;
        case DIAG_3D:
            ktilde_1 = k1 + (q1+q2)*x1 - kg;
            ktilde_2 = k2 - (q1+q2)*(1.-x2) + kg;
            A = p1*z1 - kg;
            B = (p2-q1)*z2 - (kg-q2)*(1.-z2);
            norm=1./12.*6.;
            break;
        case DIAG_6E_1:
            ktilde_1 = k1 - (q1+q2)*(1.-x1) - kg;
            ktilde_2 = k2 + (q1+q2)*x2 + kg;
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = 1./3.*6;
            break;
        case DIAG_6E_2:
            ktilde_1 = k1 - (q1+q2)*(1.-x1);
            ktilde_2 = k2 + (q1+q2)*x2 + kg;
            A = (p1-q1-q2)*z1 -kg;
            B = p2*z2 - kg*(1.-z2);
            norm = 1./3.6;
            break;
        case DIAG_6F_1:
            ktilde_1 = k1 + (q1+q2)*x1 - kg;
            ktilde_2 = k2 - (q1+q2)*(1.-x2) + kg;
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = 1./3. * 6;
            break;
        case DIAG_6F_2:
            ktilde_1 = k1 + (q1+q2)*x1 - kg;
            ktilde_2 = k2 - (q1+q2)*(1.-x2)+kg;
            A = p1*z1-kg;
            B = (p2-q1-q2)*z2 - kg*(1.-z2);
            norm = 1./3. * 6;
            break;
        case DIAG_6G_1:
            ktilde_1 = k1 + (q1+q2)*x1 - kg;
            ktilde_2 = k2 + (q1+q2)*x2 + kg;
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = CF*1./3. * 6;
            break;
        case DIAG_6G_2:
            ktilde_1 = k1 + (q1+q2)*x1 - kg;
            ktilde_2 = k2 + (q1+q2)*x2 + kg;
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = 1./3. * CF * 6;
            break;
        case DIAG_7H:
            ktilde_1 = k1 + (q1+q2)*x1 - q2 - kg;
            ktilde_2 = k2 + (q1+q2)*x2 - q1 + kg;
            A = p1*z1 - kg;
            B = (p2-q1)*z2 - kg*(1.-z2);
            norm = 1./3.*(0.5-CF) * 6;
            break;
        case DIAG_7J:
            ktilde_1 = k1 + (q1+q2)*x1 - q2 - kg;
            ktilde_2 = k2 + (q1+q2)*x2 + kg;
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = 1./3. * (CF-1./2. - 1./6.);
            break;
        case DIAG_7K:
            ktilde_1 = k1 + (q1+q2)*x1 - q1 - kg;
            ktilde_2 = k2 + (q1+q2)*x2 + kg;
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = 1./3. * (CF-1./2. - 1./6.);
            break;
        case DIAG_7L:
            ktilde_1 = k1 + (q1+q2)*x1 - kg;
            ktilde_2 = k2 + (q1+q2)*x2 - q2 + kg;
            A = p1*z1 - kg;
            B = (p2-q2)*z2 - kg*(1.-z2);
            norm = 1./3. * (CF-1./2.-1./6.);
            break;
        case DIAG_8H_1:
            ktilde_1 = k1 + (q1+q2)*x1 - kg - q2;
            ktilde_2 = k2 + (q1+q2)*x2 + kg - q1;
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = -2./9.*6;
            break;
        case DIAG_8H_2:
            ktilde_1 = k1 + (q1+q2)*x1 - kg - q2;
            ktilde_2 = k2 + (q1+q2)*x2 + kg - q1;
            A = (p1-q2)*z1-kg;
            B = (p2-q1)*z2 - kg*(1.-z2);
            norm = -2./9. * 6;
            break;
        case DIAG_8J_1:
            ktilde_1 = k1 + (q1+q2)*x1 - kg - q2;
            ktilde_2 = k2 + (q1+q2)*x2 + kg;
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = -1./3. * (CF-2./3.)*6;
            break;
        case DIAG_8J_2:
            ktilde_1 = k1 + (q1+q2)*x1 - kg - q2;
            ktilde_2 = k2 + (q1+q2)*x2 + kg;
            A = (p1-q2)*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            norm = -1./3. * (CF-2./3.);
            break;
        case DIAG_8L_1:
            ktilde_1 = k1 + (q1+q2)*x1 - kg;
            ktilde_2 = k2 + (q1+q2)*x2 + kg - q2;
            A = p1*z1 - kg;
            B = p2*z1 - kg*(1.-z2);
            norm = -1./3. * (CF-2./3);
            break;
        case DIAG_8L_2:
            ktilde_1 = k1 + (q1+q2)*x1 - kg;
            ktilde_2 = k2 + (q1+q2)*x2 + kg - q2;
            A = p1*z1 - kg;
            B = (p2-q2)*z2 - kg*(1.-z2);
            norm = -1./3. * (CF-2./3);
            break;
        default:
            cerr << "Unknown diagram " << par->diag << endl;
            exit(1);
            break;
    }
    
   
    double wf1 =par->integrator->GetProton().WaveFunction(k1, k2, x1, x2);
    double wf2 = par->integrator->GetProton().WaveFunction(ktilde_1,ktilde_2,x1-xg, x2+xg);
    
    double res = norm*wf1*wf2*f_xg*(A*B)/(A.LenSqr()*B.LenSqr());
    
    res /= xg;
    
    // Jacobian
    res /= 8.0*x1*x2*(1.-x1-x2)*std::pow(2.0*M_PI,6.0);
    
    return res;
    
}



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
    if (x1+x2 >=1) return 0;
    
    //double x3 = 1.-x1-x2;
    
    //Vec k3 = (k1+k2)*(-1);
    
    double wf1 = par->integrator->GetProton().WaveFunction( k1, k2, x1,  x2);
    
    Vec k12; Vec k22;
    Vec l; Vec l1;
    double norm=1; // Normalization factof * symmetry factor,not including g^4 / 16pi^3
    
    switch(par->diag)
    {
        case DIAG_2A:
            l=q1+q2;
            l1=Vec(0,0);
            k12 =k1 - (q1+q2)*(1.-x1);
            k22=k2 + (q1+q2)*x2;
            norm = -1./2. * 3.;
            break;
        case DIAG_3A:
            l=q1+q2;
            l1=q2;
            k12=k1 - (q1+q2)*(1.-x1);
            k22=k2 + (q1+q2)*x2;
            norm = -1./3.*3;
            break;
        case DIAG_3B:
            l=q2;
            l1=Vec(0,0);
            k12 = k1 + (q1+q2)*x1;
            k22 = k2 + (q1+q2)*x2;
            norm = -1./6. * 6;
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
            k12 = k1+(q1+q2)*x1;
            k22 = k2+(q1+q2)*x2;
            norm = 2./(3.*6.)*6;
            break;
        default:
            std::cerr << "Unknown diagram " << par->diag << " encountered!" << std::endl;
            exit(1);
    }
     

    double wf2 = par->integrator->GetProton().WaveFunction( k12, k22, x1, x2);
    
    double alpha =  par->integrator->GetX() / x1;
    double mf =  par->integrator->GetMf();
    
    if (alpha < 1e-8 or mf < 1e-8)
    {
        cerr << "Invalid alpha=" << alpha << ", mf=" << mf << " GeV!" << endl;
        exit(1);
    }
    
    
    double fintb = 0;
    if (par->integrator->UseInterpolator() == true)
        fintb = par->F_B_interpolator->Evaluate(l.Len());
    else
        fintb = par->integrator->GetF_worker()->F_int_B0(l, l1, alpha, mf*mf);
    double result = norm*wf1*wf2*fintb;
    
    if (isinf(result) or isnan(result))
    {
        cerr << "Result "<< result << " k1=" << k1 <<", k2=" << k2 << " wf1 " << wf1 << " wf2 " << wf2 << endl;
    }
    
    // Jacobian
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
    const double KLIM=9;
    double xlim=0.001;
    
    //double lower[6]={-KLIM,-KLIM,-KLIM,-KLIM,xlim,xlim};
    //double upper[6]={KLIM,KLIM,KLIM,KLIM,1.-xlim,1.-xlim};
    
    double *lower;
    double* upper;
    
    switch (diag) {
        case DIAG_2A:
        case DIAG_3A:
        case DIAG_3B:
        case DIAG_5A:
        case DIAG_5C:
        case DIAG_LO:
            F.dim=6;
            lower = new double[F.dim];
            upper = new double [F.dim];
            lower[0]=lower[1]=lower[2]=lower[3]=-KLIM;
            lower[4]=lower[5]=xlim;
            upper[0]=upper[1]=upper[2]=upper[3]=KLIM;
            upper[4]=upper[5]=1.-xlim;
            F.f=inthelperf_mc_diag2a;
            break;
        default: // Other diagrams are of this type
            F.dim=9;
            lower = new double[F.dim];
            upper = new double [F.dim];
            lower[0]=lower[1]=lower[2]=lower[3]=lower[7]=lower[8]=-KLIM;
            lower[4]=lower[5]=xlim; lower[6]=x;
            upper[0]=upper[1]=upper[2]=upper[3]=upper[7]=upper[8]=KLIM;
            upper[4]=upper[5]=upper[6]=1.-xlim;
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
    
    delete lower;
    delete upper;
    
    if (use_interpolator)
        delete F_b_interp;
    
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
    mf=0.1;
    intmethod = VEGAS;
    
    proton.SetBeta(0.55);
    proton.SetM(0.26);
    cout <<"# Initializing proton..." << endl;
    double n = proton.ComputeWFNormalizationCoefficient();
    cout << "#... done, normalization coef " << n  << endl;
    gsl_rng_env_setup ();
    
    F = new F_worker(20, 0.0001); // divisions accuracy

    const gsl_rng_type *T = gsl_rng_default;
    rng = gsl_rng_alloc (T);
    x=0.01;
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
    if (diag == DIAG_5C or diag == DIAG_7H or diag == DIAG_7J
        or diag == DIAG_7L or diag == DIAG_8H_1 or diag == DIAG_8H_2
        or diag == DIAG_8L_1 or diag == DIAG_8L_2
        or diag == DIAG_3A or diag == DIAG_3B
        or diag == DIAG_3C or diag == DIAG_3D)
        return true;
    return false;
}
