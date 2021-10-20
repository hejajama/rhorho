//
//  odderon_ft.cpp
//  rholib
//
//  Created by Heikki Mäntysaari on 14.10.2021.
//

#include "diagram_integrator.hpp"
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




///////////////
// Mixed space brute force
////
/////// Dipole ampiltude

struct mixed_space_odderon_helper
{
    DiagramIntegrator* integrator;
    Vec q12;
    Vec q23;
    Vec b;
    Diagram diag;
};

double inthelperf_mc_odderon_mixedspace(double *vec, size_t dim, void* p)
{
    /*if (dim != 4) exit(1);*/
    mixed_space_odderon_helper *par = (mixed_space_odderon_helper*)p;
    Vec q12 = par->q12;
    Vec q23 = par->q23;
    Vec b = par->b;
    Vec K;
    
    Vec qv1;
    Vec qv2;
    Vec qv3;
    inthelper_diagint momspacehelper;
    momspacehelper.integrator=par->integrator;
    
    momspacehelper.diag = par->diag;
    
    
    double momspace=0;
    if (dim == 8) // LO or type a
    {
        K = Vec (vec[6]*std::cos(vec[7]),vec[6]*std::sin(vec[7]));
        
        qv1 = (q12*2. + q23 - K)*(1./3.);
        qv2 = (q12*(-1) + q23 - K)*(1./3.);
        qv3 = (q12 + q23*2. + K)*(-1./3.);
        
        momspacehelper.q1 = qv1;
        momspacehelper.q2 = qv2;
        momspacehelper.q3 = qv3;
        
        double loparvec[6]={vec[0], vec[1], vec[2], vec[3], vec[4], vec[5]};
        
        if (par->diag == ODDERON_LO)
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
    if (qv1.LenSqr() < 1e-5 or qv2.LenSqr() < 1e-5)
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






double DiagramIntegrator::OdderonG2b(Vec b, Vec q12, Vec q23, Diagram diag)
{

    
    // Integrate over K, ktheta


    double KLIM = 12;
    double MINK=0.001;
    double xlow=x;
    double xup = 0.999;
    
    mixed_space_odderon_helper helper;
    helper.q12=q12;
    helper.q23=q23;
    helper.b=b; helper.integrator=this;
    helper.diag = diag;
    gsl_monte_function F;
       
    double *lower;
    double *upper;
    
    switch (diag) {
        case ODDERON_LO:
        case ODDERON_DIAG_14:
        case ODDERON_DIAG_15:
        case ODDERON_DIAG_16:
        case ODDERON_DIAG_17:
        case ODDERON_DIAG_20:
        case ODDERON_DIAG_21:
        case ODDERON_DIAG_22:
        case ODDERON_DIAG_36:
        case ODDERON_DIAG_37:
        case ODDERON_DIAG_18:
        case ODDERON_DIAG_29:
        case ODDERON_DIAG_32:
        case ODDERON_DIAG_38:
        case ODDERON_DIAG_39:
        case ODDERON_DIAG_40:
        case ODDERON_DIAG_19:
        case ODDERON_DIAG_30:
        case ODDERON_DIAG_33:
        case ODDERON_DIAG_41:
        case ODDERON_DIAG_42:
        case ODDERON_DIAG_48:
        case ODDERON_DIAG_31:
        case ODDERON_DIAG_34:
        case ODDERON_DIAG_43:
        case ODDERON_DIAG_49:
            F.dim=8;
            lower = new double[F.dim];
            upper = new double [F.dim];
            lower[0]=lower[1]=lower[2]=lower[3]=-KLIM;
            lower[4]=lower[5]=xlow;
            lower[6]=MINK; lower[7]=0; // minK mink theta_k
            
            upper[0]=upper[1]=upper[2]=upper[3]=KLIM;
            upper[4]=upper[5]=xup;
            upper[6]=KLIM; upper[7]=2.0*M_PI; // minK mink theta_k
            break;
        default:
            F.dim=11;
            lower = new double[F.dim];
            upper = new double [F.dim];
            lower[0]=lower[1]=lower[2]=lower[3]=lower[7]=lower[8]=-KLIM;
            lower[4]=lower[5]=xlow; lower[6]=x;
            lower[9]=0.01; lower[10]=0;
            upper[0]=upper[1]=upper[2]=upper[3]=upper[7]=upper[8]=KLIM;
            upper[4]=upper[5]=upper[6]=xup;
            upper[9]=KLIM; upper[10]=2.0*M_PI;
            
    };
    
    
    
    F.f = inthelperf_mc_odderon_mixedspace;
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

