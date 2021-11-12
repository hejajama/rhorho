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



// Integrand where all uv finite odderons are summed

/*
 * UV finite diagrams
 * Vector components are
 * [k1x,k1y,k2x,k2y,x1,x2,xg,kgx,kgy]
 */
double inthelperf_mc_finitesum(double *vec, size_t dim, void* p)
{
    if (dim != 9) exit(1);
    inthelper_diagint *par = (inthelper_diagint*)p;
    
    Vec k1(vec[0]*std::cos(vec[1]),vec[0]*std::sin(vec[1]));
    Vec k2(vec[2]*std::cos(vec[3]), vec[2]*std::sin(vec[3]));
    Vec kg(vec[7]*std::cos(vec[8]),vec[7]*std::sin(vec[8]));
    Vec q1 = par->q1;
    Vec q2 = par->q2;
    Vec q3 = par->q3;
    Vec q = q1+q2+q3;
    Vec K = q1*(-1) - q2 - q3;
    
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
    
   

    
    
    double f_xg=std::sqrt(x1*x2/((x1-xg)*(x2+xg))) * (1. - (z1+z2)/2. + z1*z2/6.);
   
    
    double sum=0;
//#pragma omp parallel for reduction(+:sum)
    for (unsigned int di=FIRST_UV_FINITE_ODDERON; di <= LAST_UV_FINITE_ODDERON; di++)
    {
        Diagram diag = DIAGRAMS[di];
        bool include_A2_B2=false;
        Vec A2(0,0);
        Vec B2(0,0); // J operator for odderon needs these
        Vec ktilde_1; Vec ktilde_2;
        Vec A,B;
        double norm=1; // normalization * symmetry factor
        
    
        switch (diag) {
        
        //// Odderon
        case ODDERON_DIAG_69:
            A = p1*z1-kg;
            B = p2*z2 - (kg-(q1+q2))*(1.-z2);
            ktilde_1 = k1 + q*x1 - q3 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q1+q2) + kg - K*xg;
            norm = -1./4. * NC/2. * 6.;
            break;
            
        case ODDERON_DIAG_72:
            A = p1*z1-kg;
            B = (p2 - q3)*z2 - (kg-(q1+q2))*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q + kg - K*xg;
            norm = -1./4. * NC/2. * 6; // TODO TR(T^a T^b D^c - T^a T^b T^c)
            break;
        case ODDERON_DIAG_75:
            A =p1*z1-kg;
            B = p2*z2 - (kg-(q1+q2))*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q1+q2) + kg - K*xg;
            norm = 1./2. * NC/2. * 6.; // TODO TR T^a T^b D^c
            break;
            
        case ODDERON_DIAG_70:
            A =p1*z1-kg;
            B=p2*z2 - (kg-(q1+q3))*(1.-z2);
            ktilde_1 = k1 + q*x1 - q2 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q1+q3) + kg - K*xg;
            norm = -1./4.* NC/2. * 6;
            break;
            
        case ODDERON_DIAG_78:
            A = p1*z1-kg;
            B = p2*z2 - (kg-q1)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg - (q2+q3);
            ktilde_2 = k2 + q*x2 - q1 + kg - K*xg;
            norm = NC/2. * 1./4. * 6.;
            break;
            
        case ODDERON_DIAG_81:
            A =p1*z1-kg;
            B = (p2-q3)*z2 - (kg-q1)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg - q2;
            ktilde_2 = k2 + q*x2 - (q1+q3) + kg - K*xg;
            norm = NC/(2.*4.) * 6.;
            break;
            
        case ODDERON_DIAG_90:
            A =p1*z1-kg;
            B = p2*z2 - (kg-q1)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg - q2;
            ktilde_2 = k2 + q*x2 - q1 + kg - K*xg;
            norm= - 6. * NC/(2.*4.);
            break;
            
            
        case ODDERON_DIAG_71:
            A =p1*z1-kg;
            B = p2*z2 - (kg - (q2+q3))*(1.-z2);
            ktilde_1 = k1 + q*x1 - q1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q2+q3) + kg - K*xg;
            norm = -1./4. * NC/2. * 6;
            break;
        
        case ODDERON_DIAG_79:
            A =p1*z1-kg;
            B = p2*z2 - (kg - q2)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg - (q1+q3);
            ktilde_2 = k2 + q*x2 - q2 + kg - K*xg;
            norm = NC/2. * 1./4. * 6.;
            break;
        
        case ODDERON_DIAG_82:
            A = p1*z1-kg;
            B = (p2-q3)*z2 - (kg-q2)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg - q1;
            ktilde_2 = k2 + q*x2 - (q2+q3) + kg - K*xg;
            norm = NC/2. * 1./4. * 6.;
            break;
            
        case ODDERON_DIAG_91:
            A = p1*z1-kg;
            B = p2*z2 - (kg-q2)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg - q1;
            ktilde_2 = k2 + q*x2 - q2 + kg - K*xg;
            norm = -NC/(2.*4.) * 6.;
            break;
            
        case ODDERON_DIAG_73:
            A = p1*z1 - kg;
            B = (p2-q2)*z2 - (kg - (q1+q3))*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q + kg - K*xg;
            norm = -1./4. * NC/2. * 6.;
            break;
            
        case ODDERON_DIAG_84:
            A = p1*z1 - kg;
            B = (p2-q2)*z2 - (kg-q1)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg - q3;
            ktilde_2 = k2 + q*x2 - (q1+q2) + kg - K*xg;
            norm = NC/2. * 1./4. * 6.;
            break;
            
        case ODDERON_DIAG_87:
            A = p1*z1 - kg;
            B = (p2-(q2+q3))*z2 - (kg-q1)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q + kg - K*xg;
            norm = NC/2. * 1./4. * 6.;
            break;
            
        case ODDERON_DIAG_96:
            A = p1*z1 -kg;
            B = (p2-q2)*z2 - (kg-q1)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q1+q2) + kg - K*xg;
            norm = (-NC/4. + 1./4. * NC/2. )* 6.;
            break;
            
            
        case ODDERON_DIAG_74:
            A = p1*z1 - kg;
            B = (p2 - q1)*z2 - (kg - (q2+q3))*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q + kg - K*xg;
            norm = -1./4. * NC/2. * 6.;
            break;
            
        case ODDERON_DIAG_85:
            A =p1*z1 - kg;
            B = (p2-q1)*z2 - (kg-q2)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg - q3;
            ktilde_2 = k2 + q*x2 - (q1+q2) + kg - K*xg;
            norm = NC/2. * 1./4. * 6.;
            break;
            
        case ODDERON_DIAG_88:
            A = p1*z1 - kg;
            B = (p2-(q1+q3))*z2 - (kg-q2)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q + kg - K*xg;
            norm = NC/2. * 1./4. * 6.;
            break;
        
        case ODDERON_DIAG_97:
            A = p1*z1 - kg;
            B = (p2-q1)*z2 - (kg-q2)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q1+q2) + kg - K*xg;
            norm = -6. * NC/(4.*2.);
            break;
            
        case ODDERON_DIAG_89:
            A = p1*z1 - kg;
            B = (p2 - (q1+q2))*z2 - (kg-q3)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q + kg - K*xg;
            norm = 6. * NC/2. * 1./4.;
            break;
        
        case ODDERON_DIAG_117:
            A = p1*z1 - kg;
            B = (p2 - (q1+q2))*z2 - kg*(1.-z2);
            ktilde_1 = k1+q*x1 - q3 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q1+q2) + kg - K*xg;
            norm = -(CF - 1./2.) * 1./4. * 6;
            break;
        
        case ODDERON_DIAG_118:
            A = p1*z1 - kg;
            B = (p2-q)*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q + kg - K*xg;
            norm = -CF*1./4. * 6;
            break;
            
        case ODDERON_DIAG_119:
            A = p1*z1 - kg;
            B = (p2 - (q1+q2))*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q1+q2) + kg - K*xg;
            norm = (2.*CF - 1./2. - NC/2.) * 1./4. * 6.;
            break;
        
        case ODDERON_DIAG_76:
            A = p1*z1 - kg;
            B = p2*z2 - (kg-(q1+q3))*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q1+q3) + kg - K*xg;
            norm = 1./2. * NC/2. * 6.;
            break;
        
        case ODDERON_DIAG_93:
            A = p1*z1 - kg;
            B = p2*z2 - (kg-q1)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg - q3;
            ktilde_2 = k2 + q*x2 - q1 + kg - K*xg;
            norm = (-NC*1./4. + 1./4. * NC/2.)*6.;
            break;
        
        case ODDERON_DIAG_99:
            A = p1*z1 - kg;
            B = (p2-q3)*z2 - (kg-q1)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q1+q3) + kg - K*xg;
            norm = (-NC/4. + 1./4.*NC/2.)*6.;
            break;
            
        case ODDERON_DIAG_77:
            A = p1*z1 - kg;
            B = p2*z2 - (kg - (q2+q3))*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q2+q3) + kg - K*xg;
            norm = 1./2. * NC/2. * 6.;
            break;
        
        case ODDERON_DIAG_94:
            A = p1*z1 - kg;
            B = p2*z2 - (kg-q2)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg - q3;
            ktilde_2 = k2 + q*x2 - q2 + kg -K*xg;
            norm = (-NC*1./4. + 1./4.*NC/2.)*6.;
            break;
            
        case ODDERON_DIAG_100:
            A = p1*z1 - kg;
            B = (p2-q3)*z2 - (kg-q2)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q2+q3) + kg - K*xg;
            norm = (-NC/4. + 1./4.*NC/2.) * 6.;
            break;
            
        case ODDERON_DIAG_86:
            A = p1*z1 - kg;
            B = (p2-q1)*z2 - (kg-q3)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg - q2;
            ktilde_2 = k2 + q*x2 - (q1+q3) + kg - K*xg;
            norm = NC/2. * 1./4. * 6.;
            break;
        
        case ODDERON_DIAG_114:
            A = p1*z1 - kg;
            B = (p2-q1)*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - (q2+q3) - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q1 + kg - K*xg;
            norm = -(CF-1./2.)*1./4.*6.;
            break;
            
        case ODDERON_DIAG_115:
            A = p1*z1 - kg;
            B = (p2-(q1+q3))*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - q2 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q1+q3) + kg - K*xg;
            norm = -(CF-1./2.)*1./4.*6.;
            break;
            
        case ODDERON_DIAG_116:
            A = p1*z1 - kg;
            B = (p2-q1)*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - q2 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q1 + kg - K*xg;
            norm = ((CF-NC/2.)*1./4. + (CF-1.)*1./4.)*6.;
            break;
            
            
        case ODDERON_DIAG_83:
            A = p1*z1 - kg;
            B = (p2-q2)*z2 - (kg-q3)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg - q1;
            ktilde_2 = k2 + q*x2 - (q2+q3) + kg - K*xg;
            norm = NC/2. * 1./4. * 6.;
            break;
            
        case ODDERON_DIAG_108:
            A = p1*z1 - kg;
            B = (p2-q2)*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - (q1+q3) - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q2 + kg - K*xg;
            norm = -(CF-1./2.)*1./4.*6.;
            break;
            
        case ODDERON_DIAG_109:
            A = p1*z1 - kg;
            B = (p2 - (q2+q3))*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - q1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q2+q3) + kg - K*xg;
            norm = -(CF-1./2.)*1./4. * 6;
            break;
            
        case ODDERON_DIAG_110:
            A = p1*z1 - kg;
            B = (p2-q2)*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - q1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q2 + kg - K*xg;
            norm = ((CF-1.)*1./4. + (CF-NC/2.)*1./4.)*6.;
            break;
            
        case ODDERON_DIAG_80:
            A = p1*z1 - kg;
            B = p2*z2 - (kg-q3)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg - (q1+q2);
            ktilde_2 = k2 + q*x2 - q3 + kg - K*xg;
            norm = NC/2. * 1./4. * 6.;
            break;
            
        case ODDERON_DIAG_105:
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            ktilde_1 = k1 - q*(1.-x1) - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            norm = -CF*1./4. * 6.;
            break;
            
        case ODDERON_DIAG_106:
            A = p1*z1 - kg;
            B = (p2-q3)*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - (q1+q2) - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q3 + kg - K*xg;
            norm = -(CF-1./2.)*1./4.*6.;
            break;
            
        case ODDERON_DIAG_107:
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - (q1+q2) - kg + K*xg;
            ktilde_2 = k2 + q*x2  + kg - K*xg;
            norm = (CF - 1./2. - 1./(2.*NC))*1./4. * 6.;
            break;
            
        case ODDERON_DIAG_92:
            A = p1*z1 - kg;
            B = p2*z2 - (kg-q3)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg - q1;
            ktilde_2 = k2 + q*x2 - q3 + kg - K*xg;
            norm = (-NC*1./4. + 1./4.*NC/2.)*6.;
            break;
            
        case ODDERON_DIAG_111:
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - (q1+q3) - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            norm = (2.*CF - 1./2. - NC/2.)*1./4. * 6.;
            break;
        
        case ODDERON_DIAG_112:
            A = p1*z1 - kg;
            B = (p2-q3)*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - q1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q3 + kg - K*xg;
            norm = ((CF-1.)*1./4. + (CF-NC/2.)*1./4.)*6.;
            break;
            
        case ODDERON_DIAG_113:
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - q1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            norm = (2.*CF - 1./2. - NC/2.)*1./4. * 6.;
            break;
            
        case ODDERON_DIAG_95:
            A = p1*z1 - kg;
            B = p2*z2 - (kg-q3)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg - q2;
            ktilde_2 = k2 + q*x2 - q3 + kg - K*xg;
            norm = (-NC/4. + 1./4.*NC/2.)*6.;
            break;
            
        case ODDERON_DIAG_123:
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - (q2+q3) - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            norm = (2.*CF - 1./2. - NC/2.)*1./4. * 6.;
            break;
            
        case ODDERON_DIAG_124:
            A = p1*z1 - kg;
            B = (p2-q3)*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - q2 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q3 + kg - K*xg;
            norm = ((CF - 1.)*1./4. + (CF-NC/2.)*1./4.)*6.;
            break;
            
        case ODDERON_DIAG_125:
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - q2 - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            norm = (2.*CF - 1./2. - NC/2.)*1./4. * 6.;
            break;
        
        case ODDERON_DIAG_129:
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - q3 - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg -K*xg;
            norm = (2.0*CF - 1./2. - NC/2.)*1./4. * 6.;
            break;
            
        case ODDERON_DIAG_130:
            A = p1*z1 - kg;
            B = (p2-q3)*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q3 + kg - K*xg;
            norm = (2.*CF-1./2. - NC/2.)*1./4. * 6.;
            break;
            
        case ODDERON_DIAG_131:
            A = p1*z1 - kg;
            B = p2*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            norm = CF*(2.-NC)*1./4.*6.;
            break;
            
        case ODDERON_DIAG_98:
            A = p1*z1 - kg;
            B = (p2-q1)*z2 - (kg-q3)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q1+q3) + kg - K*xg;
            norm = (-NC/4. + 1./4. * NC/2.) * 6.;
            break;
            
        case ODDERON_DIAG_120:
            A = p1*z1 - kg;
            B = (p2-q1)*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - q3 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q1 + kg - K*xg;
            norm = ((CF-1.)*1./4. + (CF-NC/2.)*1./4.)*6.;
            break;
            
        case ODDERON_DIAG_121:
            A = p1*z1 - kg;
            B = (p2-(q1+q3))*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q1+q3) + kg - K*xg;
            norm = (2.0*CF - 1./2. - NC/2.)*1./4. * 6.;
            break;
            
        case ODDERON_DIAG_122:
            A = p1*z1 - kg;
            B = (p2-q1)*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q1 + kg - K*xg;
            norm = (2.0*CF - 1./2. - NC/2.)*1./4. * 6.;
            break;
            
        case ODDERON_DIAG_101:
            A =p1*z1 - kg;
            B = (p2-q2)*z2 - (kg-q3)*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q2+q3) + kg - K*xg;
            norm = (-NC/4. + 1./4. * NC/2.) * 6.;
            break;
            
        case ODDERON_DIAG_126:
            A = p1*z1 - kg;
            B = (p2-q2)*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - q3 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q2 + kg - K*xg;
            norm = ((CF-1.)*1./4. + (CF-NC/2.)*1./4.)*6.;
            break;
            
        case ODDERON_DIAG_127:
            A = p1*z1 - kg;
            B = (p2-(q2+q3))*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q2+q3) + kg - K*xg;
            norm = (2.0*CF - 1./2. - NC/2.) * 1./4. * 6.;
            break;
            
        case ODDERON_DIAG_128:
            A = p1*z1 - kg;
            B = (p2-q2)*z2 - kg*(1.-z2);
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q2 + kg - K*xg;
            norm = (2.0*CF - 1./2. - NC/2.)*1./4. * 6.;
            break;
            
            
        case ODDERON_DIAG_133:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = A;
            B2 = (p1-q)*z1 - kg;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - q - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            norm = CF*1./4. * 6.;
            break;
        
        case ODDERON_DIAG_134:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = (p2-q3)*z2 - kg*(1.-z2);
            B2 = (p1-(q1+q2))*z1 - kg;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - (q1+q2) - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q3 + kg - K*xg;
            norm = (CF - (NC+1.)/2.)*1./4. * 6.;
            break;
            
        case ODDERON_DIAG_135:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = A;
            B2 = (p1-(q1+q2))*z1 - kg;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - (q1+q2) - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            norm = -(2.0*CF-(NC+1.)/2.) * 1./4.*6.;
            break;
            
        case ODDERON_DIAG_136:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = (p2-q2)*z2 - kg*(1.-z2);
            B2 = (p1-(q1+q3))*z1 - kg;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - (q1+q3) - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q2 + kg - K*xg;
            norm = (CF - (NC+1.)/2.) * 1./4. * 6.;
            break;
            
        case ODDERON_DIAG_137:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = (p2-(q2+q3))*z2 - kg*(1.-z2);
            B2 = (p1-q1)*z1 - kg;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - q1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q2+q3) + kg - K*xg;
            norm = (CF - (NC+1.)/2.) * 1./4. * 6.;
            break;
            
        case ODDERON_DIAG_138:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = (p2-q2)*z2 - kg*(1.-z2);
            B2 = (p1-q1)*z1 - kg;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - q1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q2 + kg - K*xg;
            norm = -(CF - (NC+1.)/2.)*(1./4. + 1./4.) * 6.;
            break;
            
        case ODDERON_DIAG_139:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = A;
            B2 = (p1 - (q1+q3))*z1 - kg;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - (q1+q3) - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            norm = -(2.0*CF - (NC+1.)/2.) * 1./4. * 6.;
            break;
        
        case ODDERON_DIAG_140:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = (p2-q3)*z2 - kg*(1.-z2);
            B2 = (p1-q1)*z1 - kg;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - q1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q3 + kg - K*xg;
            norm = -(CF - (NC+1.)/2.) * (1./4. + 1./4.)*6.;
            break;
        
        case ODDERON_DIAG_141:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = A;
            B2 = (p1-q1)*z1 - kg;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - q1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            norm = -(2.0*CF - (NC + 1.)/2.) * 1./4. * 6.;
            break;
            
        case ODDERON_DIAG_142:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = (p2-q1)*z2 - kg*(1.-z2);
            B2 = (p1 - (q2+q3))*z1 - kg;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - (q2+q3) - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q1 + kg - K*xg;
            norm = (CF - (NC+1.)/2.)*1./4.* 6.;
            break;
            
        case ODDERON_DIAG_143:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = (p2 - (q1+q3))*z2 - kg*(1.-z2);
            B2 = (p1-q2)*z1 - kg;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - q2 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q1+q3) + kg - K*xg;
            norm = (CF - (NC+1.)/2.)*1./4. * 6.;
            break;
            
        case ODDERON_DIAG_144:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = (p2-q1)*z2 - kg*(1.-z2);
            B2 = (p1-q2)*z1 - kg;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - q2 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q1 + kg - K*xg;
            norm = -(CF- (NC+1.)/2.) * (1./4. + 1./4.) * 6.;
            break;
            
        case ODDERON_DIAG_145:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = (p2 - (q1+q2))*z2 - kg*(1.-z2);
            B2 = (p1-q3)*z1 - kg;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - q3 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q1+q2) + kg - K*xg;
            norm = (CF - (NC+1.)/2.) * 1./4. * 6.;
            break;
            
        case ODDERON_DIAG_146:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = (p2-q)*z2 - kg*(1.-z2);
            B2 = B;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q + kg - K*xg;
            norm = CF*1./4. * 6.;
            break;
            
        case ODDERON_DIAG_147:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = (p2-(q1+q2))*z2 - kg*(1.-z2);
            B2=B;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - kg +K*xg;
            ktilde_2 = k2 + q*x2 - (q1+q2) + kg - K*xg;
            norm = -(2.*CF - (NC+1.)/2.)*1./4. * 6.;
            break;
            
        case ODDERON_DIAG_148:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = (p2-q1)*z2 - kg*(1.-z2);
            B2 = (p1-q3)*z1 - kg;
            include_A2_B2=true;
            ktilde_1 = k1+q*x1 - q3 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q1 + kg -K*xg;
            norm = -(CF - (NC+1.)/2.) * (1./4. + 1./4.)*6.;
            break;
            
            
        case ODDERON_DIAG_149:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = (p2-(q1+q3))*z2 - kg*(1.-z2);
            B2=B;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q1+q3) + kg - K*xg;
            norm = -(2.0*CF - (NC+1.)/2.) * 1./4. * 6.;
            break;
            
        
        case ODDERON_DIAG_150:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = (p2-q1)*z2 - kg*(1.-z2);
            B2 = B;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2  =k2 + q*x2 - q1 + kg - K*xg;
            norm = -(2.0*CF - (NC+1.)/2.) * 1./4. * 6.;
            break;
            
        case ODDERON_DIAG_151:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = A;
            B2 = (p1-(q2+q3))*z1 - kg;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - (q2+q3) - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            norm = -(2.0*CF - (NC+1.)/2.)*1./4. * 6.;
            break;
            
        case ODDERON_DIAG_152:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = (p2-q3)*z2 - kg*(1.-z2);
            B2 = (p1-q2)*z1 - kg;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - q2 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q3 + kg - K*xg;
            norm = -(CF - (NC+1.)/2.) * (1./4. + 1./4.)*6.;
            break;
            
        case ODDERON_DIAG_153:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = A;
            B2 = (p1-q2)*z1 - kg;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - q2 - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            norm = -(2.0*CF - (NC+1.)/2.) * 1./4. * 6.;
            break;
            
        case ODDERON_DIAG_154:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = (p2-q2)*z2 - kg*(1.-z2);
            B2 = (p1-q3)*z1 - kg;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - q3 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q2 + kg - K*xg;
            norm = -(CF - (NC+1.)/2.) * (1./4. + 1./4.)*6.;
            break;
            
        case ODDERON_DIAG_155:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = (p2 - (q2+q3))*z2 - kg*(1.-z2);
            B2 = B;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - (q2+q3) + kg - K*xg;
            norm = -(2.0 * CF - (NC+1.)/2.) * 1./4. * 6.;
            break;
            
        case ODDERON_DIAG_156:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = (p2-q2)*z2 - kg*(1.-z2);
            B2 = B;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q2 + kg - K*xg;
            norm = -(2.0*CF - (NC+1.)/2.) * 1./4. * 6.;
            break;
            
        case ODDERON_DIAG_157:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = A;
            B2 = (p1-q3)*z1 - kg;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - q3 - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            norm = -(2.0*CF - (NC+1.)/2.)*1./4. * 6.;
            break;
            
        case ODDERON_DIAG_158:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            A2 = (p2-q3)*z2 - kg*(1.-z2);
            B2 = B;
            include_A2_B2=true;
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 - q3 + kg - K*xg;
            norm = -(2.0*CF - (NC+1.)/2.) * 1./4. * 6.;
            break;
            
        case ODDERON_DIAG_159:
            A = p2*z2 - kg*(1.-z2);
            B = p1*z1 - kg;
            ktilde_1 = k1 + q*x1 - kg + K*xg;
            ktilde_2 = k2 + q*x2 + kg - K*xg;
            norm = -2.0*CF*(2.-NC)*1./4. * 6.;
            break;
            
            
        default:
            cerr << "Unknown diagram in odderon finite sum: " << par->diag << endl;
            exit(1);
            break;
        }
        double dotprod=0;
        if (A.LenSqr() < 1e-6 or B.LenSqr() < 1e-6)
            dotprod=0;
        else
        {
            if (include_A2_B2 == true)
                dotprod =(A*B)/(A.LenSqr()*B.LenSqr()) + (A2*B2)/(A2.LenSqr()*B2.LenSqr());
            else
                dotprod =(A*B)/(A.LenSqr()*B.LenSqr());
        }
        double wf2 = wf2 = par->integrator->GetProton().WaveFunction(ktilde_1,ktilde_2,x1-xg, x2+xg);
        
        sum += norm*dotprod*wf2;
        
    } // end loop over diags
    
    
    
   
    double wf1 =par->integrator->GetProton().WaveFunction(k1, k2, x1, x2);
    
    

    
   
    if (par->integrator->CollinearCutoffUVFinite())
    {
        cerr << "CollinearCutoff for UVFinite odderons not supported! " << endl;
        exit(1);
    }
       
    
    
    double res = wf1*f_xg*sum;
    
    res *= inv_xg; // same as res /= xg;
    
    // Jacobian
    res *= vec[0]*vec[2]*vec[7];
    
    res /= 8.0*x1*x2*(1.-x1-x2)*std::pow(2.0*M_PI,6.0);
    
    return res;
    
}


double intehelperf_mixed_uvsum(double* vec, size_t dim, void* p)
{
    
    
    if (dim != 6) exit(1);
    inthelper_diagint *par = (inthelper_diagint*)p;
    Vec k1(vec[0]*std::cos(vec[1]),vec[0]*std::sin(vec[1]));
    Vec k2(vec[2]*std::cos(vec[3]),vec[2]*std::sin(vec[3]));
    Vec q1 = par->q1;
    Vec q2 = par->q2;
    Vec q3 = par->q3;
    
    Vec q = q1+q2+q3;
    
    double x1=vec[4];
    double x2=vec[5];
    
    double x3 = 1.-x1-x2;
    double x = par->integrator->GetX();
    if (x3 >= 1 or x3 < x) return 0;
    if (x1+x2 >=1) return 0;
    
    double wf1 = par->integrator->GetProton().WaveFunction(k1, k2, x1,  x2);
    
    Vec k12; Vec k22;
    Vec l; Vec l1;
    double norm=1; // Normalization factor * symmetry factor,not including g^4 / 16pi^3
    
    double alpha =  par->integrator->GetX() / x1;
    // Default, changed below if necessary
    
    const double Nc=3;
    
    
    double sum=0;
//#pragma omp parallel for reduction(+:sum)
    for (unsigned int di=FIRST_UV_DIV_ODDERON; di <= LAST_UV_DIV_ODDERON; di++)
    {
        Diagram diag = DIAGRAMS[di];
        
        switch(par->diag)
        {
               
            // **** ODD DIAGRAMS, ODDERON ****
            // NOTE: In all odderon diags I factor out 2g^5 / (3 * 16pi^3) * d^abc
            // And f^abc parts are dropped
                // TODO* 1/4 in (6)
                // Todo (2pi)^3? see (160)
            case ODDERON_DIAG_14:
                l = q1 + q2 + q3;
                l1 = q3;
                k12 = k1 - (q1 + q2 + q3)*(1.-x1);
                k22 = k2 + (q1+q2+q3)*x2;
                norm = Nc/4. * 3; // 3 is symmetry factor

                // todo: alpha???
                break;
            case ODDERON_DIAG_15:
                l = q1+q2+q3;
                l1 = q2;
                k12 = k1 - (q1+q2+q3)*(1.-x1);
                k22 = k2 + (q1+q2+q3)*x2;
                norm = Nc/4. * 3;
                break;
            
            case ODDERON_DIAG_16:
                l = q1+q2+q3;
                l1 = q1;
                k12 = k1 - (q1+q2+q3)*(1.-x1);
                k22 = k2 + (q1+q2+q3)*x2;
                norm = Nc/4. * 3;
                break;
                
            case ODDERON_DIAG_17:
                l = q1+q2;
                l1=Vec(0,0);
                k12 = k1 + (q1+q2+q3)*x1 - (q1 + q2);
                k22 = k2 + (q1+q2+q3)*x2 - q3;
                norm = -Nc/8.0 * 6;
                break;
            
            case ODDERON_DIAG_20:
                l = q1+q2+q3;
                l1 = q2 + q3;
                k12 = k1 + (q1+q2+q3)*x1 - (q1+q2+q3);
                k22 = k2 + (q1+q2+q3)*x2;
                norm = -Nc/4. * 3;
                break;
                
            case ODDERON_DIAG_21:
                l = q1+q2+q3;
                l1 = q1 + q3;
                k12 = k1 + (q1+q2+q3)*x1 - (q1+q2+q3);
                k22 = k2 + (q1+q2+q3)*x2;
                norm = -Nc/4. * 3;
                break;
            
            case ODDERON_DIAG_22:
                l = q1+q2+q3;
                l1 = q1 + q2;
                k12 = k1 + (q1+q2+q3)*x1 - (q1+q2+q3);
                k22 = k2 + (q1+q2+q3)*x2;
                norm = -Nc/4. * 3;
                break;
            
            case ODDERON_DIAG_36:
                l = q1+q2+q3;
                l1 = q1+q2+q3;
                k12 = k1 - (q1+q2+q3)*(1.-x1);
                k22 = k2 + (q1+q2+q3)*x2;
                norm =  3. * 2.*CF/4.;
                break;
                
            case ODDERON_DIAG_37:
                l = q1+q2;
                l1 = q1+q2;
                k12 = k1 + (q1+q2+q3)*x1 - (q1+q2);
                k22 = k2 + (q1+q2+q3)*x2 - q3;
                norm = 6. / (2.*Nc*4.);
                break;
                
                
               ///
            case ODDERON_DIAG_18:
                l = q1+q3;
                l1 = Vec(0,0);
                k12 = k1 + q*x1-(q1+q3);
                k22 = k2 + q*x2 - q2;
                norm = -6.*Nc/(2.*4.);
                break;
            
            case ODDERON_DIAG_29:
                l = q1;
                l1 = Vec(0,0);
                k12 = k1 + q*x1 - q1;
                k22 = k2 + q*x2 - (q2+q3);
                norm = -6.*Nc/(2.*4.);
                break;
                
            case ODDERON_DIAG_32:
                l = q1;
                l1 = Vec(0,0);
                k12 = k1 + (q1+q2+q3)*x1 - q1;
                k22 = k2 + (q1+q2+q3)*x2-q2;
                norm = 6.*Nc/(2.*4.) * 2.;
                break;
            
            case ODDERON_DIAG_38:
                l=q1+q3;
                l1 = q1+q3;
                k12 = k1 + (q1+q2+q3)*x1-(q1+q3);
                k22 = k2 + (q1+q2+q3)*x2 - q2;
                norm = 6./(2.*Nc*4.);
                break;
                
            case ODDERON_DIAG_39:
                l=q1;
                l1=q1;
                k12 = k1 + (q1+q2+q3)*x1 - q1;
                k22 = k2 + (q1+q2+q3)*x2 - (q2+q3);
                norm= 6./(2.*Nc*4.);
                break;
                
            case ODDERON_DIAG_40:
                l=q1;
                l1=q1;
                k12 = k1 + (q1+q2+q3)*x1 - q1;
                k22 = k2 + (q1+q2+q3)*x2 - q2;
                norm = -6./(2.*Nc*4.) * 2.;
                break;
                
            //
                
            case ODDERON_DIAG_19:
                l=q2+q3;
                l1=Vec(0,0);
                k12 = k1 + q*x1 - (q2+q3);
                k22 = k2 + q*x2 - q1;
                norm=-6*Nc/(2*4);
                break;
            
            case ODDERON_DIAG_30:
                l=q2;
                l1=Vec(0,0);
                k12 = k1 + q*x1 - q2;
                k22 = k2 + q*x2 - (q1+q3);
                norm=-6*Nc/(2.*4.);
                break;
                
            case ODDERON_DIAG_33:
                l=q2;
                l1=Vec(0,0);
                k12 = k1 + q*x1 - q2;
                k22 = k2 + q*x2 - q1;
                norm=6.*Nc/(2.*4.)*2;
                break;
                
            case ODDERON_DIAG_41:
                l=q2+q3;
                l1=q2+q3;
                k12=k1+q*x1 - (q2+q3);
                k22 = k2 + q*x2 - q1;
                norm=6./(2.*Nc*4.);
                break;
            
            case ODDERON_DIAG_42:
                l=q2;
                l1=q2;
                k12=k1+q*x1-q2;
                k22 = k2 + q*x2 - (q1+q3);
                norm=6./(2.*Nc*4.);
                break;
            
            case ODDERON_DIAG_48:
                l=q2;
                l1=q2;
                k12=k1+q*x1-q2;
                k22=k2+q*x2-q1;
                norm=-6./(2.*Nc*4.)*2.;
                break;
                
            
            case ODDERON_DIAG_31:
                l=q3;
                l1=Vec(0,0);
                k12=k1+q*x1-q3;
                k22=k2 + q*x2-(q1+q2);
                norm=-6.*Nc/(2.*4.);
                break;
                
            case ODDERON_DIAG_34:
                l=q3;
                l1=Vec(0,0);
                k12=k1+q*x1-q3;
                k22=k2+q*x2-q1;
                norm=6.*Nc/(2.*4.)*2.;
                break;
                
            case ODDERON_DIAG_43:
                l=q3;
                l1=q3;
                k12=k1+q*x1-q3;
                k22=k2+q*x2-(q1+q2);
                norm=6./(4.*2.*Nc);
                break;
            
            case ODDERON_DIAG_49:
                l=q3;
                l1=q3;
                k12=k1+q*x1-q3;
                k22=k2+q*x2-q1;
                norm=-6./(4.*2.*Nc)*2.;
                break;
                
            
            default:
                std::cerr << "Unknown diagram " << par->diag << " encountered in inthelperf_mc_diag2a!" << std::endl;
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
            std::cerr << "Small-x limit not supported!" << endl;
            exit(1);
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
        
        //cout << l1 << " " << l << " " << wf1 << " " << wf2 << " " << fintb << endl;
       
        sum += result; // A21 gives 2pi^3
        
    }
    
    return sum* 2.0*std::pow(M_PI,3.) * vec[0]*vec[2]/( x1*x2*(1.-x1-x2)*8*std::pow(2.0*M_PI,6.) );
}


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
        qv1 = (q12*2. + q23 - K)*(1./3.);
        qv2 = (q12*(-1) + q23 - K)*(1./3.);
        qv3 = (q12 + q23*2. + K)*(-1./3.);
        
        momspacehelper.q1 = qv1;
        momspacehelper.q2 = qv2;
        momspacehelper.q3 = qv3;
        
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
            lower[0]=lower[1]=lower[2]=lower[3]=0;
            lower[4]=lower[5]=xlow;
            lower[6]=0; lower[7]=0; // minK mink theta_k
            
            upper[0]=upper[2]=KLIM;
            upper[1]=upper[3]=2.0*M_PI;
            upper[4]=upper[5]=xup;
            upper[6]=KLIM; upper[7]=2.0*M_PI; // minK mink theta_k
            break;
        default:
            F.dim=11;
            lower = new double[F.dim];
            upper = new double [F.dim];
            lower[0]=lower[1]=lower[2]=lower[3]=lower[7]=lower[8]=0;
            
            lower[4]=lower[5]=xlow; lower[6]=x;
            lower[9]=0.01; lower[10]=0;
            
            upper[0]=upper[2]=upper[7]=KLIM;
            upper[1]=upper[3]=upper[8]=2.0*M_PI;
            
            
            upper[4]=upper[5]=upper[6]=xup;
            upper[9]=KLIM; upper[10]=2.0*M_PI;
            
    };
    
    
    
    F.f = inthelperf_mc_odderon_mixedspace;
    F.params = &helper;
    

    
    double result,error;
    if (intmethod == MISER)
    {
        cerr << "You should not use MISER!" << endl;
        exit(1);
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
        } while ( (fabs( gsl_monte_vegas_chisq(s) - 1.0) > 0.2 or iter < 3) and iter < 7);
        gsl_monte_vegas_free(s);
    }
    else
        return 0;
    
    delete[] upper;
    delete[] lower;
    
    return result;
}




//// Odderon in coordinate space
///

////
/////// Dipole ampiltude


// https://arxiv.org/pdf/2001.04516.pdf (13)
double inthelperf_mc_odderon(double *vec, size_t dim, void* p)
{
    /*if (dim != 4) exit(1);*/
    dipole_helper *par = (dipole_helper*)p;
    Vec r = par->r;
    Vec b = par->b;
    
    Vec q1;
    Vec q2;
    Vec q3;
    inthelper_diagint momspacehelper;
    momspacehelper.integrator=par->integrator;
    
    momspacehelper.diag = par->diag;
    
    
    double momspace=0;
    if (dim == 12) // LO or type a
    {
        q1 = Vec(vec[6]*std::cos(vec[7]),vec[6]*std::sin(vec[7]));
        q2 = Vec(vec[8]*std::cos(vec[9]),vec[8]*std::sin(vec[9]));
        q3 = Vec(vec[10]*std::cos(vec[11]),vec[10]*std::sin(vec[11]));
        

        
        if (q1.LenSqr() < 1e-6 or q3.LenSqr() < 1e-6)
        {
            return 0; // Ward
        }
            
        
        momspacehelper.q1 = q1;
        momspacehelper.q2 = q2;
        momspacehelper.q3 = q3;
        
        double loparvec[6]={vec[0], vec[1], vec[2], vec[3], vec[4], vec[5]};
        
        if (par->diag == ODDERON_LO)
            momspace = inthelperf_mc_lo(loparvec, 6, &momspacehelper);
        else
            momspace = inthelperf_mc_diag2a(loparvec, 6, &momspacehelper);
        
    }
    else
    {
        q1 = Vec(vec[9]*std::cos(vec[10]),vec[9]*std::sin(vec[10]));
        q2 = Vec(vec[11]*std::cos(vec[12]),vec[11]*std::sin(vec[12]));
        q3 = Vec(vec[13]*std::cos(vec[14]),vec[13]*std::sin(vec[14]));
        
        if (q1.LenSqr() < 1e-6 or q3.LenSqr() < 1e-6)
            return 0; // Ward
        
        momspacehelper.q1 = q1;
        momspacehelper.q2 = q2;
        momspacehelper.q3 = q3;
        
        double parvec[9] = {vec[0],vec[1],vec[2],vec[3],vec[4],vec[5],vec[6],vec[7],vec[8]};
        momspace = inthelperf_mc_diag2b(parvec, 9, &momspacehelper);
        
    }
    
    
    Vec K = (q1+q2+q3)*(-1);
    
    
    double res = momspace / std::pow(2.0*M_PI,6.);
    
    res /= (q1.LenSqr() * q2.LenSqr() * q3.LenSqr());
    
    res *= (std::sin(r*q1 + (r*K)*0.5) - 1./3.*std::sin((r*K)*0.5));
    
    res *= -std::sin(b*K); // Imaginary part
    // So this actually computes i*T_ggg
    
    
    // Jacobian
    res *= q1.Len()*q2.Len()*q3.Len();
    
    res *= 5./18.;
    
    // Note: factor 1/4 which is the difference between <rho rho rho> and G3 is not included here
    
    if (isnan(res))
    {
        //return 0;
        //err << "NaN with K " << K << " q " << q << endl;
        //cerr << "Diag is " << diag_momentumspace << endl;
        //cerr << "Argumets" << endl;
        //cerr << qv1 << endl;
        //cerr << qv2 << endl;
        cerr << "NaN, this probably means that you need more MC integration points" << endl;
        
        
        cerr << endl;
    }
    
    return res;
}

// Color factor -g^2/2 Cf not included
double DiagramIntegrator::OdderonAmplitude(Diagram diag, Vec r, Vec b)
{

    // Integrate over q1, qtheta, q2, q2theta, q3, q3theta

    // Integrata over the same variables as in the lO diagram + q1, qtheta, q2, q2theta, q3, q3theta
    double KLIM = 12;
    double KMIN=0.001;
    double xlow=x;
    double xup = 0.999;
    
    double *lower;
    double *upper;
    
    cerr << "Päivitä napakoordinaatteihin\n";
    exit(1);
    
    gsl_monte_function Ff;
    
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
            Ff.dim=12;
            ///TODO PÄIVITÄ NAPAKOORDINAATTEIHIN
            lower = new double[Ff.dim];
            upper = new double [Ff.dim];
            lower[0]=lower[1]=lower[2]=lower[3]=-KLIM;
            lower[4]=lower[5]=xlow;
            lower[6]=KMIN; lower[7]=0; // minq min_theta_q
            lower[8]=KMIN; lower[9]=0; // minq min_theta_q
            lower[10]=0.01; lower[11]=0;
            upper[0]=upper[1]=upper[2]=upper[3]=KLIM;
            upper[4]=upper[5]=xup;
            upper[6]=KLIM; upper[7]=2.0*M_PI; // minK mink theta_k
            upper[8]=KLIM; upper[9]=2.0*M_PI;
            upper[10]=KLIM; upper[11]=2.0*M_PI;
            break;
        default:
            Ff.dim=15;
            lower = new double[Ff.dim];
            upper = new double [Ff.dim];
            lower[0]=lower[1]=lower[2]=lower[3]=lower[7]=lower[8]=-KLIM;
            lower[4]=lower[5]=xlow; lower[6]=x;
            lower[9]=KMIN; lower[10]=0;
            lower[11]=KMIN; lower[12]=0;
            lower[13]=KMIN; lower[14]=0;
            upper[0]=upper[1]=upper[2]=upper[3]=upper[7]=upper[8]=KLIM;
            upper[4]=upper[5]=upper[6]=xup;
            upper[9]=KLIM; upper[10]=2.0*M_PI;
            upper[11]=KLIM; upper[12]=2.0*M_PI;
            upper[13]=KLIM; upper[14]=2.0*M_PI;
            
    };
    
   
    
    
    dipole_helper helper;
    helper.r=r; helper.b=b; helper.integrator=this;
    helper.diag = diag;
    
       
    Ff.params = &helper;
    Ff.f = inthelperf_mc_odderon;

    
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
        } while ((fabs( gsl_monte_vegas_chisq(s) - 1.0) > 0.3 or iter < 3) and iter < 7);
        gsl_monte_vegas_free(s);
    }
    else
        return 0;
    
    delete[] upper;
    delete[] lower;
    
    return result;
}
