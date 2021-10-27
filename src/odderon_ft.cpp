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
    
    Vec k1(vec[0],vec[1]);
    Vec k2(vec[2], vec[3]);
    Vec kg(vec[7],vec[8]);
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
    res /= 8.0*x1*x2*(1.-x1-x2)*std::pow(2.0*M_PI,6.0);
    
    return res;
    
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
        } while ((fabs( gsl_monte_vegas_chisq(s) - 1.0) > 0.3 or iter < 2) and iter < 7);
        gsl_monte_vegas_free(s);
    }
    else
        return 0;
    
    delete[] upper;
    delete[] lower;
    
    return result;
}

