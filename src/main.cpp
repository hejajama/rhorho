#include <iostream>
#include <string>
#include <sstream>
#include <gsl/gsl_errno.h>
#include <cmath>
#include "vector.hpp"
#include "gitsha1.h"

#include "functions.hpp"
#include "diagram_integrator.hpp"
using namespace std;

double StrToReal(std::string str)
{
    std::stringstream buff(str);
    double tmp;
    buff >> tmp;
    return tmp;
}
long int StrToInt(std::string str)
{
    std::stringstream buff(str);
    double  tmp;
    buff >> tmp;
    return static_cast<long int>(tmp);
}

enum MODE
{
    ONEDIM,
    TWODIM, 
    WARD,
    FOURDIM,
    DIPOLE_BRUTEFORCE,
    DIPOLE_BRUTEFORCE_ANGLEDEP,
    MIXED_SPACE_BRUTEFORCE,
    MIXED_SPACE_BRUTEFORCE_ANGLEDEP
};

void handler (const char * reason,
const char * file,
int line,
int gsl_errno)
{
    cerr << "# Error " << gsl_errno << " line " << line << " file " << file << " reason " << reason << endl;
    
}

int main(int argc, char* argv[])
{
    gsl_set_error_handler(handler);
    
    cout << "# Command: ";
    for (int i=1; i<argc; i++)
        cout << argv[i] << " ";
    cout << endl;
    cout << "# Git version " << g_GIT_SHA1 << " local repo " << g_GIT_LOCAL_CHANGES << " main build " << __DATE__  << " " << __TIME__ << endl;
    
    
    DiagramIntegrator *integrator = new DiagramIntegrator;
    long int mcintpoints = 1e6;
    string diagram = "LO";
    MODE mode = ONEDIM;
    
    double q12=0.5;
    double theta_b_q = 0;
    double b=0;
    double r=1;
    double alpha_s=0;
    double x=0;

    bool scale_wf_params=false;
    
    for (int i=1; i< argc; i++)
    {
        if (string(argv[i])=="-x")
        {
            integrator->SetX(StrToReal(argv[i+1]));
            x=StrToReal(argv[i+1]);
        }
        else if (string(argv[i])=="-beta")
            integrator->GetProton().SetBeta(StrToReal(argv[i+1]));
        else if (string(argv[i])=="-wavef_mass")
            integrator->GetProton().SetM(StrToReal(argv[i+1]));
        else if (string(argv[i])=="-perturbative_mass")
            integrator->SetPerturbativeMass(StrToReal(argv[i+1]));
        else if (string(argv[i])=="-mcintpoints")
            mcintpoints = StrToInt(argv[i+1]);
        else if (string(argv[i])=="-diag")
            diagram = argv[i+1];
        else if (string(argv[i])=="-powerlaw")
            integrator->GetProton().SetWaveFunction(Power);
        else if (string(argv[i])=="-harmonic_oscillator")
            integrator->GetProton().SetWaveFunction(HarmoinicOscillator);
        else if (string(argv[i])=="-1d")
            mode = ONEDIM;
        else if (string(argv[i])=="-2d")
            mode = TWODIM;
        else if (string(argv[i])=="-4d")
            mode = FOURDIM;
        else if (string(argv[i])=="-ward")
            mode = WARD;
        else if (string(argv[i])=="-q12")
            q12 = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-theta_b_q")
            theta_b_q = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-dipole_bruteforce")
            mode = DIPOLE_BRUTEFORCE;
        else if (string(argv[i])=="-dipole_bruteforce_angledep")
            mode = DIPOLE_BRUTEFORCE_ANGLEDEP;
        else if (string(argv[i])=="-mixed_space_bruteforce")
            mode = MIXED_SPACE_BRUTEFORCE;
        else if (string(argv[i])=="-mixed_space_bruteforce_angledep")
        mode = MIXED_SPACE_BRUTEFORCE_ANGLEDEP;
        else if (string(argv[i])=="-b")
            b = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-r")
            r = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-smallx")
            integrator->SetSmallX(true);
        else if (string(argv[i])=="-regulate_uv_finite")
            integrator->SetCollinearCutoffUVFinite(true);
        else if (string(argv[i])=="-scale_beta")
        {
            scale_wf_params=true;
            alpha_s = StrToReal(argv[i+1]);
        }
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unknown parameter " << argv[i] << endl;
            exit(1);
        }
    
    }
    
    if (scale_wf_params)
    {
        cout <<"# Scaling wave function parameters: x and alpha_s dependent beta" << endl;
        cout <<"# Note that this data now only supports alphas=" << alpha_s << endl;
        integrator->GetProton().ScaleWaveFunctionParameters(alpha_s, x);
    }
    integrator->GetProton().ComputeWFNormalizationCoefficient();
    
    
    
   
    Diagram diag = integrator->DiagramType(diagram);
    
   
    integrator->UseInterpolator(false);
    
    cout << integrator->InfoStr();
    
    if ( diag==DIAG_2A or diag == DIAG_3A or diag == DIAG_3A_2 or diag == DIAG_3B or diag == DIAG_3B_2 or diag == DIAG_5A or diag == DIAG_5C or diag == DIAG_5C_1)
        mcintpoints /= 70; // there is one more intergal
                            // within the MC integral
    
    cout << "# Diagram " << diagram << " id " << diag << " mcintpoints " << mcintpoints << endl << "#" << endl;
    integrator->SetMCIntPoints(mcintpoints);
    
    
    
    // Test limit q1->0, q2->K, should vanish
    if (mode == WARD)
    {
        Vec K(1,0);
        for (double q=0.49; q<=0.5; q+=0.001)
        {
            Vec q1(q,0);
            double d = integrator->IntegrateDiagram(diag, q1-K*0.5, q1*(-1)-K*0.5);
            cout << (q1-K*0.5).Len() << " " << d << endl;
        }
    }
    
    if (mode == ONEDIM)
    {
    // q1=q2

	const double MINQ=0.05; //0.02*2;
	const double MAXQ = 1.06; //3;
	const double QSTEP = 0.10;
	const int QPOINTS = static_cast<int>((MAXQ-MINQ)/QSTEP);
	double *res = new double[QPOINTS];

	#pragma omp parallel for
        for (int qi=0; qi <= QPOINTS; qi++)
        {
	   double q = MINQ + qi*QSTEP;

	    // Parallel
        // Vec q1(q/2.,0);
       //  Vec q2(q/2.,0);
            
	    // Perpend
	    Vec q1(q/std::sqrt(2.),0);
	    Vec q2(0, q/std::sqrt(2.));
            double d = integrator->IntegrateDiagram(diag, q1, q2);
            if (integrator->Add_Q1Q2_exchange(diag))
            {
                cout << "#... adding cross graph q1<->q2" << endl;
                d +=integrator->IntegrateDiagram(diag, q1, q2);
            }
           
           res[qi]=d; 
        }

	for (int qi=0; qi <= QPOINTS; qi++)
        {
	   double q = MINQ + qi*QSTEP;
	   cout << q << " " << res[qi] << endl;
   	}
    }
    if (mode == TWODIM)
    {
    
    // LO paper reproduce at finite q12
    
        Vec q12v(q12*cos(theta_b_q), q12*sin(theta_b_q));
        
        cout << "# q12 = " << q12v << " orientation " << theta_b_q << endl;
        
        const double MAXK = 5;
        const int KPOINTS = 49;
        const double kstep =static_cast<double>(2*MAXK)/(KPOINTS-1);
        
        double *result = new double[KPOINTS*KPOINTS];
#pragma omp parallel for collapse(2)
        for (int kyi = 0; kyi < KPOINTS; kyi++)
        {
            for (int kxi = 0; kxi < KPOINTS; kxi ++ )
             {
                 double kx = -MAXK + kxi*kstep;
                 double ky = -MAXK + kyi*kstep;
                 Vec K(kx,ky);
                 Vec q1 = (q12v - K)*0.5;
                 Vec q2 = (q12v + K)*(-0.5);
                 
                 double d = integrator->IntegrateDiagram(diag, q1, q2);
                 if (integrator->Add_Q1Q2_exchange(diag))
                 {
                     cout << "#... adding cross graph q1<->q2" << endl;
                     d +=integrator->IntegrateDiagram(diag, q1, q2);
                 }
                 
                 result[kyi*KPOINTS+kxi] = d;
             }
        }
        
        cout <<"# kx   ky   result" << endl;
       for (int kxi = 0; kxi < KPOINTS; kxi++)
        {
            for (int kyi = 0; kyi < KPOINTS; kyi ++ )
             {
                 double kx = -MAXK + kxi*kstep;
                 double ky = -MAXK + kyi*kstep;
                 cout << kx << " " << ky << " " << result[kyi*KPOINTS+kxi] << endl;
             }
        }
        
        delete[] result;
    }
    
    if (mode == FOURDIM)
    {
        const double MAXK = 3.0;
        const int KPOINTS = 35;
        const double kstep =static_cast<double>(2*MAXK)/KPOINTS;
        cout << "# G(q - 1/2*k, -q - 1/2*k) /( (q-1/2k)^2(q+1/2k)^2" << endl;
        cout << "# kx  ky  qx  qy  diagram" << endl;
        
        double *result = new double[KPOINTS*KPOINTS*KPOINTS*KPOINTS];
        cerr << "NOTE: 4D GRID CALCULATION IS NOT THREAD SAFE! ONLY LO and NLO TYPE B WORK!" << endl;
#pragma omp parallel for collapse(3)
        for (int kyi=0; kyi < KPOINTS; kyi++)
        {
            for (int kxi=0; kxi < KPOINTS; kxi++)
            {
                for (int qyi=0; qyi < KPOINTS; qyi++)
                {
                    for (int qxi=0; qxi < KPOINTS; qxi++)
                    {
                        double ky = -MAXK + kyi*kstep;
                        double kx = -MAXK + kxi*kstep;
                        double qy = -MAXK + qyi*kstep;
                        double qx = -MAXK + qxi*kstep;
                
                        Vec q(qx,qy);
                        Vec K(kx,ky);
                        Vec q1 = q - K*0.5;
                        Vec q2 = q*(-0.5) + K*(-0.5);
                        
                        double res = 0;
                        
                        // Ward->0 if one of the momenta vanishes
                        if (q1.Len() > 1e-4 and q2.Len() > 1e-4)
                        {
                        
                            double d = integrator->IntegrateDiagram(diag, q1, q2);
                            if (integrator->Add_Q1Q2_exchange(diag))
                            {
                                cout << "#... adding cross graph q1<->q2" << endl;
                                d +=integrator->IntegrateDiagram(diag, q1, q2);
                            }
                             res =d/(q1.LenSqr()*q2.LenSqr());
                        }
                        
                        if (isnan(res))
                        {
                            cerr << "Note NaN at kx=" << kx << ", ky=" << ky << ", qx=" << qx <<", qy=" << qy << endl;
                            cerr << "q-K/2:" << endl << q1 << endl;
                            cerr << "q+K/2:" << endl << q2 << endl;
                            res=0;
                        }
                        
                        int i = kyi*pow(KPOINTS,3) + kxi*pow(KPOINTS,2) + qyi*KPOINTS + qxi;
                        result[i] = res;
                       
                        
                    }
                }
            }
        }
        
        for (int kyi=0; kyi < KPOINTS; kyi++)
        {
            for (int kxi=0; kxi < KPOINTS; kxi++)
            {
                for (int qyi=0; qyi < KPOINTS; qyi++)
                {
                    for (int qxi=0; qxi < KPOINTS; qxi++)
                    {
                        double ky = -MAXK + kyi*kstep;
                        double kx = -MAXK + kxi*kstep;
                        double qy = -MAXK + qyi*kstep;
                        double qx = -MAXK + qxi*kstep;
                        
                        int i = kyi*pow(KPOINTS,3) + kxi*pow(KPOINTS,2) + qyi*KPOINTS + qxi;
                        cout << ky << " " << kx << " " << qy << " " << qx << " " << result[i] << endl;
                    }
                }
            }
        }
        
        delete[] result;
    }
    
    else if (mode == DIPOLE_BRUTEFORCE)
    {
        Vec bv(b*std::cos(theta_b_q), b*std::sin(theta_b_q));
        cout <<"# Dipole amplitude, b=" << bv << endl;
        cout << "# r = (r,0)" << endl;
        const double MAXR = 5;
        const double MINR = 0.1;
        const int rpoints = 10;
        const double RSTEP = (MAXR-MINR)/(rpoints-1);
        double dipoles[rpoints];
        
#pragma omp parallel for
        for (int i=0; i<rpoints; i++)
        {
            double r = MINR + i*RSTEP;
            Vec rv(r,0);
            double d = integrator->DipoleAmplitudeBruteForce(diag, rv, bv);
            dipoles[i]=d;
        }
        for (int i=0; i<rpoints; i++)
        {
            double r = MINR + i*RSTEP;
            cout << r << " " << dipoles[i] << endl;
        }
        
        
    }
    else if (mode == DIPOLE_BRUTEFORCE_ANGLEDEP)
        {
            
            Vec bv(b,0);
            cout <<"# Dipole amplitude, b=" << b << endl;
            cout << "# r = " << r  << endl;
            const double MINTH = 0;
            const double MAXTH = M_PI/2.;
            const int THPOINTS = 3;
            const double THSTEP = (MAXTH-MINTH)/(THPOINTS-1);
            double *dipoles = new double[THPOINTS];
            cout <<"# th(r,b)   N(r,b,thrb)" << endl;
    #pragma omp parallel for
            for (int i=0; i<THPOINTS; i++)
            {
                double th = MINTH + i*THSTEP;
                Vec rv(r*std::cos(th),r*std::sin(th));
                double d = integrator->DipoleAmplitudeBruteForce(diag, rv, bv);
                dipoles[i]=d;
            }
            for (int i=0; i<THPOINTS; i++)
            {
                double th = MINTH + i*THSTEP;
                cout << th << " " << dipoles[i] << endl;
            }
            
            delete[] dipoles;
            
        }
    
    
    else if (mode == MIXED_SPACE_BRUTEFORCE_ANGLEDEP)
        {
            
            Vec bv(b,0);
            cout <<"#  Mixed space, b=" << b << ", q12=" << q12 << endl;
            
            const double MINTH = 0; //-M_PI;
            const double MAXTH = M_PI/2.0;
            const int THPOINTS = 2;
            const double THSTEP = (MAXTH-MINTH)/(THPOINTS-1);
            double *dipoles = new double[THPOINTS];
            cout <<"# th(r,b)   G2" << endl;
    #pragma omp parallel for
            for (int i=0; i<THPOINTS; i++)
            {
                double th = MINTH + i*THSTEP;
                Vec q12v(q12*std::cos(th),q12*std::sin(th));
                double d = integrator->MixedSpaceBruteForce(diag, q12v, bv);
                dipoles[i]=d;
            }
            for (int i=0; i<THPOINTS; i++)
            {
                double th = MINTH + i*THSTEP;
                cout << th << " " << dipoles[i] << endl;
            }
            
            delete[] dipoles;
            
        }
    
    else if (mode == MIXED_SPACE_BRUTEFORCE)
           {
               
               Vec q12v(q12,0);
               cout <<"#  Mixed space, q12=" << q12 << " angle " << theta_b_q << endl;
               
               const double MINB = 0;
               const double MAXB = 10;
               const int BPOINTS = 20;
               const double BSTEP = (MAXB-MINB)/(BPOINTS-1);
               double *dipoles = new double[BPOINTS];
               cout <<"# th(r,b)   G2" << endl;
       #pragma omp parallel for
               for (int i=0; i<BPOINTS; i++)
               {
                   double b = MINB + i*BSTEP;
                   Vec bv(b*std::cos(theta_b_q),b*std::sin(theta_b_q));
                   double d = integrator->MixedSpaceBruteForce(diag, q12v, bv);
                   dipoles[i]=d;
               }
               for (int i=0; i<BPOINTS; i++)
               {
                   double b = MINB + i*BSTEP;
                   cout << b << " " << dipoles[i] << endl;
               }
               
               delete[] dipoles;
               
           }
    
    
    
    delete integrator;
    return 0;
}
