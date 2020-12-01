#include <iostream>
#include <string>
#include <sstream>
#include <gsl/gsl_errno.h>

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
    WARD
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
    
    DiagramIntegrator *integrator = new DiagramIntegrator;
    long int mcintpoints = 1e6;
    string diagram = "LO";
    MODE mode = ONEDIM;
    
    double q12=0.5;
    double theta_b_q = 0;
    
    for (int i=1; i< argc; i++)
    {
        if (string(argv[i])=="-x")
            integrator->SetX(StrToReal(argv[i+1]));
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
        else if (string(argv[i])=="-ward")
            mode = WARD;
        else if (string(argv[i])=="-q12")
            q12 = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-theta_b_q")
            theta_b_q = StrToReal(argv[i+1]);
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unknown parameter " << argv[i] << endl;
            exit(1);
        }
    
    }
    
    integrator->GetProton().ComputeWFNormalizationCoefficient();
    
    
    
   
    Diagram diag = integrator->DiagramType(diagram);
    
   
    integrator->UseInterpolator(false);
    
    cout << integrator->InfoStr();
    
    if ( diag == DIAG_3A or diag == DIAG_3A_2 or diag == DIAG_5A or diag == DIAG_5C or diag == DIAG_5C_1)
        mcintpoints /= 70; // there is one more intergal
                            // within the MC integral
    
    cout << "# Diagram " << diagram << " id " << diag << " mcintpoints " << mcintpoints << endl << "#" << endl;
    integrator->SetMCIntPoints(mcintpoints);
    
    
    // Test limit q1->0, q2->K, should vanish
    if (mode == WARD)
    {
        Vec K(1,0);
        for (double q=0.47; q<0.5; q+=0.002)
        {
            Vec q1(q,0);
            double d = integrator->IntegrateDiagram(diag, q1-K*0.5, q1*(-1)-K*0.5);
            cout << (q1-K*0.5).Len() << " " << d << endl;
        }
    }
    
    if (mode == ONEDIM)
    {
    // q1=q2

        for (double q=0.02; q<5; q+=0.05)
        {
            Vec q1(q/2.,0);
            Vec q2(q/2.,0);
            
            double d = integrator->IntegrateDiagram(diag, q1, q2);
            if (integrator->Add_Q1Q2_exchange(diag))
            {
                cout << "#... adding cross graph q1<->q2" << endl;
                d +=integrator->IntegrateDiagram(diag, q1, q2);
            }
           
            
            cout << q << " " << d << endl;
            
        }
    }
    if (mode == TWODIM)
    {
    
    // LO paper reproduce at finite q12
    
        Vec q12v(q12*cos(theta_b_q), q12*sin(theta_b_q));
        
        cout << "# q12 = " << q12v << " orientation " << theta_b_q << endl;
        
        const double MAXK = 5;
        const int KPOINTS = 50;
        const double kstep =static_cast<double>(2*MAXK)/KPOINTS;
        
        for (double kx = -MAXK; kx <= MAXK + kstep/2.; kx += kstep  )
        {
            for (double ky = -MAXK; ky <= MAXK + kstep/2.; ky += kstep )
             {
                 Vec K(kx,ky);
                 Vec q1 = (q12v - K)*0.5;
                 Vec q2 = (q12v + K)*(-0.5);
                 
                 double d = integrator->IntegrateDiagram(diag, q1, q2);
                 if (integrator->Add_Q1Q2_exchange(diag))
                 {
                     cout << "#... adding cross graph q1<->q2" << endl;
                     d +=integrator->IntegrateDiagram(diag, q1, q2);
                 }
                 
                 cout << kx << " " << ky << " " << d << endl;
             }
        }
    }
    
    delete integrator;
    return 0;
}
