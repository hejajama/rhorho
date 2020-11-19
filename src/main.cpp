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
   
    Diagram diag = integrator->DiagramType(string(argv[1]));
    
  
    
    
    int mcintpoints = StrToReal(argv[2]);
   
    integrator->UseInterpolator(false);
    
    cout << "# Computing diagram " << diag << " (" << argv[1] << "). m = " << integrator->GetProton().GetM() <<" GeV, beta = " << integrator->GetProton().GetBeta() << " GeV" << endl;
    
    if ( diag == DIAG_3A or diag == DIAG_3A_2 or diag == DIAG_5A or diag == DIAG_5C or diag == DIAG_5C_1)
        mcintpoints /= 50; // there is one more intergal
                            // within the MC integral
    
    
    integrator->SetMCIntPoints(mcintpoints);
    
    
    // Test limit q1->0, q2->K, should vanish
    Vec K(1,0);
    for (double q=0.47; q<0.5; q+=0.002)
    {
        Vec q1(q,0);
        double d = integrator->IntegrateDiagram(diag, q1-K*0.5, q1*(-1)-K*0.5);
        cout << (q1-K*0.5).Len() << " " << d << endl;
    }
    
    //cout <<integrator->IntegrateDiagram(diag,Vec(1,0), Vec(0,0)) << endl;
    
    
    // q1=q2
    
    
    /*
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
   */
    
    // LO paper reproduce at finite q12
    /*
    double q12 = StrToReal(argv[3]);
    double theta_b_q= StrToReal(argv[4]);
    Vec q12v(q12*cos(theta_b_q), q12*sin(theta_b_q));
    
    cout << "# q12 = " << q12v << " orientation " << theta_b_q << endl;
    
    const double MAXK = 3;
    const int KPOINTS = 30;
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
     */
    
    delete integrator;
    return 0;
}
