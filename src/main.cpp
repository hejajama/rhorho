#include <iostream>
#include <string>
#include <sstream>
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


int main(int argc, char* argv[])
{
    DiagramIntegrator *integrator = new DiagramIntegrator;
   
    Diagram diag = integrator->DiagramType(string(argv[1]));
    
    /*for (double q12=0.05; q12<3; q12+=0.05)
    {
        Vec q1(q12/2.,0);
        Vec q2(q12/2.,0);
        double diag2a = integrator->IntegrateDiagram(DIAG_2A, q1, q2);
        cout << q12 << " " << diag2a << endl;
    }
    */
    
    
    
    int mcintpoints = StrToReal(argv[2]);
   
    integrator->UseInterpolator(false);
    
    cout << "# Computing diagram " << diag << " (" << argv[1] << "). m = " << integrator->GetProton().GetM() <<" GeV, beta = " << integrator->GetProton().GetBeta() << " GeV" << endl;
    
    if ( diag == DIAG_3A or diag == DIAG_5A or diag == DIAG_5C)
        mcintpoints /= 50; // there is one more intergal
                            // within the MC integral
    
    integrator->SetMCIntPoints(mcintpoints);
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
    
    
    delete integrator;
    return 0;
}
