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
   
    cout << "# Computing diagram " << diag << " (" << argv[1] << ")" << endl;
   
    
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
    
    
    integrator->SetMCIntPoints(mcintpoints);
    for (double q=0.05; q<4; q+=0.1)
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
