#include <iostream>
#include "functions.hpp"
#include "diagram_integrator.hpp"
using namespace std;
int main()
{
    DiagramIntegrator *integrator = new DiagramIntegrator;
   
    
    
    /*for (double q12=0.05; q12<3; q12+=0.05)
    {
        Vec q1(q12/2.,0);
        Vec q2(q12/2.,0);
        double diag2a = integrator->IntegrateDiagram(DIAG_2A, q1, q2);
        cout << q12 << " " << diag2a << endl;
    }
    */
    
   
    int mcintpoints = 1e7;
    integrator->SetMCIntPoints(mcintpoints);
    for (double q=0.05; q<10; q+=0.05)
    {
        Vec q1(0,0);
        Vec q2(q,0);
        
        integrator->UseInterpolator(false);
        double diag2a = integrator->IntegrateDiagram(DIAG_2A, q1, q2);
       
        cout << q << " " << diag2a << endl;;
        
    }
    
    delete integrator;
    return 0;
}
