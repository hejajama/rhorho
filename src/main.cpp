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
    
    
    
   
    int mcintpoints_anal = 1e7;
    int mcintpoints_num = 1e6;  // When F is computed numerically
    integrator->UseInterpolator(false);
    
    
    
    /*
    for (double q=0.05; q<4; q+=0.1)
    {
        Vec q1(q/2.,0);
        Vec q2(q/2.,0);
        cout << "### Diagram 2a" << endl;
        integrator->SetMCIntPoints(mcintpoints_anal);
        double diag2a = integrator->IntegrateDiagram(DIAG_2A, q1, q2);
        
        cout << "### Diagram 3a" << endl;
        integrator->SetMCIntPoints(mcintpoints_num);
        double diag3a_1 = integrator->IntegrateDiagram(DIAG_3A, q1, q2);
        // here in this test q1=q2
        double diag3a_2 = diag3a_1;
        //double diag3a_2 = integrator->IntegrateDiagram(DIAG_3A, q2, q1);
        
        cout << "### Diagram 3b" << endl;
        integrator->SetMCIntPoints(mcintpoints_anal);
        double diag3b_1 = integrator->IntegrateDiagram(DIAG_3B, q1, q2);
        double diag3b_2 = diag3b_1;
        //double diag3b_2 = integrator->IntegrateDiagram(DIAG_3B, q2, q1);
       
        
        cout << q << " " << diag2a << " " << diag3a_1 + diag3a_2 << " " << diag3b_1 + diag3b_2 << endl;;
        
    }
     */
    
    
    for (double q=0.05; q<4; q+=0.1)
    {
        Vec q1(q/2.,0);
        Vec q2(q/2.,0);
        cout << "### Diagram 5a" << endl;
        integrator->SetMCIntPoints(mcintpoints_num);
        double diag5a = integrator->IntegrateDiagram(DIAG_5A, q1, q2);
        
        cout << "### Diagram 5c" << endl;
        integrator->SetMCIntPoints(mcintpoints_num);
        double diag5a_1 = integrator->IntegrateDiagram(DIAG_5C, q1, q2);
        // here in this test q1=q2
        double diag5a_2 = diag5a_1;
        
       
        
        cout << q << " " << diag5a << " " << diag5a_1 + diag5a_2 << endl;;
        
    }
    
    
    /*
    for (int mcintpoints=1e6; mcintpoints<1e9; mcintpoints*=2)
    {
        Vec q1(-0.5/2.,0); Vec q2(0.5/2.,0);
        integrator->SetMCIntPoints(mcintpoints);
        double diag =integrator->IntegrateDiagram(DIAG_2A, q1, q2);
        cout << mcintpoints << " " << diag << endl;
    }
    */
    delete integrator;
    return 0;
}
