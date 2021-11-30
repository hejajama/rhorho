#include "functions.hpp"
#include "diagram_integrator.hpp"
//#define BOOST_TEST_MODULE Helpers
#include <boost/test/unit_test.hpp>

double inthelperf_mc_lo(double *vec, size_t dim, void* p);
double inthelperf_mc_diag2b(double *vec, size_t dim, void* p);

BOOST_AUTO_TEST_SUITE(HelperFunctions)
BOOST_AUTO_TEST_CASE(B0tests)
{
    BOOST_CHECK_CLOSE(B0(1, 2), 0.462098, 0.0001f);
    BOOST_CHECK_CLOSE(B0(10, 50), 0.0193589, 0.0001f);
    BOOST_CHECK_CLOSE(B0(0.01,0.001), 418.839640629419, 0.0001f);
}

BOOST_AUTO_TEST_CASE(Ftest)
{
    F_worker F(20,0.00001);
    
    Vec l(2,1);
    Vec l1(-1,3);
    double alpha=0.1; double m2=0.1*0.1;
    BOOST_CHECK_CLOSE(F.F_int_B0(l, l1, alpha, m2), -0.399215, 0.0001);
    
    l1=Vec(0,0,0);
    l=Vec(1,0);
    alpha=0.01;
    m2=0.05*0.05;
    BOOST_CHECK_CLOSE(F.F_int_B0(l, l1, alpha, m2), -1.03158, 0.0001);
}


BOOST_AUTO_TEST_CASE(LO_DIAG)
{
    double k1 = 0.1;
    double k2 = 0.2;
    double th1 = 0;
    double th2 = 0;
    
    double x1=0.3; double x2=0.4;
    double intvec[6] = {k1, th1, k2, th2, x1,x2}; // k1x  k1y k2x k2y x1 x2
    
    DiagramIntegrator *integrator = new DiagramIntegrator;
    integrator->GetProton().SetBeta(0.55);
    integrator->GetProton().SetM(0.26);
    integrator->GetProton().ComputeWFNormalizationCoefficient();
    
    Diagram diag = integrator->DiagramType("LO");
    //integrator->SetMCIntPoints(1e5);
    inthelper_diagint helper;
    helper.q1=Vec(1,0); helper.q2=Vec(0.5,0);
    helper.integrator=integrator;
    helper.diag = diag;
    
    double integrand =inthelperf_mc_lo(intvec, 6, &helper);
    
    // Jacobian and normalizatoin factors should be added to the mathematica result
    double mathematica =-1.0004e6;
    mathematica *= 16.*std::pow(M_PI,3.)*k1*k2  / (8.0*x1*x2*(1.-x1-x2)*std::pow(2.0*M_PI,6.0));
    
    BOOST_CHECK_CLOSE(integrand, mathematica, 0.1);
    
    
    
    delete integrator;
}



BOOST_AUTO_TEST_CASE(Q1_0_SET_2)
{
    // When q_1=0, we should get 6e_2 + 8h_2 + 8j_2=0 at all x_i, k_i
    double x1=0.3; double x2=0.3;
    double xg = 0.1;
    //[k1x,k1y,k2x,k2y,x1,x2,xg,kgx,kgy]
    if (xg > std::min(x1,1.-x2)) cerr << "Invalid x1,x2,xg in q1=0 test!" << endl;
    double intvec[9] = {0.1, -0.4, 0.2,0.3, x1,x2, xg, 0.2, -0.05};
    
    DiagramIntegrator *integrator = new DiagramIntegrator;
    integrator->GetProton().SetBeta(0.55);
    integrator->GetProton().SetM(0.26);
    integrator->GetProton().ComputeWFNormalizationCoefficient();
    
    Diagram diag1 = integrator->DiagramType("6e_2");
    Diagram diag2 = integrator->DiagramType("8h_2");
    Diagram diag3 = integrator->DiagramType("8j_2");
    
    inthelper_diagint helper;
    Vec q1(0,0); Vec q2(0.2,0.1);
    helper.q1=q1; helper.q2=q2;
    helper.integrator=integrator;
    
    helper.diag = diag1;
    double integrand1 = inthelperf_mc_diag2b(intvec, 9, &helper);
    
    helper.diag = diag2;
    double integrand2 = inthelperf_mc_diag2b(intvec, 9, &helper);
    
    helper.diag = diag3;
    double integrand3 = inthelperf_mc_diag2b(intvec, 9, &helper);
    
    BOOST_CHECK_SMALL(integrand1+integrand2+integrand3, 1e-6);
    
    // Check also integrated
    /*
    integrator->SetMCIntPoints(4e6);
    integrator->UseInterpolator(false);
    double int1 = integrator->IntegrateDiagram(diag1, q1, q2);
    double int2 = integrator->IntegrateDiagram(diag2, q1, q2);
    double int3 = integrator->IntegrateDiagram(diag3, q1, q2);
    
    
    BOOST_CHECK_SMALL(int1+int2+int3, 5.);
     */
    
    // Test ingerand in set 1 of finite diagrams, 3c' + 7h + 7j + 8h' + 8j'
    double sum=0;
    Diagram diags[5] = {integrator->DiagramType("3c_2"), integrator->DiagramType("7h"), integrator->DiagramType("7j"),
        integrator->DiagramType("8h_1"), integrator->DiagramType("8j_1") };
    double color_factors[5] = {3., 0.5, 0.5, 0.5, 0.5};
    for (int i=0; i<5; i++)
    {
        helper.diag = diags[i];
        double integrand = inthelperf_mc_diag2b(intvec, 9, &helper);
        sum += color_factors[i]*integrand;
    }
    BOOST_CHECK_SMALL(sum, 0.001);
    delete integrator;
}



BOOST_AUTO_TEST_SUITE_END()
